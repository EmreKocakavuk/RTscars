###########################################
################ Settings #################
library(DBI)
library(tidyverse)
library(ggpubr)
library(survival)
library(survminer)
library(gridExtra)

############  ######  ######  ######  ######  ############
# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq")) %>% select(aliquot_barcode, coverage_adj_mut_freq)
hm_freq <- mut_freq %>% filter(coverage_adj_mut_freq >=10, !grepl("-TP-", aliquot_barcode)) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  select(case_barcode) %>% distinct() %>% 
  filter(!case_barcode %in% c("GLSS-MD-0027", "GLSS-SF-0024", "GLSS-SU-0002", "GLSS-SU-0270", "TCGA-DU-6407")) #HM is only the 3rd surgery (exclude, since we focus only on surgery 1 and 2)

scars_set = dbReadTable(con,  Id(schema="analysis", table="scars_set")) %>% filter(!grepl("-DK-", aliquot_barcode)) %>% #filter(treatment_combi != "NA") %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"))
scars_pairs = dbReadTable(con,  Id(schema="analysis", table="scars_pairs")) %>% filter(!grepl("-DK-", case_barcode))
cases = dbReadTable(con,  Id(schema="clinical", table="cases")) %>% 
  filter(case_barcode %in% scars_set$case_barcode) %>% select(case_barcode, age = case_age_diagnosis_years)
scars_final = scars_pairs %>% inner_join(scars_set) %>% select(-aliquot_barcode, -surgery_number) %>% distinct() %>% 
  inner_join(cases)

surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
surgical_interval_2 <- surgeries %>% filter(surgery_number == 2, case_barcode %in% scars_final$case_barcode) %>% 
  select(case_barcode, surgical_interval_mo) %>% filter(!grepl("GLSS-DK-", case_barcode))

rt_dose <- surgeries %>% 
  group_by(case_barcode, surgery_number) %>% mutate(dose = sum(treatment_radiation_dose_gy, na.rm = T)) %>% 
  ungroup() %>% 
  select(case_barcode, dose, surgery_number) %>% 
  filter(dose != 0, surgery_number == 1) %>% distinct()

subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))

########################
## Load data
########################
mf_query <- dbGetQuery(con, read_file("sql/mf_snv_indel.sql"))

mf <- mf_query %>% 
  filter(variant_type == "DEL") %>%
  filter(is.na(reptype)) %>%
  filter(case_barcode %in% scars_final$case_barcode) %>% left_join(scars_set) %>% 
  mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-")) %>% left_join(surgical_interval_2) %>% 
  select(-aliquot_barcode, -surgery_number) %>% distinct() #join with scars set

cdkn2a_query <- dbGetQuery(con, read_file("sql/cell_cycle_status.sql"))
cdkn2a <- cdkn2a_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% filter(variant == "CDKN2A del")

cdkn2a_prim <- cdkn2a %>% select(case_barcode, aliquot_barcode = tumor_barcode_a, status) %>% 
  mutate(status = ifelse(status %in% c("P", "S"), "homdel", "wt"))
cdkn2a_rec <- cdkn2a %>% select(case_barcode, aliquot_barcode = tumor_barcode_b, status) %>% 
  mutate(status = ifelse(status %in% c("R", "S"), "homdel", "wt"))
cdkn2a_both <- rbind(cdkn2a_prim, cdkn2a_rec)

########################################################  
library(magrittr)

p_p16 <- cdkn2a_both %>% left_join(scars_final) %>% filter(idh_codel_subtype != "IDHwt") %>% #To clarify: CDKN2A is coding the tumor suppresor gene p16
  mutate(class = ifelse(aliquot_barcode%in%scars_final$tumor_barcode_a, "P", "R"))

p4 <- p_p16 %>% filter(class == "P") %$% table(variant=.$status) %>% as.data.frame() %>% mutate(class = "P")
p5 <- p_p16 %>% filter(class == "R") %$% table(variant=.$status, rt = .$received_rt) %>% as.data.frame() %>% mutate(class = "R")
p6 <- full_join(p4, p5) %>% unite("class_rt", class, rt) %>% spread(class_rt, Freq) %>% 
  rename("Primary" = "P_NA", "Recurrence\nRT+" = R_1, "Recurrence\nRT-" = R_0, "Variant" = variant) %>% 
  mutate(Variant = ifelse(Variant == "wt", "WT", "Homdel")) %>% 
  column_to_rownames("Variant")

# CDKN2A
plot_p16 <- 
  p6 %>% rownames_to_column("Variant") %>% 
  gather(key, value, -Variant) %>% spread(Variant, value) %>% mutate(prop = Homdel/(Homdel+WT)) %>% 
  ggplot(aes(x=key, y=prop, fill = key)) + 
  geom_bar(stat="identity") + 
  theme_bw() + theme(legend.position = "none") +
  labs(x="", y="Proportion\nwith homdel", title = "CDKN2A homdel", fill = "") +
  stat_compare_means(comparisons = list( c("Primary", "Recurrence\nRT+"), c("Recurrence\nRT-", "Recurrence\nRT+")),
                     label.y = c(.4, .35), 
                     label = "p.signif") + 
  geom_text(aes(label=paste(round(prop*100, 0), "%"), vjust=-.3)) + 
  scale_fill_brewer(palette = 6) + coord_cartesian(ylim=c(0,0.45))


fisher.test(p6[,-2])
fisher.test(p6[,-1])
fisher.test(p6)

pdf(file = "figures/main/Fig4_GLASS_cdkn2a_homdel.pdf", height = 5, width = 5, bg = "transparent", useDingbats = FALSE)
ggarrange(plot_p16, tableGrob(p6, theme = ttheme_minimal()),  align = "hv", ncol=1, heights = c(10,1))
dev.off()


dd <- mf %>% filter(idh_codel_subtype != "IDHwt") %>% 
  group_by(received_rt) %>% 
  mutate(del_cutoff = ifelse(mf_rs > quantile(mf_rs, 0.75), "Q4", 
                             ifelse(mf_rs > quantile(mf_rs, 0.5), "Q3", 
                                    ifelse(mf_rs > quantile(mf_rs, 0.25), "Q2", "Q1"))),
         quantile = ntile(mf_r, 3)) %>%
  ungroup()


surv_p16 <- 
  cdkn2a_rec %>% left_join(surgical_interval_2) %>% left_join(dd) %>% select(-received_rt) %>% left_join(scars_final) %>% 
  select(case_barcode, cdkn2a = status, idh_codel_subtype, received_rt, received_tmz, overall_survival, 
         vital_status, HM, age, interval = surgical_interval_mo,  del_cutoff, mf_r, mf_rs, quantile) %>% 
  mutate(vital_status = ifelse(vital_status == "dead", 2, 1), post_rec_surv = overall_survival-interval, 
         cdkn2a = factor(cdkn2a, levels = c("wt", "homdel")), 
         HM = factor(HM, levels = c("Non-HM", "HM"))) %>% 
  left_join(rt_dose)

km_del <- ggsurvplot(survfit(Surv(overall_survival, vital_status) 
                   ~ quantile, data = surv_p16 %>% filter(received_rt==1)), pval = T, risk.table = T, 
           pval.coord = c(200, 0.5),
           palette = c("royalblue3","#fc8d59","tomato3"),
           surv.median.line = "v", ylab = "Overall survival \n survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=17)", "Intermediate (n=16)", "High (n=16)"), legend.title = "Small deletion burden")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))


# Main figure separating overall survival (5a_left) into surgical interval (5b_middle) and post-recurrence survival (5b_right) 
fig5a_left <- ggsurvplot(survfit(Surv(overall_survival, vital_status) 
                   ~ quantile, data = surv_p16 %>% filter(received_rt==1)), pval = T, risk.table = T, 
           pval.coord = c(0, 0.2),
           palette = c("royalblue3","#fc8d59","tomato3"),
           surv.median.line = "v", ylab = "Overall survival \n survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=17)", "Intermediate (n=16)", "High (n=16)"), legend.title = "Newly acquired\nSmall deletion burden")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

fig5a_middle <- ggsurvplot(survfit(Surv(surgical_interval_mo, pfs_status) 
                   ~ quantile, data = surv_p16 %>% filter(received_rt==1) %>% mutate(pfs_status=1) %>% 
                     left_join(surgical_interval_2)), pval = T, risk.table = T, 
           pval.coord = c(150, 0.2),
           palette = c("royalblue3","#fc8d59","tomato3"),
           surv.median.line = "v", ylab = "Surgical interval \n survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=17)", "Intermediate (n=16)", "High (n=16)"), legend.title = "Newly acquired\nSmall deletion burden")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

fig5a_right <- ggsurvplot(survfit(Surv(post_rec_surv, vital_status) 
                   ~ quantile, data = surv_p16 %>% filter(received_rt==1)), pval = T, risk.table = T, 
           pval.coord = c(100, 0.2),
           palette = c("royalblue3","#fc8d59","tomato3"),
           surv.median.line = "v", ylab = "Post-recurrence \n survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=17)", "Intermediate (n=16)", "High (n=16)"), legend.title = "Newly acquired\nSmall deletion burden")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

pdf(file = "figures/main/Fig5A.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
ggarrange(fig5a_left, fig5a_middle, fig5a_right, nrow = 1, common.legend = T, legend = "top")
dev.off()

##
cox_del <- forestmodel::forest_model(coxph(Surv(overall_survival, vital_status) 
                                ~ mf_r +factor(received_tmz)  + idh_codel_subtype + age + cdkn2a,
                                data = surv_p16 %>% filter(received_rt==1)))

pdf(file = "figures/supp/FigS5C.pdf", height = 3, width = 9, bg = "transparent", useDingbats = FALSE)
cox_del
dev.off()

# Survival effects of small deletion burden are independent of RT dose (for rebuttal)
dose_reb1 <- 
surv_p16 %>% filter(received_rt == 1) %>% 
  ggplot(aes(x=dose)) +
  geom_histogram(binwidth=1, alpha=.5) + theme_bw() + labs(x="Dose (Gy)", y= "Sample count")

dose_reb2<- 
surv_p16 %>% filter(received_rt == 1) %>% 
  ggplot(aes(x=dose, y=mf_r)) +
  geom_point() + stat_cor() + geom_smooth(method = "lm") + theme_bw() + scale_y_log10() +
  labs(x="Dose (Gy)", y="Newly acquired\nsmall deletion burden\n(del/mb)") #+ylim(c(0,10))
  

dose_reb3 <- 
forestmodel::forest_model(coxph(Surv(overall_survival, vital_status) 
                                ~ log(mf_r) + age + dose,
                                data = surv_p16 %>% filter(received_rt==1)))

pdf(file = "figures/rebuttal/FigR7.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(dose_reb1, dose_reb2), dose_reb3, align = "hv", nrow=2)
dev.off()

# CDKN2A homdel analysis
km_cdkn2a <- ggsurvplot(survfit(Surv(overall_survival, vital_status) 
                   ~ cdkn2a , data = surv_p16),
           pval = T, risk.table = T, 
           pval.coord = c(150, 0.6),
           palette = c("#01665e", "#8c510a"), 
           surv.median.line = "v", ylab = "Overall \n survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("CDKN2A wt", "CDKN2A homdel"))$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_blank(), plot.title = element_text(face = "bold"))

cox_cdkn2a <- forestmodel::forest_model(coxph(Surv(overall_survival, vital_status) 
                                ~  cdkn2a +factor(received_rt) + factor(received_tmz) + idh_codel_subtype + age,
                                data = surv_p16), format_options = list(text_size = 3, shape = 15))


pdf(file = "figures/supp/FigS5a_GLASS_cdkn2a.pdf", height = 4, width = 13, bg = "transparent", useDingbats = FALSE)
ggpubr::ggarrange(km_cdkn2a, cox_cdkn2a, align = "hv",  widths = c(1,1.5))
dev.off()

### END ###