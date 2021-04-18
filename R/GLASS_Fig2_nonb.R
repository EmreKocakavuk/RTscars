###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(parallel)
library(ggpubr)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq")) %>% select(aliquot_barcode, coverage_adj_mut_freq)
hm_freq <- mut_freq %>% filter(coverage_adj_mut_freq >=10, !grepl("-TP-", aliquot_barcode)) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  select(case_barcode) %>% distinct() %>% 
  filter(!case_barcode %in% c("GLSS-MD-0027", "GLSS-SF-0024", "GLSS-SU-0002", "GLSS-SU-0270", "TCGA-DU-6407")) #HM is only the 3rd surgery (exclude, since we focus only on surgery 1 and 2)

scars_set = dbReadTable(con,  Id(schema="analysis", table="scars_set")) %>% filter(!grepl("-DK-", aliquot_barcode)) %>% #filter(treatment_combi != "NA") %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"))
scars_pairs = dbReadTable(con,  Id(schema="analysis", table="scars_pairs")) %>% filter(!grepl("-DK-", case_barcode))
scars_final = scars_pairs %>% inner_join(scars_set) %>% select(-aliquot_barcode, -surgery_number) %>% distinct()

surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
surgical_interval_2 <- surgeries %>% filter(surgery_number == 2, case_barcode %in% scars_final$case_barcode) %>% 
  select(case_barcode, surgical_interval_mo) %>% filter(!grepl("GLSS-DK-", case_barcode))

subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))


################
################ Detailed analysis of genomic correalations with non-B DNA structures
################
nonb_random = dbReadTable(con,  Id(schema="analysis", table="random_variants_nonbdb"))

qres_nonb <- dbGetQuery(con, read_file("sql/del_nonbdb.sql"))


nonb <- qres_nonb %>% filter(is.na(reptype)) %>% 
  mutate(del_length = nchar(ref)-1) %>% filter(del_length <=20) %>% 
  mutate(chrom = ifelse(chrom == 23, "X",chrom)) %>% 
  inner_join(scars_final) %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"), 
         received_rt = ifelse(received_rt == 1, "RT+", "RT-"), 
         fraction_new = ifelse(fraction %in% c("P", "S"), "P/S", "R")) %>% 
  mutate(cutoff_0 = ifelse(dist == 0, 1, 0), 
         cutoff_100 = ifelse(dist <= 100, 1, 0)) %>% 
  filter(idh_codel_subtype != "IDHwt")  %>% 
  filter(is.na(reptype)) # Filter mutations in regions defined by repeatmasker

b_prim <- nonb %>% filter(fraction %in% c("P", "S")) %>% mutate(fraction = "Primary")
b_rec <- nonb %>% filter(fraction %in% c("R", "S")) %>% mutate(fraction = "Recurrence")
b_both <- full_join(b_prim, b_rec)

random_mean <- nonb_random %>% mutate(cutoff_0 = ifelse(dist == 0, 1, 0), 
                                      cutoff_100 = ifelse(dist <= 100, 1, 0)) %>% 
  group_by(feature_type) %>% mutate(mean = mean(dist), 
                                    median = median(dist),
                                    mean_cutoff_0 = mean(cutoff_0), 
                                    mean_cutoff_100 = mean(cutoff_100)) %>% ungroup() %>% 
  select(feature_type, mean_random_100 = mean_cutoff_100, mean_random_0 = mean_cutoff_0, mean_random = mean, median_random = median) %>% distinct()

mean_dist <- nonb %>% group_by(variant_type, case_barcode, fraction, feature_type) %>% mutate(mean = mean(dist), 
                                                                                                  median = median(dist),
                                                                                                  mean_cutoff_0 = mean(cutoff_0), 
                                                                                                  mean_cutoff_100 = mean(cutoff_100)) %>% ungroup() %>% 
  select(case_barcode,variant_type, fraction, feature_type, received_rt, HM, mean, mean_cutoff_0, mean_cutoff_100, received_tmz) %>% distinct()

b_dist <- b_both %>% group_by(variant_type, fraction,case_barcode, feature_type) %>% mutate(mean = mean(dist), 
                                                                                                  median = median(dist),
                                                                                                  mean_cutoff_0 = mean(cutoff_0), 
                                                                                                  mean_cutoff_100 = mean(cutoff_100)) %>% ungroup() %>% 
  select(case_barcode,variant_type, fraction, feature_type, received_rt, HM, mean, mean_cutoff_0, mean_cutoff_100, received_tmz) %>% distinct()

# Compare mean distances for each feature
fig_s2c <- 
  b_dist %>% filter(variant_type == "DEL") %>% left_join(random_mean) %>% 
  ggplot(aes(x=fraction, y=mean/1000)) + 
  geom_boxplot() + 
  theme_bw() + scale_y_log10() +
  stat_compare_means(paired=T, 
                     label.y.npc = 0.1,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + 
  facet_grid(received_rt + HM~feature_type, scales="free_y") + labs(y="Mean distance (kb)", x = "") + 
  theme(strip.background = element_rect(fill="white"))

fig2c_left <-
  mean_dist %>% filter(variant_type == "DEL", fraction == "R", HM == "Non-HM") %>%left_join(random_mean) %>% 
  ggplot(aes(x=reorder(feature_type, -mean_cutoff_100), y=-log10(mean_cutoff_100/mean_random_100), color = received_rt)) + 
  stat_summary(fun.data = mean_cl_normal, geom="pointrange", position = position_dodge(width=.5)) +
  theme_bw() +
  stat_compare_means(label.y=.5, aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_color_brewer(palette = 3, type = "qual") +
  geom_hline(yintercept = 0, linetype="dashed") + facet_wrap(~fraction) + 
  coord_flip(ylim=c(-1,1)) + labs(x="", y = "Background distribution", color = "") + 
  theme(legend.position = c(.9, 0.1), strip.background = element_rect(fill="white"))

# Plot distribution for Non-B DNA genomic features 
fig2c_right <- 
  b_both %>% filter(variant_type == "DEL", dist < 10000, HM == "Non-HM", received_rt=="RT+") %>% 
  ggplot(aes(x=dist/1000, color =feature_type, linetype=fraction)) +
  geom_line(aes(y=..y..), stat = "ecdf") +
  theme_bw() + 
  #scale_x_sqrt(expand=c(0,0), breaks = c(0,0.1, 0.5, 1,2,3,4,5,6,7,8,9,10)) +
  scale_x_sqrt(expand=c(0,0), breaks = c(0,1,2.5,5,7.5,10)) +
  scale_color_brewer(type = "qual", palette = 2) +
  labs(y="ECDF", x="Distance to feature (kb)", color = "", linetype = "") + coord_cartesian(xlim=c(0,10))

##

### Export figures
pdf(file = "figures/main/Fig2_GLASS_nonb.pdf", height = 3, width = 10, bg = "transparent", useDingbats = FALSE)
ggarrange(fig2c_left, fig2c_right, align = "hv", widths = c(1,1.5))
dev.off()

pdf(file = "figures/supp/FigS2_GLASS_nonb.pdf", height = 5, width = 10, bg = "transparent", useDingbats = FALSE)
fig_s2c
dev.off()

### END ###