###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(parallel)
library(ggpubr)
library(magrittr)

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

## Import queries
mf_query <- dbGetQuery(con, read_file("sql/mf_snv_indel.sql"))

mf <- mf_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(scars_set) %>% 
  mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-")) %>% left_join(surgical_interval_2) %>% 
  select(-aliquot_barcode, -surgery_number) %>% distinct() %>%  #join with scars set
  filter(is.na(reptype)) #filter out mutations in regions defined by repeatmasker

mf_sv_query <- dbGetQuery(con, read_file("sql/mf_sv.sql"))

mf_sv <- mf_sv_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(scars_set) %>% 
  mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-")) %>% left_join(surgical_interval_2) %>% 
  select(-aliquot_barcode, -surgery_number) %>% distinct() %>% 
  filter(!is.na(mf_r)) #join with scars set and remove samples without sufficient quality for SV calling
  

aneuploidy_query <- dbGetQuery(con, read_file("sql/cnv_missegregation.sql"))

aneuploidy <- aneuploidy_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(scars_set) %>% 
  mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-")) %>% left_join(surgical_interval_2) %>% 
  select(-aliquot_barcode, -surgery_number) %>% distinct() #join with scars set

cdkn2a_query <- dbReadTable(con,  Id(schema="analysis", table="gatk_cnv_by_gene")) %>% filter(gene_symbol == "CDKN2A") %>% 
  mutate(call = ifelse(hlvl_call == -2, "homdel", "wt")) %>% 
  select(aliquot_barcode, call) 

cdkn2a <- scars_final %>% rename(aliquot_barcode = tumor_barcode_a) %>% left_join(cdkn2a_query) %>% 
  rename(call_a = call, tumor_barcode_a= aliquot_barcode, aliquot_barcode = tumor_barcode_b) %>% left_join(cdkn2a_query) %>% 
  rename(call_b = call, tumor_barcode_b=aliquot_barcode)

#################
#Small indels and SNVs
#################
# Plot deletion burden across RT categories

f1a <- 
  mf %>% filter(variant_type == "DEL") %>%  
  ggplot(aes(x=received_rt, y=mf_r)) + 
  geom_boxplot(aes(fill=received_rt), width = 0.75) + 
  theme_classic(base_size = 18) + #
  stat_compare_means(label.x = 1.4, label.y = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_y_log10() +
  labs(y = "Newly acquired deletions\n(deletions/mb)", x = "", fill = "", alpha = "") +
  scale_fill_manual(values = alpha(c("white", "black"), .7)) + theme(legend.position = "none") 

f1a

mf %>% filter(variant_type == "DEL") %>% 
  select(case_barcode, received_rt, mf_is, mf_rs, mf_r, idh_codel_subtype, HM) %>% 
  mutate(fc = mf_rs/mf_is) %>% 
  mutate(diff=mf_rs - mf_is) %>% 
  mutate(double = ifelse(fc >=2, "yes", "no")) %>% 
  mutate(increase = ifelse(diff >=1, "yes", "no")) %>% 
  mutate(mf_r_cut = ifelse(mf_r >=1, "yes", "no")) %$%
  table(.$received_rt, .$double)
  #table(.$received_rt, .$double, .$HM)
  #table(.$received_rt, .$mf_r_cut, .$HM)
  #fisher.test(table(.$received_rt, .$double))
  #fisher.test(table(.$received_rt, .$mf_r_cut)) ### newly acquired del burden > 1/mb 


fig_s1a <-
mf %>% 
  filter(variant_type %in% c("DEL", "INS", "SNP")) %>% 
  ggplot(aes(x=received_rt, y=mf_r)) + 
  geom_boxplot(aes(fill=received_rt), width = 0.75) + 
  #geom_jitter(aes(alpha=0.5), width = 0.1)+
  theme_classic(base_size=18) +
  stat_compare_means(aes(label = sprintf("P = %2.1e", as.numeric(..p.format..))),
                     label.x = 1.2, label.y.npc = 0.8) +
  scale_y_log10() +
  labs(y = "Newly acquired mutations\n(mutations/mb)", x = "", fill = "", alpha = "") +
  scale_fill_manual(values = alpha(c("white", "black"), .7)) + theme(legend.position = "none") + 
  facet_grid(~variant_type, scales = "free_y") + 
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white")) 


fig_s1b <- mf %>% filter(variant_type == "DEL") %>%  mutate(subtype = ifelse(idh_codel_subtype == "IDHwt", "IDHwt", "IDHmut")) %>% 
  ggplot(aes(x=received_rt, y=mf_r)) + 
  geom_boxplot(aes(fill=received_rt), width = 0.75) + 
  theme_classic(base_size = 18) + #
  stat_compare_means(label.x = 1.4, 
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_y_log10() +
  labs(y = "Newly acquired deletions\n(deletions/mb)", x = "", fill = "", alpha = "") +
  scale_fill_manual(values = alpha(c("white", "black"), .7)) + theme(legend.position = "none") + facet_wrap(~subtype) + 
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white") ) 


f1b <- 
  mf %>% filter(variant_type == "DEL") %>% 
  gather(key, value, mf_is, mf_rs) %>% 
  mutate(key = ifelse(key == "mf_is", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin((aes(fill=key))) +
  geom_line(aes(group=case_barcode), alpha=0.3) +
  geom_boxplot(width=0.1) +
  facet_grid(received_rt~ HM, scales = "free_y") +
  scale_y_log10() +
  theme_classic(base_size = 18) + #(base_size=20
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9, paired = T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  labs(y = "Small deletion\nburden (del/mb)", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

f1b

t_mod <-  mf %>% filter(!is.na(received_tmz)) %>% 
  select(case_barcode, variant_type, mf_r, idh_codel_subtype, received_rt, received_tmz, surgical_interval_mo, HM) %>% 
  spread(variant_type, mf_r, fill = NA) %>% 
  mutate(TMZ = factor(received_tmz), Hypermutation = factor(HM, levels =rev( c("HM", "Non-HM"))), RT = factor(received_rt), `Subtype` = idh_codel_subtype, `Surgical Interval` = surgical_interval_mo, 
         MF = SNP + DEL + INS, 
         INDEL = DEL + INS)
  
f1c <- forestmodel::forest_model(lm(formula = log(DEL+1) ~ TMZ + Hypermutation + RT + Subtype  + `Surgical Interval`, data=t_mod), 
                                 exponentiate = T) #limits = c(-.13,.18)

f1c

reb_t_mod <- mf %>%
  mutate(project = substr(tumor_barcode_a, 1, 4)) %>% 
  mutate(center = ifelse(project == "GLSS", substr(case_barcode, 6, 7), "TCGA"), seq = substr(tumor_barcode_a, 21, 23)) %>% 
  mutate(center = ifelse(center %in% c("JX", "HK", "MD"), "JX", center)) %>% 
  select(center, case_barcode, variant_type, mf_r, mf_rs, idh_codel_subtype, received_rt, received_tmz, surgical_interval_mo, HM, seq) %>% 
  spread(variant_type, mf_r, fill = NA) %>% 
  mutate(TMZ = factor(received_tmz), Hypermutation = factor(HM, levels =rev( c("HM", "Non-HM"))), RT = factor(received_rt), `Subtype` = idh_codel_subtype, `Surgical Interval` = surgical_interval_mo, 
         MF = SNP + DEL + INS, Center = center, Seq = seq,
         INDEL = DEL + INS) #Effects independent of batch effects

rf1c <- forestmodel::forest_model(lm(formula = log(DEL+1) ~ RT + Center + Hypermutation + Subtype + `Surgical Interval` ,   #Subtype + TMZ + Center + RT +   Hypermutation  + `Surgical Interval`, 
                                     data=reb_t_mod), #%>% filter(!center %in% c("DF")), 
                                 exponentiate = T) #limits = c(-.13,.18)

fig_s1d <- 
  jtools::plot_summs(lm(formula = log(DEL + 1) ~  TMZ + Hypermutation + RT + Subtype + `Surgical Interval`, data=t_mod),
                             lm(formula = log(INS + 1) ~  TMZ  + Hypermutation + RT + Subtype + `Surgical Interval`, data=t_mod),
                             lm(formula = log(INDEL + 1) ~  TMZ + Hypermutation + RT + Subtype + `Surgical Interval`, data=t_mod),
                             lm(formula = log(SNP + 1) ~  TMZ + Hypermutation + RT + Subtype + `Surgical Interval`, data=t_mod),
                             lm(formula = log(MF + 1) ~  TMZ + Hypermutation + RT + Subtype + `Surgical Interval`, data=t_mod),
                             model.names = c("Deletion burden", "Insertion burden", "Indel burden","SNV burden", "Tumor mutational burden (TMB)"), 
                             plot.distributions = F, 
                             confint =T, exp = F) + theme_bw(base_size = 20) + labs(y="", x="Estimate") #+ coord_cartesian(xlim=c(0.9, 5.5))#+ scale_x_sqrt(breaks= c(1,2,4,6)) #+ coord_cartesian(xlim=c(0.75,6))
#################
#SVs
#################
# sanity check
mf_sv %>% mutate(si = ifelse(s + i ==0, NA, s+i)) %>% 
  group_by(idh_codel_subtype,received_rt, svtype) %>% 
  mutate(mean = mean(r/si, na.rm = T), 
         median=median(r/si, na.rm=T)) %>% 
  select(received_rt, svtype, mean, median) %>% distinct()
  
### Distribution fold change
supp_inc50_distr <- 
mf_sv %>% 
  mutate(fc = (r)/(i+s)) %>% filter(fc <= 2) %>% 
  ggplot(aes(x=fc)) + geom_histogram() + facet_grid(svtype~received_rt) + theme_bw() + 
  geom_vline(xintercept = 0.5, color="tomato") + labs(x = "Percent increase\n(Post- vs. pretreatment)", y = "Count") + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))


###Increase 50% !
inc50 <- mf_sv %>% mutate(increase = ifelse((r+s) > 1.5*(i+s), "Yes", "No"))

#Figure 4a
main_sv_inc50 <- 
inc50 %>%  filter(idh_codel_subtype != "IDHwt") %>% select(-n) %>%  
  add_count(received_rt, svtype, increase) %>% 
  select(received_rt, svtype, increase, n) %>% distinct() %>% #spread(increase, n) %>% 
  mutate(svtype = factor(svtype, levels=c("BND", "DUP", "DEL", "INV"), 
                         labels = c("TRA", "DUP", "DEL", "INV"))) %>% 
  ggplot(aes(x=received_rt, y=n, fill = increase)) + 
  geom_bar(position = "fill", stat = "identity") + #, stat="identity"
  geom_text(aes(label=paste("n =",n)), position = position_fill(vjust = 0.5), color = "grey20") +
  facet_grid(~svtype) + theme_classic2() + 
  labs(fill = "Increase > 50%?",  y = "Proportion of samples", x = "") + 
  scale_fill_manual(values = c("grey95", "#377eb8"))



# Apply two-sided fisher's exact test
inc50 %>% filter(svtype == "DEL", idh_codel_subtype != "IDHwt")%$% #Excluding IDHwt due to large imbalance betweeen RT+ and RT-
  table(.$increase, .$received_rt) %>% fisher.test()

# DEL: P=0.03206
# INV: P=0.02068
# DUP: P=0.7651
# BND: P=1

# Multivariate logistic regression model
t_log_del <- inc50 %>% filter(svtype%in%c("DEL"), !is.na(received_tmz)) %>%
  mutate(increase = ifelse(increase == "Yes", 1, 0), increase = as.factor(increase))

t_log_inv <- inc50 %>% filter(svtype%in%c("INV"), !is.na(received_tmz)) %>%
  mutate(increase = ifelse(increase == "Yes", 1, 0), increase = as.factor(increase))

t_log_dup <- inc50 %>% filter(svtype%in%c("DUP"), !is.na(received_tmz)) %>%
  mutate(increase = ifelse(increase == "Yes", 1, 0), increase = as.factor(increase))

t_log_bnd <- inc50 %>% filter(svtype%in%c("BND"), !is.na(received_tmz)) %>%
  mutate(increase = ifelse(increase == "Yes", 1, 0), increase = as.factor(increase))


supp_log_del <- 
forestmodel::forest_model(glm(increase ~ received_rt +  received_tmz +  surgical_interval_mo + idh_codel_subtype, 
                              family=binomial(link='logit'),data = t_log_del), exponentiate = T) + labs(title = "DEL")

supp_log_inv <- 
  forestmodel::forest_model(glm(increase ~ received_rt +  received_tmz +  surgical_interval_mo + idh_codel_subtype, 
                                family=binomial(link='logit'),data = t_log_inv), exponentiate = T)+ labs(title = "INV")

supp_log_dup <- 
  forestmodel::forest_model(glm(increase ~ received_rt +  received_tmz +  surgical_interval_mo + idh_codel_subtype, 
                                family=binomial(link='logit'),data = t_log_dup), exponentiate = T)+ labs(title = "DUP")

supp_log_bnd <- 
  forestmodel::forest_model(glm(increase ~ received_rt +  received_tmz +  surgical_interval_mo + idh_codel_subtype, 
                                family=binomial(link='logit'),data = t_log_bnd), exponentiate = T)+ labs(title = "BND")



#################
#Aneuploidy
#################
#important here to do paired analysis
aneuploidy_paired <- 
aneuploidy %>% 
  mutate(seq = substr(tumor_barcode_a, 21, 23)) %>% 
  filter(idh_codel_subtype != "IDHwt") %>% 
  mutate(loss_b = n_simple_missegregation_loss_b + n_complex_missegregation_loss_neut_b + n_complex_missegregation_gain_loss_b, 
         loss_a = n_simple_missegregation_loss_a + n_complex_missegregation_loss_neut_a + n_complex_missegregation_gain_loss_a)


main_aneuploidy_whole_chr <-
  aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_a, n_simple_missegregation_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Whole chromosome\naneuploidy", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_distr <- 
  aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_a, n_simple_missegregation_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_a", "Primary", "Recurrence")) %>% 
    ggplot(aes(x=value, fill = key)) +
  geom_density(alpha=0.75, aes(y=..count..)) + 
  geom_bar(position = position_dodge2(preserve = "single",padding = 0), width = 1, alpha=0.75) +
    facet_grid(~received_rt) + theme_bw() +
  labs(x="Whole chromosome loss", y = "Count", fill = "") + 
  theme(legend.position = c(0.9, 0.8)) + 
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(strip.background =element_rect(fill="white"))

main_aneuploidy_whole_chr_del <-
  aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_loss_a, n_simple_missegregation_loss_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_loss_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Whole chromosome\nloss", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_distr_del <- 
  aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_loss_a, n_simple_missegregation_loss_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_loss_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=value, fill = key)) +
  geom_density(alpha=0.5) + facet_grid(~received_rt) + theme_bw() +
  labs(x="Whole chromosome loss", y = "Density", fill = "") + 
  theme(legend.position = c(0.9, 0.8)) + 
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(strip.background =element_rect(fill="white"))


main_aneuploidy_whole_chr_gain <-
  aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_gain_a, n_simple_missegregation_gain_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_gain_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
   facet_grid(~received_rt , scales="free")  + 
  labs(y = "Whole chromosome\ngain", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_distr_gain <- 
aneuploidy_paired %>%  gather(key, value, n_simple_missegregation_loss_a, n_simple_missegregation_loss_b) %>% 
  mutate(key = ifelse(key == "n_simple_missegregation_loss_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=value, fill = key)) +
  geom_density(alpha=0.5) + facet_grid(~received_rt) + theme_bw() +
  labs(x="Whole chromosome loss", y = "Density", fill = "") + 
  theme(legend.position = c(0.9, 0.8)) + 
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(strip.background =element_rect(fill="white"))

supp_aneuploidy_partial <-
  aneuploidy_paired %>%  gather(key, value, n_complex_missegregation_a, n_complex_missegregation_b) %>% 
  mutate(key = ifelse(key == "n_complex_missegregation_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Partial\naneuploidy", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_partial_gain_neut <-
  aneuploidy_paired %>%  gather(key, value, n_complex_missegregation_gain_neut_a, n_complex_missegregation_gain_neut_b) %>% 
  mutate(key = ifelse(key == "n_complex_missegregation_gain_neut_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Chromosome arm\ngain/neutral", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_partial_loss_neut <-
  aneuploidy_paired %>%  gather(key, value, n_complex_missegregation_loss_neut_a, n_complex_missegregation_loss_neut_b) %>% 
  mutate(key = ifelse(key == "n_complex_missegregation_loss_neut_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Chromosome arm\nloss/neutral", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_partial_gain_loss <-
  aneuploidy_paired %>%  gather(key, value, n_complex_missegregation_gain_loss_a, n_complex_missegregation_gain_loss_b) %>% 
  mutate(key = ifelse(key == "n_complex_missegregation_gain_loss_a", "Primary", "Recurrence")) %>% 
  ggplot(aes(x=key, y=value)) + 
  geom_violin(aes(fill = key)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(key, received_rt))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  geom_line(aes(group=case_barcode), linetype="solid", alpha = 0.2) + theme_bw() + 
  facet_grid(~received_rt , scales="free")  + 
  labs(y = "Chromosome arm\ngain/loss", x = "", fill = "", alpha = "") +
  scale_fill_brewer(palette = 6, type = "qual", direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

##########
########### Check aneuploidy association with RT and CDKN2A
# Poisson model
a_model <- 
  aneuploidy_paired %>% 
  mutate(diff = n_simple_missegregation_loss_b - n_simple_missegregation_loss_a, 
         fc = n_simple_missegregation_loss_b/n_simple_missegregation_loss_a)  %>% 
  select(-received_rt) %>% 
  left_join(cdkn2a) %>% mutate(call_b = factor(call_b, levels = c("wt", "homdel")))


pdf(file = "figures/supp/Fig4_GLASS_aneuploidy_poisson.pdf", height = 4, width = 8.5, bg = "transparent", useDingbats = FALSE)
forestmodel::forest_model(glm(formula = n_simple_missegregation_loss_b ~  idh_codel_subtype + factor(received_rt)*call_b + factor(received_tmz) + surgical_interval_mo, 
                             data=a_model, family = "poisson"), 
                          exponentiate = T)
dev.off()

############## Check characteristics of small del high/low at primary
del_hilo <- mf %>% 
  filter(variant_type == "DEL") %>% filter(received_rt=="RT+") %>% 
  #group_by(idh_codel_subtype) %>% 
  mutate(tile = ntile(mf_r, 2)) %>% 
  ungroup() %>% 
  select(case_barcode, tile) %>% distinct() 

reb_prim_1 <- 
mf %>% left_join(del_hilo) %>% mutate(tile = ifelse(tile == 1, "low", "high")) %>% 
  mutate(variant_type = ifelse(variant_type == "SNP", "SNV", variant_type)) %>% 
  filter(variant_type %in% c("DEL", "INS", "SNV"), !is.na(tile)) %>% 
  ggplot(aes(x=factor(tile), y = mf_is)) + geom_boxplot(aes(fill = tile)) + 
  theme_bw() + facet_grid(HM~variant_type, scales ="free") + 
  stat_compare_means(label = "p.format", label.x = 1.4, label.y.npc = 0.9) + 
  labs(y = "Mutational burden\nprimary tumor (mut/mb)", x = "Post-treatment\nsmall deletion burden", fill = "", alpha = "") +
  scale_fill_brewer(palette = 9, direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

reb_prim_2 <- 
aneuploidy %>% left_join(del_hilo) %>% mutate(tile = ifelse(tile == 1, "low", "high")) %>% 
  filter(!is.na(tile)) %>% 
  ggplot(aes(x=factor(tile), y = n_simple_missegregation_a)) +
  geom_boxplot(aes(fill = tile)) + 
  theme_bw() + facet_grid(~HM, scales ="free") +
  stat_compare_means(label = "p.format", label.x = 1.4, label.y.npc = 0.9) + 
  labs(y = "Aneuploidy score\nprimary tumor", x = "Post-treatment\nsmall deletion burden", fill = "", alpha = "") +
  scale_fill_brewer(palette = 9, direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

ggarrange(reb_prim_1, reb_prim_2)

######### Export figures
pdf(file = "figures/main/Fig1A+B+C_GLASS.pdf", height = 8, width = 10, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(f1a, f1b, widths = c(1.5,3), legend = "none", align = "hv"), 
          ggarrange(f1c), align = "hv", nrow = 2, heights = c(1.4,1))
dev.off()

pdf(file = "figures/supp/FigS1A_GLASS.pdf", height = 3.5, width = 7, bg = "transparent", useDingbats = FALSE)
fig_s1a
dev.off()

pdf(file = "figures/supp/FigS1B_GLASS.pdf", height = 4, width = 5, bg = "transparent", useDingbats = FALSE)
fig_s1b
dev.off()

pdf(file = "figures/supp/FigS1D_GLASS.pdf", height = 6, width = 12, bg = "transparent", useDingbats = FALSE)
fig_s1d
dev.off()

pdf(file = "figures/rebuttal/FigR4.pdf", height = 5, width = 9, bg = "transparent", useDingbats = FALSE)
rf1c
dev.off()

pdf(file = "figures/rebuttal/FigR5.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
ggarrange(reb_prim_1, reb_prim_2)
dev.off()



pdf(file = "figures/supp/FigS4_GLASS_sv_inc50_distr.pdf", height = 4, width = 5, bg = "transparent", useDingbats = FALSE)
supp_inc50_distr
dev.off()

pdf(file = "figures/main/Fig4_GLASS_SV.pdf", height = 5, width = 8, bg = "transparent", useDingbats = FALSE)
main_sv_inc50
dev.off()

pdf(file = "figures/supp/FigS4_GLASS_SV_logistic.pdf", height = 6, width = 16, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(supp_log_del, supp_log_inv), ggarrange(supp_log_dup, supp_log_bnd), ncol=1)
dev.off()

pdf(file = "figures/main/Fig4_GLASS_aneuploidy.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
main_aneuploidy_whole_chr
dev.off()

pdf(file = "figures/main/Fig4_GLASS_aneuploidy_loss.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
main_aneuploidy_whole_chr_del
dev.off()

pdf(file = "figures/main/Fig4_GLASS_aneuploidy_gain.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
main_aneuploidy_whole_chr_gain
dev.off()

pdf(file = "figures/supp/Fig4_GLASS_aneuploidy_partial.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
supp_aneuploidy_partial
dev.off()

pdf(file = "figures/supp/Fig4_GLASS_aneuploidy_partial_gain.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
supp_aneuploidy_partial_gain_neut
dev.off()

pdf(file = "figures/supp/Fig4_GLASS_aneuploidy_partial_loss.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
supp_aneuploidy_partial_loss_neut
dev.off()

pdf(file = "figures/supp/Fig4_GLASS_aneuploidy_partial_gain_loss.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
supp_aneuploidy_partial_gain_loss
dev.off()


pdf(file = "figures/supp/FigS4_GLASS_aneuploidy_del_distr.pdf", height = 3, width = 6, bg = "transparent", useDingbats = FALSE)
supp_aneuploidy_distr_del
dev.off()


### END ###