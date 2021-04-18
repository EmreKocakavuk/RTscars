#### Majority of analyses for HMF (Hartwig Medical Foundation) dataset ####
###############################################################################

library(tidyverse)
library(DBI)
library(pracma)
library(ggpubr)
library(survminer)
library(survival)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

res_query <- dbGetQuery(con, read_file(file = 'sql/mf_hartwig.sql'))
hmf_meta <- dbReadTable(con, Id(schema="public", table="hmf_meta")) %>% 
  mutate(status_hrd_msi = factor(status_hrd_msi, levels = c("HR_proficient", "HRD", "MSI")),
         RT = factor(RT, levels = c("RT-", "RT+ pal", "RT+ cur")))

res <- res_query %>% left_join(hmf_meta)

### Load in sequence context for HMF
qres_delseq_query <- dbGetQuery(con, read_file("sql/del_sequence_context_hw.sql"))

qres_delseq <- qres_delseq_query %>% left_join(hmf_meta)

# Load in ID signatures
ID_sig_relative <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig")) %>% filter(sampleid %in% hmf_meta$sampleid)
ID_sig_absolute <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig_absolute")) %>% filter(sampleid %in% hmf_meta$sampleid)

## Define functions for statistical testing
testWilcoxGroup <- function(df) {
  wtest = wilcox.test(df$value ~ df$RT, paired = F, conf.int = TRUE)
  data.frame(n = nrow(df),
             median_a = median(df$value[df$RT == 1], na.rm = T),
             median_b = median(df$value[df$RT == 0], na.rm = T),
             mean_a = mean(df$value[df$RT == 1], na.rm = T),
             mean_b = mean(df$value[df$RT == 0], na.rm = T),
             statistic = wtest$statistic,
             estimate = wtest$estimate,
             or = exp(wtest$estimate), 
             lcl = wtest$conf.int[1],
             ucl = wtest$conf.int[2],
             wilcox_p = wtest$p.value,
             #fdr = p.adjust(wtest$p.value, method = "fdr"), 
             test_str = sprintf("n=%s\n%s",
                                nrow(df),
                                case_when(wtest$p.value < 0.0001 ~ "P<0.0001",
                                          wtest$p.value > 0.05 ~ sprintf("P=%s", format(round(wtest$p.value, 2),scientific=F)),
                                          TRUE ~ sprintf("P=%s", format(round(wtest$p.value, 4),scientific=F)))),
             stringsAsFactors = FALSE)
}



########### Load data from db

#### CDKN2A by gene
library(magrittr)
cnv_by_gene_query = dbGetQuery(con, read_file("sql/hw_cnv_by_gene.sql"))

#########################################################################################
### Aneuploidy per chromosome arm
hwav_query <- dbGetQuery(con, "SELECT * FROM hwav")

#########################################################################################
### Detailed aneuploidy query
aneuploidy_query <- dbGetQuery(con, read_file("sql/cnv_missegregation_hw.sql"))

#########################################################################################
# Load survival info 
hmf_survival <- dbGetQuery(con, 
                           "SELECT (DATE_PART('year', death_date) - DATE_PART('year', biopsy_date)) * 12 *30 +
       (DATE_PART('month', death_date) - DATE_PART('month', biopsy_date)) *30 +
	   (DATE_PART('day', death_date) - DATE_PART('day', biopsy_date))AS survival_time_dead,
	   (DATE_PART('year', treatment_end_date) - DATE_PART('year', biopsy_date)) * 12 * 30 +
       (DATE_PART('month', treatment_end_date) - DATE_PART('month', biopsy_date)) *30 + 
	   (DATE_PART('day', treatment_end_date) - DATE_PART('day', biopsy_date)) AS survival_time_alive, 
	   sampleid, death_date, biopsy_date, treatment_end_date, 
extract(year from age(biopsy_date,to_date(birth_year::varchar, 'YYYY')))::integer AS age
FROM hwm") %>% 
  mutate(status = ifelse(is.na(death_date), "alive", "dead"), 
         survival_mo = ifelse(status == "dead", survival_time_dead, survival_time_alive))

###################
## Some exploratory numbers on the treatment cohorts
hmf_meta %$%
  table(RT)

hmf_meta %$%
  table(RT, status_hrd_msi)


########################################
######### Small deletion burden analysis 
variant_mb <- res

variant_mb %>% mutate(mutfreq = nmut/2589) %>% 
  group_by(variant_type, RT) %>% 
  summarise(mean(mutfreq), median(mutfreq))

variant_mb %>% mutate(mutfreq = nmut/2589) %>% filter(variant_type == "DEL") %>% 
  group_by(primary_tumor_location, variant_type, RT) %>% 
  summarise(mean(mutfreq), median(mutfreq))

own_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", 
                      "#B09C85FF", "#ff7f00", "#6a3d9a")

# Small deletion burden separated by HRD/MSI and RT 
variant_status <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  mutate(RT = ifelse(RT=="RT-","RT-", "RT+")) %>% 
  ggplot(aes(x=status_hrd_msi, y = nmut/2589, fill = RT)) + 
  theme_bw() + scale_y_log10() +
  geom_boxplot() + 
  scale_fill_brewer() + labs(x="", y="Small deletion\nburden (del/mb)") +
  theme(legend.position = "none") + EnvStats::stat_n_text() + 
  stat_compare_means(label.x = 1.25, label.y.npc = 0.9, 
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..))))

variant_mb_fig1 <- 
variant_mb %>% filter(variant_type == "DEL") %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589)) + 
  facet_grid(~primary_tumor_location, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(fill = primary_tumor_location, alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="Small deletion\nburden (del/mb)") +
  scale_fill_manual(values = own_colors) +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.63,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  #stat_compare_means(comparisons = list(c("RT-", "RT+ pal"), c("RT-", "RT+ cur"), c("RT+ cur", "RT+ pal")))  +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

rebuttal_variant_mb_fig1 <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589)) + 
  facet_grid(~primary_tumor_location, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(fill = primary_tumor_location, alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="Small deletion\nburden (del/mb)") +
  scale_fill_manual(values = own_colors) +
  stat_compare_means(label.x = 1.2, label.y.npc = 0.85,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  EnvStats::stat_n_text(angle = 90, hjust=-.1)
  

supp_variant_lung <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  #filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Lung") %>% filter(primary_tumor_subtype %in% c("Non-Small Cell", "Small Cell")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title= "Lung") +
  scale_fill_manual(values = "#8491B4FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.85,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()
  
rebuttal_variant_lung <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Lung") %>% filter(primary_tumor_subtype %in% c("Non-Small Cell", "Small Cell")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title= "Lung") +
  scale_fill_manual(values = "#8491B4FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.85,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

supp_variant_breast <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  #filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Breast") %>% filter(!primary_tumor_subtype %in% c("ER-positive/HER2 unknown", "Subtype unknown")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Breast") +
  scale_fill_manual(values = "#4DBBD5FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

rebuttal_variant_breast <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Breast") %>% filter(!primary_tumor_subtype %in% c("ER-positive/HER2 unknown", "Subtype unknown")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Breast") +
  scale_fill_manual(values = "#4DBBD5FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

supp_variant_sarcoma <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  #filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Bone/Soft tissue") %>% filter(primary_tumor_subtype %in% c("Ewing sarcoma", "Sarcoma", "Liposarcoma")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Bone/Soft tissue") +
  scale_fill_manual(values = "#E64B35FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

rebuttal_variant_sarcoma <- 
variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Bone/Soft tissue") %>% filter(primary_tumor_subtype %in% c("Ewing sarcoma", "Sarcoma", "Liposarcoma")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Bone/Soft tissue") +
  scale_fill_manual(values = "#E64B35FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

rebuttal_variant_skin <- 
variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Skin") %>% filter(primary_tumor_subtype %in% c("Melanoma", "Skin squamous cell carcinoma")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Skin") +
  scale_fill_manual(values = "#B09C85FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

rebuttal_variant_uterus <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Uterus") %>%  
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  facet_grid(~primary_tumor_subtype, scales = "free_y") +
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="", title = "Uterus") +
  scale_fill_manual(values = "#6a3d9a") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + EnvStats::stat_n_text()

rebuttal_variant_glioma <- 
  variant_mb %>% filter(variant_type == "DEL") %>% 
  mutate(RT = ifelse(RT=="RT-", "RT-", "RT+")) %>% 
  #filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location == "Nervous system") %>%  
  filter(!primary_tumor_subtype %in% c("Meningeoma", "Solitair fibreuze tumor in brain", "Medulloblastoma")) %>% 
  ggplot(aes(x=factor(RT), y = nmut/2589, fill = primary_tumor_location)) + 
  theme_bw() + scale_y_log10() +
  geom_boxplot(aes(alpha = factor(RT))) + #, outlier.shape = NA) + 
  labs(x="", y="Small deletion burden\n(del/mb)", title = "Glioma") +
  scale_fill_manual(values = "#91D1C2FF") +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_alpha_manual(values = c(0.1,0.5, 1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #coord_cartesian(ylim = c(0.01,11)) +
  theme(axis.title.x=element_blank()) + EnvStats::stat_n_text()

variant_mb_fig2 <- 
variant_mb %>% filter(variant_type == "DEL") %>%  
  left_join(hmf_meta) %>% 
  select(sampleid, primary_tumor_location, RT) %>% distinct() %>% 
  add_count(primary_tumor_location,RT) %>% 
  select(-sampleid) %>% distinct() %>%  
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=RT, y = n)) + 
  geom_bar(aes(fill = primary_tumor_location, alpha =RT), stat = "identity", width = 0.7) + 
  geom_text(aes(label=paste("n =",n)), vjust=0.2,hjust=-.1, size=4, angle=90) +
  coord_cartesian(ylim=c(0,650)) +
  facet_grid(~primary_tumor_location) +
  theme_bw() +
  #theme_bw(base_size =20) + 
  labs(x="", y="Number of samples") +
  #ggsci::scale_fill_lancet() +
  scale_fill_manual(values = own_colors) +
  scale_alpha_manual(values = c(0.1, 0.5,1)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

pdf(file = "figures/main/Fig1D.pdf", height = 8, width = 14, bg = "transparent", useDingbats = FALSE)
ggarrange(variant_mb_fig1, variant_mb_fig2, ncol = 1, align = "hv", heights = c(1,0.8))
dev.off()

pdf(file = "figures/supp/FigS1F.pdf", height = 5, width = 10, bg = "transparent", useDingbats = FALSE)
ggarrange(supp_variant_breast, 
          ggarrange(supp_variant_lung, supp_variant_sarcoma), nrow=2) %>% 
  annotate_figure(left = text_grob("Small deletion burden\n(del/mb)", rot = 90))
dev.off()

pdf(file = "figures/supp/FigS1G.pdf", height = 5, width = 5, bg = "transparent", useDingbats = FALSE)
variant_status
dev.off()

pdf(file = "figures/rebuttal/FigR1.pdf", height = 5, width = 13, bg = "transparent", useDingbats = FALSE)
rebuttal_variant_mb_fig1
dev.off()

pdf(file = "figures/rebuttal/FigR2.pdf", height = 8, width = 15, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(rebuttal_variant_lung, rebuttal_variant_breast, widths = c(2,4), align = "hv"),
          ggarrange(rebuttal_variant_sarcoma, rebuttal_variant_skin, rebuttal_variant_uterus, widths = c(3,2,2), nrow=1, align = "hv"), 
          nrow=2, align="hv") %>%  
  annotate_figure(left = text_grob("Small deletion burden\n(del/mb)", rot = 90))
dev.off()

pdf(file = "figures/rebuttal/FigR3.pdf", height = 4, width = 4, bg = "transparent", useDingbats = FALSE)
rebuttal_variant_glioma
dev.off()


variant_mb_supp_del <- 
forestmodel::forest_model(glm(log(nmut/2589) ~
      age +
      primary_tumor_location +
      RT +
      status_hrd_msi + 
      factor(received_pyrimidine_antagonist) +
      factor(received_taxane) +
      factor(received_anthracyclines) +
      factor(received_alkylating) +
      factor(received_platinum) +
      factor(received_hormonal_therapy) +
      factor(received_targeted_therapy) +
      factor(received_immunotherapy),
    data = variant_mb %>% filter(variant_type == "DEL") %>% left_join(hmf_meta)), 
    exponentiate = T) + labs(title = "Small deletion burden")

variant_mb_supp_ins <- 
forestmodel::forest_model(glm(log(nmut/2589) ~
                                age +
                                primary_tumor_location +
                                RT +
                                status_hrd_msi + 
                                factor(received_pyrimidine_antagonist) +
                                factor(received_taxane) +
                                factor(received_anthracyclines) +
                                factor(received_alkylating) +
                                factor(received_platinum) +
                                factor(received_hormonal_therapy) +
                                factor(received_targeted_therapy) +
                                factor(received_immunotherapy),
                              data = variant_mb %>% filter(variant_type == "INS")%>%left_join(hmf_meta)), 
                          exponentiate = T) +
  labs(title = "Small insertion burden")

variant_mb_supp_snp <- 
forestmodel::forest_model(glm(log(nmut/2589) ~
                                age +
                                primary_tumor_location +
                                RT +
                                status_hrd_msi + 
                                factor(received_pyrimidine_antagonist) +
                                factor(received_taxane) +
                                factor(received_anthracyclines) +
                                factor(received_alkylating) +
                                factor(received_platinum) +
                                factor(received_hormonal_therapy) +
                                factor(received_targeted_therapy) +
                                factor(received_immunotherapy),
                              data = variant_mb %>% filter(variant_type == "SNP") %>% left_join(hmf_meta)), 
                          exponentiate = T) + 
  labs(title = "SNV burden")

pdf(file = "figures/supp/FigS1H.pdf", height = 20, width = 10, bg = "transparent", useDingbats = FALSE)
ggarrange(variant_mb_supp_del, variant_mb_supp_ins,variant_mb_supp_snp, ncol = 1, align = "hv")
dev.off()

################################################
######### Small deletion length analysis 
del_length <- 
  qres_delseq %>% mutate(length = nchar(seq_del)) %>% 
  filter(length <=20) %>% 
  group_by(sampleid) %>% 
  mutate(mean_length = mean(length)) %>% 
  ungroup() %>% 
  select(sampleid, mean_length) %>% distinct()

hmf_del_length <- 
  hmf_meta %>% 
  left_join(del_length) 

fig_del_length <- 
hmf_del_length %>% 
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=RT, y =mean_length, fill = RT)) + 
  geom_boxplot() +  theme_bw() +
  #scale_y_log10() + 
  labs(x="", y="Mean deletion length (bp)", alpha = "") +
  scale_fill_manual(values = alpha(c("white","#636363", "black"), .7)) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) +
  stat_compare_means(label.x = 1.4, label="p.format")

fig_del_length_supp <- 
  hmf_del_length %>% 
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=RT, y =mean_length,  alpha = RT)) + 
  geom_boxplot(aes(fill = primary_tumor_location)) +  theme_bw() +
  labs(x="", y="Mean deletion length (bp)", alpha = "") +
  stat_compare_means(label.x = 1.4, 
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  scale_fill_manual(values = own_colors) +
  facet_grid(~primary_tumor_location) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) 


hmf_del_length_distr <- 
qres_delseq %>%  mutate(deletion_length = nchar(seq_del)) %>% 
  filter(deletion_length <=20) %>% select(sampleid, deletion_length) %>% 
  left_join(hmf_meta) %>% filter(!is.na(RT))
  
fig_length_distr <- 
hmf_del_length_distr %>% filter(deletion_length >1) %>% 
  add_count(sampleid, deletion_length, name = "n") %>% 
  add_count(sampleid, name = "total") %>% 
  mutate(prop = n/total) %>% 
  select(sampleid, deletion_length, RT, n, total, prop,primary_tumor_location) %>% 
  distinct() %>% 
  ggplot(aes(x=factor(deletion_length), y=prop, color = factor(RT))) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", position = position_dodge2(width=0.5)) +
  stat_compare_means(label="p.signif", hide.ns = T) + 
  scale_color_brewer(palette = 1, type = "seq", direction = 1) +
  scale_y_log10(breaks=c(0.01, 0.03, 0.1, 0.3)) + 
  labs(y = "Proportion", x = "Deletion length (bp)", color = "") +
  theme_bw()

pdf(file = "figures/main/Fig2B.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
ggarrange(fig_del_length, fig_length_distr, align = "hv", widths = c(1,2))
dev.off()

pdf(file = "figures/supp/FigS2B.pdf", height = 5, width = 13, bg = "transparent", useDingbats = FALSE)
fig_del_length_supp
dev.off()

####### Small deletion sequence context analysis
hmf_detailed <- 
  qres_delseq %>% 
  left_join(hmf_meta) %>% 
  mutate(deletion_length = nchar(seq_del)) %>% 
  filter(deletion_length <=20) %>%  
  mutate(category_length = case_when(deletion_length == 1 ~ "1bp",
                                     deletion_length == 2 ~ "2bp", 
                                     deletion_length == 3 ~ "3bp", 
                                     deletion_length == 4 ~ "4bp", 
                                     deletion_length >= 5 ~ "+5bp",
                                     TRUE ~ "NA")) %>% 
  mutate(category_simple = case_when(deletion_length == 1 ~"1bp", 
                                     deletion_length > 1 & homology_length > 1 ~ ">1bp MH", 
                                     deletion_length > 1 & homology_length <= 1 ~ ">1bp Non-MH"))

pancan_homology_panel <- 
hmf_detailed %>% 
  #filter(status_hrd_msi == "HR_proficient") %>% 
  select(sampleid, category_simple, RT, primary_tumor_location) %>% 
  add_count(sampleid, name = "t") %>% 
  add_count(sampleid, category_simple, name = "n") %>% 
  mutate(prop = n/t) %>% 
  select (sampleid, category_simple, prop, RT, primary_tumor_location) %>% distinct() %>% 
  spread(category_simple, prop, fill= 0) %>% 
  gather(category_simple, prop, `1bp`:`>1bp Non-MH`) %>% 
  mutate(category_simple = factor(category_simple, levels = c(">1bp MH", "1bp", ">1bp Non-MH"))) %>% 
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=factor(RT), y = prop)) +
  geom_boxplot(aes(color = category_simple)) + 
  scale_color_manual(values = c("1bp" = "#5d5f60", ">1bp Non-MH" = "#2d004b", ">1bp MH" = "#f6711f")) +
  facet_grid(~category_simple) +
  facet_grid(~category_simple~primary_tumor_location) +
  stat_compare_means(label.x = 1.4, label.y.npc = 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"), legend.position = "top") + labs(x="", y = "Proportion", color = "Deletion\ncharacteristic")
  

pdf(file = "figures/supp/FigS2F.pdf", height = 5, width = 13, bg = "transparent", useDingbats = FALSE)
pancan_homology_panel
dev.off()



####### Analysis of association between RT and chromosome-arm aneuploidy
hmf_hwav <- 
  hwav_query%>% select(-variant_type, -abs_arm_count, -rel_arm_freq) %>% distinct() %>% 
  group_by(sampleid) %>% mutate(ploidy= mean(arm_ploidy)) %>% ungroup() %>% 
  select(-arm_call, - arm, -arm_ploidy, - covered_arm_size) %>% distinct()

# Poisson model
poisson_form <- hmf_meta %>% left_join(hmf_survival) %>% 
  left_join(aneuploidy_query) %>% 
  left_join(cnv_by_gene_query) %>% filter(gene_symbol == "CDKN2A") %>% 
  left_join(hmf_hwav) %>% 
  #filter(ploidy <= 2.5) %>% 
  mutate(call = ifelse(ploidy_based_cn == -1, "homdel", "wt"),
         call = factor(call, levels= rev(c("homdel", "wt"))))

poisson_figure <- 
  forestmodel::forest_model(glm(formula = n_simple_missegregation_loss ~ RT + call + 
                                  primary_tumor_location, 
                                data = poisson_form, family = "poisson"), exp=T)

aneuploidy_hartwig_1 <- 
  hmf_meta %>% 
  left_join(aneuploidy_query) %>% 
  left_join(cnv_by_gene_query) %>% filter(gene_symbol == "CDKN2A") %>% 
  mutate(call = ifelse(ploidy_based_cn == -1, "homdel", "wt")) %>% 
  ggplot(aes(x=call, y = n_simple_missegregation_loss)) + theme_bw() + 
  geom_violin(aes(fill = call)) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(call, RT))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.4, 
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  labs(x="CDKN2A", y = "Whole chromosome loss") + 
  theme(legend.position = "none") + 
  scale_fill_brewer(type = "seq", palette = 7, direction = -1) + facet_grid(~RT) +
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_distr_1 <- 
  hmf_meta %>% 
  left_join(aneuploidy_query) %>% 
  left_join(cnv_by_gene_query) %>% filter(gene_symbol == "CDKN2A") %>% 
  mutate(call = ifelse(ploidy_based_cn == -1, "homdel", "wt")) %>% 
  ggplot(aes(x=n_simple_missegregation_loss, fill=call)) + theme_bw() + 
  geom_density(alpha=0.5) + facet_grid(~RT) + 
  #stat_compare_means(label.x = 1.4, 
  #                   aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  labs(x="Whole chromosome loss", y = "Density", fill = "CDKN2A") + 
  theme(legend.position = c(0.9, 0.8)) + 
  scale_fill_brewer(type = "seq", palette = 7, direction = -1) + facet_grid(~RT) +
  theme(strip.background =element_rect(fill="white"))

aneuploidy_hartwig_2 <- 
  hmf_meta %>% left_join(hmf_survival) %>% 
  left_join(aneuploidy_query) %>% 
  left_join(cnv_by_gene_query) %>% filter(gene_symbol == "CDKN2A") %>% 
  mutate(call = ifelse(ploidy_based_cn == -1, "homdel", "wt")) %>% 
  mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
  ggplot(aes(x=RT, y = n_simple_missegregation_loss)) + theme_bw() + 
  geom_violin(aes(fill = call)) + 
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.05) + 
  geom_count(alpha=0.5,aes(size = after_stat(prop), group = interaction(call, RT))) + scale_size_area(max_size = 7) +
  stat_compare_means(label.x = 1.4, 
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  labs(x="", y = "Whole chromosome loss") + 
  theme(legend.position = "none") + 
  scale_fill_brewer(type = "seq", palette = 7, direction = -1) + facet_grid(~call) +
  theme(legend.position = "none", strip.background =element_rect(fill="white"))

supp_aneuploidy_distr_2 <- 
  hmf_meta %>% 
  left_join(aneuploidy_query) %>% 
  left_join(cnv_by_gene_query) %>% filter(gene_symbol == "CDKN2A") %>% 
  mutate(call = ifelse(ploidy_based_cn == -1, "CDKN2A homdel", "CDKN2A WT")) %>% 
  ggplot(aes(x=n_simple_missegregation_loss, fill=RT)) + theme_bw() + 
  geom_density(alpha=0.5) +
  #stat_compare_means(label.x = 1.4, 
  #                   aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
  labs(x="Whole chromosome loss", y = "Density", fill = "RT") + 
  theme(legend.position = c(0.9, 0.8)) + 
  scale_fill_brewer(type = "seq", palette = 6, direction = 1) + facet_grid(~call) +
  theme(strip.background =element_rect(fill="white"))

pdf(file = "figures/main/Fig4_pancan_aneuploidy_poisson.pdf", height = 9, width = 11, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(aneuploidy_hartwig_1, aneuploidy_hartwig_2, align="hv"), 
          poisson_figure, align = "v", nrow=2, heights = c(1.25,1))
dev.off()

pdf(file = "figures/supp/Fig4_aneuploidy_distr.pdf", height = 5, width = 6, bg = "transparent", useDingbats = FALSE)
ggarrange(supp_aneuploidy_distr_1, 
          supp_aneuploidy_distr_2, nrow=2)
dev.off()

#########
#########
#########
######### Survival

test_hmf <- 
  ID_sig_absolute %>% left_join(hmf_meta) %>% left_join(hmf_survival) %>% 
  left_join(poisson_form) %>%  
  filter(survival_mo >=0) %>% 
  mutate(status = ifelse(status == "dead", 2, 1)) %>% 
  filter(RT!= "RT-") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location!="Others") %>% 
  mutate(quantile = ntile(ID8, 3)) %>% 
  mutate(chr_loss = ntile(n_simple_missegregation_loss,3))

### ID8 KM-plot
km_rtsig <- 
ggsurvplot(survfit(Surv(survival_mo/30, status)
                   ~ quantile , data = test_hmf), pval = T, risk.table = T, 
           pval.coord = c(20, 0.5),
           palette = c("#4575b4", "#fc8d59", "#d73027"), 
           surv.median.line = "v", ylab = "Survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=320)", "Intermediate (n=319)", "High (n=319)"), legend.title = "ID8")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

### CDKN2A KM-Plot
km_cdkn2a <- 
ggsurvplot(survfit(Surv(survival_mo/30, status)
                   ~ call , data = test_hmf), pval = T, risk.table = T, 
           pval.coord = c(20, 0.5),
           palette = c("#01665e", "#8c510a"), 
           surv.median.line = "v", ylab = "Survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("CDKN2A wt (n=829)", "CDKN2A homdel (n=129)"))$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_blank(), plot.title = element_text(face = "bold"))

### Aneuploidy KM-Plot
km_aneuploidy <- 
ggsurvplot(survfit(Surv(survival_mo/30, status)
                   ~ chr_loss , data = test_hmf), pval = T, risk.table = T, 
           pval.coord = c(20, 0.5),
           palette = c("#4575b4", "#fc8d59", "#d73027"), 
           surv.median.line = "v", ylab = "Survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=320)", "Intermediate (n=319)", "High (n=319)"), legend.title = "Aneuploidy")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

#### Export survival figures
pdf(file = "figures/supp/FigS5_pancan_survival.pdf", height = 4, width = 14, bg = "transparent", useDingbats = FALSE)
ggarrange(km_cdkn2a, km_aneuploidy, km_rtsig, align = "hv", nrow = 1)
dev.off()

########## Survival analysis with small deletion burden
test_del_surv <- 
  variant_mb %>%
  spread(variant_type, nmut) %>% 
  mutate(del_indel = DEL/(DEL+INS), del_mut = DEL/(DEL+INS+SNP), INDEL = DEL+INS) %>% 
  left_join(hmf_meta) %>% left_join(hmf_survival) %>% 
  filter(survival_mo >=0 ) %>% 
  mutate(survival_mo = survival_mo/30) %>% 
  mutate(status = ifelse(status == "dead", 2, 1)) %>% 
  filter(RT!= "RT-") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location!="Others") %>% 
  mutate(quantile = ntile(DEL,3))

pdf(file = "figures/main/Fig5_pancan_survival_del_burden.pdf", height = 4, width = 5, bg = "transparent", useDingbats = FALSE)
ggsurvplot(survfit(Surv(survival_mo, status)
                   ~ quantile , data = test_del_surv), pval = T, risk.table = T, 
           pval.coord = c(20, 0.5),
           palette = c("#4575b4", "#fc8d59", "#d73027"), 
           surv.median.line = "v", ylab = "Survival probability", xlab = "Time (months)", size=2,  legend=c(0.7, 0.8), 
           legend.labs = c("Low (n=320)", "Intermediate (n=319)", "High (n=319)"), legend.title = "Small deletion burden")$plot +
  theme(legend.text = element_text(size = 14, color = "black"), legend.title = element_text(size = 14, color="black", face = "bold"))

dev.off()

### END ###