# Clonality Analysis in GLASS using the CCF (cancer cell fraction) metric
# Visualized in Extended Data Figure 1c

###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(ggpubr)
library(EnvStats)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq")) %>% select(aliquot_barcode, coverage_adj_mut_freq)
hm_freq <- mut_freq %>% filter(coverage_adj_mut_freq >=10, !grepl("-TP-", aliquot_barcode)) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  select(case_barcode) %>% distinct() %>% 
  filter(!case_barcode %in% c("GLSS-MD-0027", "GLSS-SF-0024", "GLSS-SU-0002", "GLSS-SU-0270", "TCGA-DU-6407")) #HM is only the 3rd surgery (exclude, since we focus only on surgery 1 and 2)

scars_set = dbReadTable(con,  Id(schema="analysis", table="scars_set")) %>% filter(!grepl("-DK-", aliquot_barcode)) %>% #filter(treatment_combi != "NA") %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"))
scars_pairs = dbReadTable(con,  Id(schema="analysis", table="scars_pairs")) %>% filter(!grepl("-DK-", case_barcode))
scars_final = scars_pairs %>% inner_join(scars_set) %>% select(-aliquot_barcode, -surgery_number) %>% 
  distinct() %>% filter(idh_codel_subtype != "IDHwt")

ccf_query <- dbGetQuery(con, read_file(file="sql/ccf_del.sql"))

#Compare mean ccf per patient

pdf(file = "figures/supp/SuppFig1c.pdf", height = 5, width = 5, bg = "transparent", useDingbats = FALSE)
ccf_query %>% inner_join(scars_final) %>%
  group_by(case_barcode, fraction) %>% 
  mutate(ccf_mean = mean(ccf), RT= ifelse(received_rt==1, "RT+", "RT-"), 
         fraction = factor(fraction, levels = c("P", "S", "R"))) %>% 
  ungroup() %>% 
  select(case_barcode, fraction, RT, ccf_mean, 
         received_tmz, idh_codel_subtype, HM) %>% 
  distinct() %>% 
  ggplot(aes(x=RT, y=ccf_mean)) + geom_boxplot(aes(fill=RT), width = 0.75) + 
  stat_compare_means(label.x = 1.2, label.y=0.1,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + 
  theme_bw() + facet_grid(HM~fraction) + 
  labs(y = "Mean CCF\n(Cancer Cell Fraction)", x = "", fill = "", alpha = "") +
  scale_fill_manual(values = alpha(c("white", "black"), .7)) + 
  theme(legend.position = "none", strip.background = element_rect(fill="white"))
dev.off()


### END ###