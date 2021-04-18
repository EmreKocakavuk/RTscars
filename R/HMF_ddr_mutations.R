## Effects of DNA damage response (DDR) gene mutations on RT-induced small deletion burden ##

library(tidyverse)
library(DBI)
library(ggpubr)
library(survminer)
library(survival)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

hmf_meta <- dbReadTable(con, Id(schema="public", table="hmf_meta")) %>% 
  mutate(status_hrd_msi = factor(status_hrd_msi, levels = c("HR_proficient", "HRD", "MSI")),
         RT = factor(RT, levels = c("RT-", "RT+ pal", "RT+ cur")))

res_query <- dbGetQuery(con, read_file(file = 'sql/mf_hartwig.sql'))
res <- res_query %>% select(sampleid, variant_type, nmut) %>% 
  inner_join(hmf_meta) %>% filter(variant_type%in% c("DEL", "INS", "SNP"))

mut_query <- dbGetQuery(con, read_file(file = 
"WITH 
anno_meta AS
(
SELECT ha.sampleid, variant_gene, sum(CASE WHEN variant_impact IN ('HIGH', 'MODERATE') THEN 1 ELSE 0 END) AS mut_n
FROM hwg_anno ha
INNER JOIN hmf_meta hm ON ha.sampleid = hm.sampleid
GROUP BY 1,2
), 
ddr_genes AS
( 
SELECT gene_symbol
FROM ref.ddr_gene_list
)
SELECT sampleid, gene_symbol, mut_n
FROM anno_meta am
INNER JOIN ddr_genes dr ON am.variant_gene = dr.gene_symbol"))


mut <- hmf_meta %>% left_join(mut_query) %>% 
  mutate(mut_n= ifelse(mut_n==0, 0, 1)) %>% 
  spread(gene_symbol, mut_n, , fill=0)

ID_sig_relative <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig")) %>% filter(sampleid %in% hmf_meta$sampleid)
ID_sig_absolute <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig_absolute")) %>% filter(sampleid %in% hmf_meta$sampleid)

#########################################
pdf(file = "figures/rebuttal/FigR8.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
mut %>% left_join(res) %>% gather(key, value, ATM, ATR, TP53, CHEK1, CHEK2, WEE1, PARP1, PRKDC) %>%
  filter(variant_type == "DEL") %>% 
  mutate(value = ifelse(value == 1, "mut", "WT"), radiotherapy_pretreatment = ifelse(radiotherapy_pretreatment == 0, "RT-", "RT+")) %>% 
  ggplot(aes(x=factor(value), y=nmut/2589)) + geom_boxplot() +
  EnvStats::stat_n_text() + 
  scale_y_log10() + facet_grid(~key) + theme_bw() + labs(x="", y="Small deletion burden\n(del/mb)") + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) +
  stat_compare_means(label.x = 1.1, label.y.npc = 0.63,label="p.format")
dev.off()  

pdf(file = "figures/rebuttal/FigR9.pdf", height = 6, width = 7, bg = "transparent", useDingbats = FALSE)
forestmodel::forest_model(lm(log(nmut/2589) ~
                                age +
                                primary_tumor_location +
                                RT +
                                status_hrd_msi +
                                factor(ATM) + 
                                factor(ATR) +
                                factor(CHEK1) +
                                factor(CHEK2) +
                                factor(PARP1) + 
                                factor(PRKDC) +
                                factor(TP53) +
                                factor(WEE1) +
                                factor(received_pyrimidine_antagonist) +
                                factor(received_taxane) +
                                factor(received_anthracyclines) +
                                factor(received_alkylating) +
                                factor(received_platinum) +
                                factor(received_hormonal_therapy) +
                                factor(received_targeted_therapy) +
                                factor(received_immunotherapy),
                              data = mut %>% left_join(res) %>% filter(variant_type == "DEL")), 
                          exponentiate = T) + labs(title = "Small deletion burden")
dev.off()

####### Check whether survival association is independent of these mutations
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

test_del_surv <- 
  res %>% mutate(nmut=nmut/2589) %>% 
  spread(variant_type, nmut) %>% 
  left_join(hmf_meta) %>% left_join(hmf_survival) %>% left_join(mut) %>% 
  filter(survival_mo >=0 ) %>% 
  mutate(survival_mo = survival_mo/30) %>% 
  mutate(status = ifelse(status == "dead", 2, 1)) %>% 
  filter(RT!= "RT-") %>% 
  filter(status_hrd_msi == "HR_proficient") %>% 
  filter(primary_tumor_location!="Others") %>% 
  mutate(quantile = ntile(DEL,3))

pdf(file = "figures/rebuttal/FigR10.pdf", height = 5, width = 9, bg = "transparent", useDingbats = FALSE)
forestmodel::forest_model(coxph(Surv(survival_mo, status) ~   
                                  age + log(DEL) + #primary_tumor_location + 
                                  #status_hrd_msi + 
                                  factor(ATM) + 
                                  factor(ATR) +
                                  factor(CHEK1) +
                                  factor(CHEK2) +
                                  factor(PARP1) +
                                  factor(PRKDC) +
                                  factor(TP53) +
                                  factor(WEE1),
                                  #factor(received_pyrimidine_antagonist) +
                                  #factor(received_taxane) +
                                  #factor(received_anthracyclines) +
                                  #factor(received_alkylating) +
                                  #factor(received_platinum) +
                                  #factor(received_hormonal_therapy) +
                                  #factor(received_targeted_therapy) +
                                  #factor(received_immunotherapy),
                                data = test_del_surv))
dev.off()

### END ###