library(dndscv)
library(tidyverse)
library(DBI)
library(ggthemes)
library(gridExtra)
library(ggplot2)

data("cancergenes_cgc81", package="dndscv")

###################
## DB connection
###################

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

res_fraction_query <- dbGetQuery(con, read_file("sql/dndscv_input_by_fraction.sql"))

res_fraction <- res_fraction_query %>% inner_join(scars_final) %>% mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-"))


###################
## Run dNdS for all samples in cohort
###################

dnds_all <- res_fraction %>%
  select(case_barcode,chrom,pos,ref,mut) %>% 
  distinct() %>%
  dndscv(refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)

message("Ran dNdSCV for ", n_distinct(res_fraction$case_barcode) - length(dnds_all$exclsamples))

dnds_all_global = dnds_all$globaldnds
dnds_all_sel_cv = dnds_all$sel_cv
dnds_all_gene_ci = geneci(dnds_all, gene_list = known_cancergenes)

###################
## Examine results
###################
print(dnds_all$globaldnds)
print(dnds_all$nbreg$theta)

sel_cv = dnds_all$sel_cv
print(head(sel_cv, sum(sel_cv$qglobal_cv < 0.05),n=25), digits = 3)

###################
## Run dNdS seperately for private/shared variants and RT-treatment
###################

result_list <- lapply(na.omit(unique(res_fraction$received_rt)), function(st) {
  result_list <- lapply(na.omit(unique(res_fraction$fraction)), function(fr) {
    message("Computing dNdS for ", fr, " and ", st)
    qres_subset = res_fraction %>% filter(fraction == fr, received_rt == st) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    message(".. Found ", nrow(dnds_subset$annotmuts), " mutations")
    globaldnds <- cbind(fraction = fr, RT = st, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, received_rt = st, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, received_rt = st, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list,'[[',3))
  list(res1,res2,res3)
})
dnds_fraction_global <- data.table::rbindlist(lapply(result_list,'[[',1))
dnds_fraction_sel_cv <- data.table::rbindlist(lapply(result_list,'[[',2)) %>% 
  mutate(qind_cv = p.adjust(pind_cv, method = "BH"))
dnds_fraction_gene_ci <- data.table::rbindlist(lapply(result_list,'[[',3))
print(dnds_fraction_global)
print(dnds_fraction_sel_cv)

dnds_fraction_sel_cv %>% filter(fraction == "R", qind_cv < 0.05)

###################
## Run dNdS seperately for private/shared variants and RT-treatment (ONLY FOR IDHmut)
###################
res_fraction_idhmut <- res_fraction %>% filter(idh_codel_subtype != "IDHwt")


  result_list_idhmut <- lapply(na.omit(unique(res_fraction_idhmut$fraction)), function(fr) {
    result_list_idhmut <- lapply(na.omit(unique(res_fraction_idhmut$received_rt)), function(rt) {
    message("Computing dNdS for ", fr, "and", rt)
    qres_subset = res_fraction_idhmut %>% filter(fraction == fr, received_rt == rt ) %>% select(case_barcode,chrom,pos,ref,mut) %>% distinct()
    dnds_subset = dndscv(qres_subset, refdb = "hg19", outmats = TRUE, max_coding_muts_per_sample = NULL)
    message(".. Ran dNdSCV for ", n_distinct(qres_subset$case_barcode) - length(dnds_subset$exclsamples))
    message(".. Found ", nrow(dnds_subset$annotmuts), " mutations")
    globaldnds <- cbind(fraction = fr, RT = rt, dnds_subset$globaldnds)
    sel_cv <- cbind(fraction = fr, received_rt = rt, dnds_subset$sel_cv)
    gene_ci <- cbind(fraction = fr, received_rt = rt, geneci(dnds_subset, gene_list = known_cancergenes))
    list(globaldnds, sel_cv, gene_ci)
  })
  res1 <- data.table::rbindlist(lapply(result_list_idhmut,'[[',1))
  res2 <- data.table::rbindlist(lapply(result_list_idhmut,'[[',2))
  res3 <- data.table::rbindlist(lapply(result_list_idhmut,'[[',3))
  list(res1,res2,res3)
})

dnds_fraction_idhmut_global <- data.table::rbindlist(lapply(result_list_idhmut,'[[',1))
dnds_fraction_idhmut_sel_cv <- data.table::rbindlist(lapply(result_list_idhmut,'[[',2)) %>% 
  mutate(qind_cv = p.adjust(pind_cv, method = "BH"))
dnds_fraction_idhmut_sgene_ci <- data.table::rbindlist(lapply(result_list_idhmut,'[[',3))
print(dnds_fraction_idhmut_global)
print(dnds_fraction_idhmut_sel_cv)

dnds_fraction_idhmut_sel_cv %>% filter(fraction == "R", qind_cv < 0.05)

###################
## Plot dNdS for all data
###################

p <- dnds_all_global %>%
  ggplot(aes(x = name, y = mle, ymin = cilow, ymax = cihigh, color = name)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 18) +
  labs(y = "Global dN/dS", color = "Variant Type")
plot(p)


###################
## Plot dNdS-CV by RT - by fraction
###################

p<-dnds_fraction_global %>%
  filter(name == 'wall') %>%
  mutate(received_rt = factor(RT, levels = c("RT+", "RT-")),
         fraction = factor(fraction, levels = c("P", "S", "R")),
         fraction = fct_recode(fraction, "Primary only" = "P", "Shared" = "S", "Recurrence only" = "R")) %>%
  ggplot(aes(x = RT, y = mle, ymin = cilow, ymax = cihigh, color = fraction)) +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary only" = "#CA2FB4", "Shared" = "#CA932F", "Recurrence only" = "#2fb3ca"), drop=F) +
  #facet_wrap( ~ subtype) +
  geom_hline(yintercept = 1) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim=c(0.75,2)) +
  labs(x = "RT", y = "Global dN/dS", color = "Variant Fraction")

plot(p)


###################
## Plot dNdS-CV genes by RT - fraction
###################

dat <- dnds_fraction_sel_cv %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,received_rt, fill = list(qind_cv=1)) %>%
  group_by(fraction,received_rt) %>%
  arrange(qind_cv, pind_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots <- lapply(split(dat, paste(dat$received_rt, dat$fraction))[c(1,5,3,2,6,4)], function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:7,], aes(x=gene_name, y=-log10(qind_cv), fill=fraction), width=0.75) +
    geom_bar(aes(y=-log10(pind_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(FDR)", x="") +
    guides(fill = F) +
    facet_wrap(~ received_rt, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#CA2FB4","S" = "#CA932F", "R" = "#2fb3ca"), drop=F) +
    theme_bw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0)) 
  
  return(p)
})

#### Export
pdf(file = "figures/supp/GLASS_dndscv_IDHmut+IDHwt.pdf", height = 5, width = 10, bg = "transparent", useDingbats = FALSE)
do.call(grid.arrange,(c(myplots, ncol=3)))
dev.off()
###################
## Plot dNdS-CV genes by RT - fraction (IDHmut only)
###################

dat <- dnds_fraction_idhmut_sel_cv %>%
  mutate(fraction = factor(fraction, levels = c("P", "S", "R"))) %>%
  complete(gene_name,fraction,received_rt, fill = list(qind_cv=1)) %>%
  group_by(fraction,received_rt) %>%
  arrange(qind_cv, pind_cv) %>%
  ungroup() %>%
  arrange(fraction)

myplots <- lapply(split(dat, paste(dat$received_rt, dat$fraction))[c(1,5,3,2,6,4)], function(df){
  df$gene_name = factor(df$gene_name, levels = unique(df$gene_name))
  p <- ggplot(df[1:7,], aes(x=gene_name, y=-log10(qind_cv), fill=fraction), width=0.75) +
    geom_bar(aes(y=-log10(pind_cv)), fill = "gray90", stat="identity", position = position_dodge()) + 
    geom_bar(stat="identity", position = position_dodge()) + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2) +
    labs(y="-log10(FDR)", x="") +
    guides(fill = F) +
    facet_wrap(~ received_rt, scales = "free_y") + 
    coord_flip(ylim = c(0,10)) +
    scale_fill_manual(values = c("P" = "#CA2FB4","S" = "#CA932F", "R" = "#2fb3ca"), drop=F) +
    theme_bw(base_size = 16) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_line(size = 0))
  
  return(p)
})

#### Export
pdf(file = "figures/supp/GLASS_dndscv_IDHmut.pdf", height = 5, width = 10, bg = "transparent", useDingbats = FALSE)
do.call(grid.arrange,(c(myplots, ncol=3)))
dev.off()

