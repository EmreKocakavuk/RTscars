# Mutational signature analysis using SBS (single base substitution) and ID (indel) signatures from COSMIC/Sigprofiler

###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(ggpubr)
library(corrplot)
library(magrittr)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

hmf_meta <- dbReadTable(con, Id(schema="public", table="hmf_meta")) %>% 
  mutate(status_hrd_msi = factor(status_hrd_msi, levels = c("HR_proficient", "HRD", "MSI")),
         RT = factor(RT, levels = c("RT-", "RT+ pal", "RT+ cur")))

ID_sig_relative <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig")) %>% filter(sampleid %in% hmf_meta$sampleid)
ID_sig_absolute <- dbReadTable(con,  Id(schema="analysis", table="hmf_id_sig_absolute")) %>% filter(sampleid %in% hmf_meta$sampleid)

SBS_sig_relative <- dbReadTable(con,  Id(schema="analysis", table="hmf_sbs_sig")) %>% filter(sampleid %in% hmf_meta$sampleid)
SBS_sig_absolute <- dbReadTable(con,  Id(schema="analysis", table="hmf_sbs_sig_absolute")) %>% filter(sampleid %in% hmf_meta$sampleid)

## Exploratory analysis
ID_sig_relative %>% select(-sampleid) %>% cor() %>% corrplot()
SBS_sig_relative %>% select(-sampleid) %>% cor() %>% corrplot()
SBS_sig_relative %>% left_join(ID_sig_relative) %>% select(-sampleid) %>% cor() %>% corrplot()
SBS_sig_absolute %>% left_join(ID_sig_absolute) %>% select(-sampleid) %>% cor() %>% corrplot()

#############
#
#Functions
#

testWilcoxGroup <- function(df) {
  wtest = wilcox.test(df$value ~ df$obs, paired = F, conf.int = TRUE)
  data.frame(n = nrow(df),
             median_a = median(df$value[df$obs == 1], na.rm = T),
             median_b = median(df$value[df$obs == 0], na.rm = T),
             mean_a = mean(df$value[df$obs == 1], na.rm = T),
             mean_b = mean(df$value[df$obs == 0], na.rm = T),
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

do_wilcox <- function (df) 
{ df %>% 
  group_by(key) %>% 
  do(testWilcoxGroup(.)) %>%
  ungroup() %>% mutate(fdr = p.adjust(wilcox_p, method = "fdr"), fc_mean = mean_a/mean_b, 
                       log_or = log(or), 
                       log_fc_mean = log(fc_mean)) %>% 
  mutate(signif = ifelse(fdr > 0.1, NA, -log10(fdr))) %>% 
  mutate(label = ifelse(fdr > 0.1, NA, key))
}

level_id <- c("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7", "ID8", "ID9", "ID10", 
              "ID11", "ID12", "ID13", "ID14", "ID15", "ID16", "ID17")

level_sbs <- c("SBS1", "SBS2", "SBS3", "SBS4", "SBS5", "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d", 
               "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS11", "SBS12", "SBS13", "SBS14", "SBS15", "SBS16", 
               "SBS17a", "SBS17b", "SBS18", "SBS19", "SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25", 
               "SBS26", "SBS27", "SBS28", "SBS29", "SBS30", "SBS31", "SBS32", "SBS33", "SBS34", "SBS35", 
               "SBS36", "SBS37", "SBS38", "SBS39", "SBS40", "SBS41", "SBS42", "SBS43", "SBS44", "SBS45", 
               "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", "SBS52", "SBS53", "SBS54", "SBS55",
               "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", "SBS84", "SBS85")

selected_sbs <- c("SBS1", "SBS2", "SBS3", "SBS6", "SBS8", "SBS11", "SBS13", "SBS14", "SBS15", "SBS18", 
                  "SBS20", "SBS21", "SBS23", "SBS25", "SBS26", "SBS31", "SBS35", "SBS39", "SBS44")

level_sbs <- selected_sbs %>% as.data.frame() %>% mutate(dat=gsub(pattern = "SBS", "", .)) %>% pull(dat)

cp <- coord_polar()
cp$is_free <- function() TRUE

theme_test <- 
  theme(legend.position = c(0.5, 0), 
        legend.justification = c(0.5, 1), 
        legend.direction = "horizontal", 
        plot.title = element_text(hjust = 0.5),
        #legend.text=element_text(size=rel(0.6)), 
        #legend.key.size = unit(0.3, 'cm'), 
        legend.text = element_text(size = 11),
        legend.background=element_blank(), 
        axis.title.y = element_text(hjust=1, vjust=-2))

plot_id <- function(df) {
  df %>% 
    mutate(key=gsub(pattern = "ID", "", key), key=factor(key, levels = seq(1:17))) %>% 
    ggplot(aes(x = factor(key), y = log2(or))) +
    geom_bar(stat="identity", width = 1, aes(fill = signif)) + 
    coord_polar() + theme_light(base_size=12) + 
    scale_fill_gradient2(low = "darkgrey", high = "darkred") + 
    labs(y = expression(atop(paste(log[2]("OR")))), x="", fill= "") + theme_test
}

plot_sbs <- function(df) {
  df %>% filter(key %in% selected_sbs) %>% 
    mutate(key=gsub(pattern = "SBS", "", key), key = factor(key, levels = level_sbs)) %>% 
    ggplot(aes(x = factor(key), y = log2(or))) +
    geom_bar(stat="identity", width = 1, aes(fill = signif)) + 
    coord_polar() + theme_light(base_size=12) + 
    scale_fill_gradient2(low = "darkgrey", high = "darkred") + 
    labs(y = expression(atop(paste(log[2]("OR")))), x="", fill= "") + theme_test
}

ID_long <- ID_sig_relative %>% left_join(hmf_meta) %>% 
  gather(key, value, ID1:ID17)

SBS_long <- SBS_sig_relative %>% left_join(hmf_meta) %>% 
  gather(key, value,`SBS1`:`SBS85`)
  
#RT
rt_id <- ID_long %>% 
  filter(RT %in% c("RT+ cur", "RT-")) %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox
 
rt_sbs <- SBS_long %>% 
  filter(RT %in% c("RT+ cur", "RT-")) %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox

#MSI
msi_id <- ID_long %>% 
  mutate(obs = ifelse(status_hrd_msi=="MSI", 0, 1)) %>% do_wilcox

msi_sbs <- SBS_long %>% 
  mutate(obs = ifelse(status_hrd_msi=="MSI", 0, 1)) %>% do_wilcox

# HRD
hrd_id <- ID_long %>% 
  mutate(obs = ifelse(status_hrd_msi=="HRD", 0,1)) %>% do_wilcox

hrd_sbs <- SBS_long %>% 
  mutate(obs = ifelse(status_hrd_msi=="HRD", 0, 1)) %>% do_wilcox

# RT in MSI
rt_in_msi_id <- ID_long %>% 
  mutate(RT= ifelse(RT=="RT-", "RT-", "RT+")) %>% filter(status_hrd_msi=="MSI") %>% 
  #ilter(RT %in% c("RT+ cur", "RT-")status_hrd_msi == "MSI") %>% 
  group_by(key) %>% filter(sum(value)>0) %>% ungroup() %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox

rt_in_msi_sbs <- SBS_long %>% 
  mutate(RT= ifelse(RT=="RT-", "RT-", "RT+")) %>% filter(status_hrd_msi=="MSI") %>% 
  #filter(RT %in% c("RT+ cur", "RT-"), status_hrd_msi == "MSI") %>% 
  group_by(key) %>% filter(sum(value)>0) %>% ungroup() %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox

# RT in HRD
rt_in_hrd_id <- ID_long %>% 
  mutate(RT= ifelse(RT=="RT-", "RT-", "RT+")) %>% filter(status_hrd_msi=="HRD") %>% 
  #filter(RT %in% c("RT+ cur", "RT-"), status_hrd_msi == "HRD") %>% 
  group_by(key) %>% filter(sum(value)>0) %>% ungroup() %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox

rt_in_hrd_sbs <- SBS_long %>% 
  mutate(RT= ifelse(RT=="RT-", "RT-", "RT+")) %>% filter(status_hrd_msi=="HRD") %>% 
  #filter(RT %in% c("RT+ cur", "RT-"), status_hrd_msi == "HRD") %>% 
  group_by(key) %>% filter(sum(value)>0) %>% ungroup() %>% 
  mutate(obs = ifelse(RT=="RT-", 1, 0)) %>% do_wilcox

margin <- theme(plot.margin = unit(c(-3,0,-3,0), "lines"))
p1 <- rt_id %>% plot_id() + labs(title = "RT", y = "")  + margin
p2 <- msi_id %>% plot_id() + labs(title = "MSI", y="") + margin
p3 <- hrd_id %>% plot_id() + labs(title = "HRD", y="") + margin

t1 <- rt_sbs %>% plot_sbs() + labs(y = "") + margin
t2 <- msi_sbs %>% plot_sbs() + labs(y = "") + margin
t3 <- hrd_sbs %>% plot_sbs() + labs(y = "") + margin

#### Rebuttal: Show Platinum signatures for validation purposes
pdf(file = "figures/rebuttal/FigR11.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
SBS_sig_absolute %>% left_join(hmf_meta) %>% 
  gather(key, value, `SBS1`:`SBS85`) %>% 
  filter(key %in% c("SBS31", "SBS35"), !is.na(received_platinum)) %>% 
  mutate(platinum = ifelse(received_platinum == 1, "Platinum +", "Platinum -")) %>% 
  ggplot(aes(x=platinum, y=value, fill=factor(received_platinum))) + 
  geom_boxplot() + facet_grid(~key) + theme_bw() + 
  stat_compare_means(label.x = 1.5, label.y.npc = 0.05, label = "p.format") + 
  scale_fill_manual(values = c("white", "gray50")) + scale_y_log10() + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) + #EnvStats::stat_n_text() + 
  labs(y="Absolute contribution", x="")
dev.off()


#### Suppl. Figure 3f - Separation of absolute/relative ID8 by tumor type
own_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", 
                "#B09C85FF", "#ff7f00", "#6a3d9a")

supp_join <- function (x) {
  x %>% left_join(hmf_meta) %>% 
    mutate(RT = factor(RT,labels = c("RT-" = "RT-", "RT+ pal" = "RT+\npal", "RT+ cur" = "RT+\ncur"))) %>% 
    ggplot(aes(x=RT, y = ID8)) + 
    geom_boxplot(aes(fill = primary_tumor_location, alpha = factor(RT))) + #, outlier.shape = NA)
    theme_bw() + scale_y_log10() + facet_grid(~primary_tumor_location) +
    labs(x="", y="ID8") +
    stat_compare_means(label.x = 1.25, label.y.npc = 0.2,
                       aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) +
    scale_fill_manual(values = own_colors) +
    theme(legend.position = "none", strip.background =element_rect(fill="white"))
}

supp_id8_absolute <- 
  ID_sig_absolute %>%  supp_join() + labs(title = "Absolute")

supp_id8_relative <- 
  ID_sig_relative %>% supp_join() + labs(title = "Relative")

pdf(file = "figures/supp/FigS3F.pdf", height = 6, width = 14, bg = "transparent", useDingbats = FALSE)
ggarrange(supp_id8_absolute, supp_id8_relative, align = "hv", ncol=1)
dev.off()

###
### GLASS
###
# Create connection to database and load necessary tables
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq")) %>% select(aliquot_barcode, coverage_adj_mut_freq)
hm_freq <- mut_freq %>% filter(coverage_adj_mut_freq >=10, !grepl("-TP-", aliquot_barcode)) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  select(case_barcode) %>% distinct() %>% 
  filter(!case_barcode %in% c("GLSS-MD-0027", "GLSS-SF-0024", "GLSS-SU-0002", "GLSS-SU-0270", "TCGA-DU-6407")) #HM is only the 3rd surgery (exclude, since we focus only on surgery 1 and 2)

scars_set = dbReadTable(con,  Id(schema="analysis", table="scars_set")) %>% filter(!grepl("-DK-", aliquot_barcode)) %>% #filter(treatment_combi != "NA") %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"))
scars_pairs = dbReadTable(con,  Id(schema="analysis", table="scars_pairs")) %>% filter(!grepl("-DK-", case_barcode))
scars_final = scars_pairs %>% inner_join(scars_set) %>% select(-aliquot_barcode, -surgery_number) %>% 
  distinct() %>% filter(idh_codel_subtype != "IDHwt")

############################### ID/SBS SIGNATURE ANALYSIS #########################
subtypes = dbReadTable(con,  Id(schema="clinical", table="subtypes"))

scars_id_sig_query = dbReadTable(con,  Id(schema="analysis", table="scars_id_sig"))
id_sig <- scars_id_sig_query %>% filter(aliquot_barcode %in% c(scars_final$tumor_barcode_a, scars_final$tumor_barcode_b)) %>% 
  mutate(case_barcode = substr(aliquot_barcode,1,12)) %>% inner_join(scars_set) 

scars_id_sig_fraction_query = dbReadTable(con,  Id(schema="analysis", table="scars_id_sig_fraction"))
id_sig_fraction <- scars_id_sig_fraction_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(subtypes)


scars_id_sig_absolute_query = dbReadTable(con,  Id(schema="analysis", table="scars_id_sig_absolute"))
id_sig_abs <- scars_id_sig_absolute_query %>% filter(aliquot_barcode %in% c(scars_final$tumor_barcode_a, scars_final$tumor_barcode_b)) %>% 
  mutate(case_barcode = substr(aliquot_barcode,1,12)) %>% inner_join(scars_set)

scars_id_sig_fraction_absolute_query = dbReadTable(con,  Id(schema="analysis", table="scars_id_sig_fraction_absolute"))
id_sig_fraction_abs <- scars_id_sig_fraction_absolute_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(subtypes)

#
scars_sbs_sig_query = dbReadTable(con,  Id(schema="analysis", table="scars_sbs_sig"))
sbs_sig <- scars_sbs_sig_query %>% filter(aliquot_barcode %in% c(scars_final$tumor_barcode_a, scars_final$tumor_barcode_b)) %>% 
  mutate(case_barcode = substr(aliquot_barcode,1,12)) %>% inner_join(scars_set)

scars_sbs_sig_fraction_query = dbReadTable(con,  Id(schema="analysis", table="scars_sbs_sig_fraction"))
sbs_sig_fraction <- scars_sbs_sig_fraction_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(subtypes)


scars_sbs_sig_absolute_query = dbReadTable(con,  Id(schema="analysis", table="scars_sbs_sig_absolute"))
sbs_sig_abs <- scars_sbs_sig_absolute_query %>% filter(aliquot_barcode %in% c(scars_final$tumor_barcode_a, scars_final$tumor_barcode_b))

scars_sbs_sig_fraction_absolute_query = dbReadTable(con,  Id(schema="analysis", table="scars_sbs_sig_fraction_absolute"))
sbs_sig_fraction_abs <- scars_sbs_sig_fraction_absolute_query %>% filter(case_barcode %in% scars_final$case_barcode) %>% left_join(subtypes) 

###
ID_long <- id_sig_fraction %>% filter(fraction == "R") %>% 
  left_join(scars_final) %>% 
  gather(key, value, ID1:ID17)

SBS_long <- sbs_sig_fraction %>% filter(fraction == "R") %>% 
  left_join(scars_final) %>% 
  gather(key, value,`SBS1`:`SBS85`)

#RT
rt_id <- ID_long %>% 
  mutate(obs = ifelse(received_rt==0, 1, 0)) %>% do_wilcox

rt_sbs <- SBS_long %>% 
  mutate(obs = ifelse(received_rt==0, 1, 0)) %>% do_wilcox

#HM
hm_id <- ID_long %>% 
  mutate(obs = ifelse(HM=="HM", 0, 1)) %>% do_wilcox

hm_sbs <- SBS_long %>% 
  mutate(obs = ifelse(HM=="HM", 0, 1)) %>% do_wilcox

q1 <- rt_id %>% plot_id() + labs(title = "RT") + margin
q2 <- hm_id %>% plot_id() + labs(title = "HM", y = "") + margin

s1 <- rt_sbs %>% plot_sbs() + margin
s2 <- hm_sbs %>% plot_sbs() + labs(y = "") + margin


pdf(file = "figures/main/Fig3.pdf", height = 7, width = 13, bg = "transparent", useDingbats = FALSE)
ggpubr::ggarrange(q1, q2, NULL, p1, p2, p3, s1, s2, NULL, t1, t2, t3, 
                  align="hv", nrow=2, ncol=6, widths = c(1,1,0.1,1,1,1))
dev.off()

### END ###