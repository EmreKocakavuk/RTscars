## Sequence context and deletion length of small deletions in GLASS ##

###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(parallel)
library(ggpubr)
library(ggridges)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
mut_freq = dbReadTable(con,  Id(schema="analysis", table="mut_freq")) %>% select(aliquot_barcode, coverage_adj_mut_freq)
hm_freq <- mut_freq %>% filter(coverage_adj_mut_freq >=10, !grepl("-TP-", aliquot_barcode)) %>% mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  select(case_barcode) %>% distinct() %>% 
  filter(!case_barcode %in% c("GLSS-MD-0027", "GLSS-SF-0024", "GLSS-SU-0002", "GLSS-SU-0270", "TCGA-DU-6407")) #HM is only the 3rd surgery (exclude, since we focus only on surgery 1 and 2)

scars_set = dbReadTable(con,  Id(schema="analysis", table="scars_set")) %>% filter(!grepl("-DK-", aliquot_barcode)) %>% #filter(treatment_combi != "NA") %>% 
  mutate(HM = ifelse(case_barcode %in% hm_freq$case_barcode, "HM", "Non-HM"))
scars_pairs = dbReadTable(con,  Id(schema="analysis", table="scars_pairs")) %>% filter(!grepl("-DK-", case_barcode))
scars_final = scars_pairs %>% inner_join(scars_set) %>% select(-aliquot_barcode, -surgery_number) %>% distinct() %>% 
  filter(idh_codel_subtype != "IDHwt")

surgeries = dbReadTable(con,  Id(schema="clinical", table="surgeries"))
surgical_interval_2 <- surgeries %>% filter(surgery_number == 2, case_barcode %in% scars_final$case_barcode) %>% 
  select(case_barcode, surgical_interval_mo) %>% filter(!grepl("GLSS-DK-", case_barcode))

################
################ Detailed homology sequence with repeat annotation 
################
qres_delseq <- dbGetQuery(con, read_file("sql/del_sequence_context.sql"))

rep_dna <- c("DNA", "DNA?", "DNA/hAT", "DNA/hAT?", "DNA/hAT-Ac", "DNA/hAT-Blackjack", "DNA/hAT-Charlie", "DNA/hAT-Tag1", 
             "DNA?/hAT-Tip100?", "DNA/hAT-Tip100", "DNA/hAT-Tip100?", "DNA/Kolobok", "DNA/Merlin", "DNA/MULE-MuDR", 
             "DNA/PIF-Harbinger", "DNA?/PiggyBac?", "DNA/PiggyBac", "DNA/TcMar", "DNA/TcMar?", "DNA/TcMar-Mariner", 
             "DNA/TcMar-Pogo", "DNA/TcMar-Tc1", "DNA/TcMar-Tc2", "DNA/TcMar-Tigger")

rep_line <- c("LINE/CR1", "LINE/Dong-R4", "LINE/Jockey", "LINE/L1", "LINE/L1-Tx1", "LINE/L2", 
              "LINE/Penelope", "LINE/RTE-BovB", "LINE/RTE-X")

rep_ltr <- c("LTR", "LTR?", "LTR/ERV1", "LTR/ERV1?", "LTR/ERVK", "LTR/ERVL", 
              "LTR/ERVL?", "LTR/ERVL-MaLR", "LTR/Gypsy", "LTR/Gypsy?")

rep_rc <- c("RC?/Helitron?", "RC/Helitron")

rep_rna <- c("RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA")

rep_satelite <- c("Satellite", "Satellite/acro", "Satellite/centr", "Satellite/telo")

rep_sine <- c("SINE/5S-Deu-L2", "SINE/Alu", "SINE/MIR", "SINE?/tRNA", 
              "SINE/tRNA", "SINE/tRNA-Deu", "SINE/tRNA-RTE")


qres_delseq_clin <- qres_delseq %>% 
  rename(homology_sequence = homology) %>% 
  mutate(chrom = ifelse(chrom == 23, "X",chrom), 
         pos = as.numeric(pos), 
         deletion_length = nchar(actual_deletion), 
         deletion_length_without_homology = deletion_length-homology_length, 
         non_microhomology_sequence = substring(actual_deletion,nchar(homology_sequence) + 1,nchar(actual_deletion))) %>% 
  inner_join(scars_final) %>% 
  mutate(category_simple = case_when(deletion_length == 1 ~"1bp", 
                                     deletion_length > 1 & homology_length > 1 ~ "MH", 
                                     deletion_length > 1 & homology_length <= 1 ~ "Non-MH"), 
         homology_simple = case_when(category_simple %in% c("1bp", "Non-MH") ~ "Non-MH", 
                                     category_simple == "MH" & !is.na(mh_repeat_sequence) ~ mh_repeat_sequence, 
                                     category_simple == "MH" & is.na(mh_repeat_sequence) & nchar(homology_sequence) < 4 ~ homology_sequence, 
                                     TRUE ~ "complex"), 
         homology_simple_length = case_when(homology_simple == "Non-MH" ~ "Non-MH", 
                                            homology_simple == "complex" ~ "complex", 
                                            category_simple == "MH" & nchar(homology_simple) == 1 ~ "1NR", 
                                            category_simple == "MH" & nchar(homology_simple) == 2 ~ "2NR", 
                                            category_simple == "MH" & nchar(homology_simple) == 3 ~ "3NR")) %>% 
  mutate(reptype_simple = case_when(reptype %in% rep_dna ~ "DNA", 
                                    reptype %in% rep_line ~ "LINE", 
                                    reptype %in% rep_ltr ~ "LTR",
                                    reptype %in% rep_rc ~ "RC",
                                    reptype %in% rep_rna ~ "RNA", 
                                    reptype %in% rep_satelite ~"Satelite",
                                    reptype %in% rep_sine ~ "SINE", 
                                    TRUE ~ reptype)) %>% 
  filter(is.na(reptype))

########
######## Basic numbers
########
detailed <- qres_delseq_clin %>% 
  mutate(received_rt = ifelse(received_rt == 1, "RT+", "RT-")) %>% 
  filter(idh_codel_subtype != "IDHwt", 
         deletion_length <=20) %>% 
  mutate(category_length = case_when(deletion_length == 1 ~ "1bp",
                                     deletion_length == 2 ~ "2bp", 
                                     deletion_length == 3 ~ "3bp", 
                                     deletion_length == 4 ~ "4bp", 
                                     deletion_length >= 5 ~ "+5bp",
                                     TRUE ~ "NA")) %>% 
  filter(is.na(reptype))

detailed_tmp <- detailed %>% 
  group_by(case_barcode, fraction) %>% 
  mutate(mean_length = mean(deletion_length), 
         median_length = median(deletion_length)) %>% 
  ungroup() %>% 
  select(case_barcode, fraction, mean_length, median_length, idh_codel_subtype, received_rt, HM) %>% 
  distinct()

# Separate into Primary vs Recurrence -------------------------------------

d_primary <- detailed %>% filter(fraction %in% c("P", "S")) %>% mutate(fraction = "Primary")
d_rec <- detailed %>% filter(fraction %in% c("R", "S")) %>% mutate(fraction = "Recurrence")
d_both <- full_join(d_primary, d_rec) %>% filter(idh_codel_subtype != "IDHwt")

fig_s2a <- 
  detailed_tmp %>% filter(fraction == "R") %>% 
  ggplot(aes(x=received_rt, y = mean_length)) + 
  geom_boxplot(aes(fill=received_rt)) + 
  labs(x="", y="Mean length of newly\nacquired deletions (bp)", alpha = "") +
  theme_classic() + scale_fill_manual(values = alpha(c("white", "black"), .7)) + theme(legend.position = "none", strip.background =element_rect(fill="white")) +
  stat_compare_means(label.x = 1.4, label.y.npc= 0.9,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + 
  facet_grid(~fraction) # --> put into supplement

# Fig 2a long
fig2a_left <- 
  d_both %>% 
  group_by(case_barcode, fraction) %>% 
  mutate(mean_length = mean(deletion_length), 
         median_length = median(deletion_length)) %>% 
  ungroup() %>% 
  select(case_barcode, fraction, mean_length, median_length, idh_codel_subtype, received_rt, HM) %>% 
  distinct() %>% 
  ggplot(aes(x=fraction, y = mean_length)) + 
  geom_line(aes(group=case_barcode), alpha =0.5, size=2) +
  geom_point(aes(color=fraction), alpha=.5, size = 3) + 
  facet_grid(~received_rt, scales = "free") +
  labs(x="", y="Mean deletion length (bp)", alpha = "") +
  theme_bw() + 
  scale_color_brewer(type = "qual", palette = 6, direction = -1) + 
  theme(legend.position = "none", strip.background =element_rect(fill="white")) +
  stat_compare_means(label.x = 1.4, label.y.npc= 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..))))

fig2a_right <- 
  d_both %>% filter(deletion_length <=20, deletion_length >1, HM == "Non-HM", received_rt == "RT+") %>% 
  add_count(case_barcode, fraction, deletion_length,  name = "n") %>% 
  add_count(case_barcode, fraction, name = "total") %>% 
  mutate(prop = n/total) %>% 
  select(case_barcode, fraction, deletion_length, received_rt, HM, n, total, prop, idh_codel_subtype) %>% 
  distinct() %>% 
  ggplot(aes(x=factor(deletion_length), y=prop, color = fraction)) +
  stat_summary_bin(fun.data = mean_cl_normal, geom = "pointrange", fun.args = list(mult = 1), position = position_dodge2(width=0.5)) +
  #stat_summary(fun = mean, geom = "line", aes(group=factor(fraction)), position = position_dodge2(width=0.5)) +
  stat_compare_means(label="p.signif" , hide.ns = T, paired=T) + 
  scale_color_brewer(palette = 6, type = "qual", direction = -1) +
  scale_y_log10() + labs(y = "Proportion", x = "Deletion length (bp)", color = "") +
  theme_bw() #+ #facet_wrap(~reptype_simple)

# Numbers for piechart ----------------------------------------------------

library(magrittr)
detailed %$% table(.$category_simple)
detailed %$% table(.$category_simple, .$fraction)

# With Microhomology
detailed %>% filter(homology_length > 1) # Deletions with Microhomology 

detailed %>% filter(homology_length > 1, mh_num_repeats > 0) # Deletions with Microhomology and MH-repeats 
detailed %>% filter(homology_length > 1, is.na(mh_num_repeats)) # Deletions with Microhomology no MH-repeats 

detailed %>% filter(homology_length > 1, mh_num_repeats > 0, nchar(mh_repeat_sequence) == 1) # Deletions with Microhomology and 1nt MH-repeats 
detailed %>% filter(homology_length > 1, mh_num_repeats > 0, nchar(mh_repeat_sequence) == 2) # Deletions with Microhomology and 2nt MH-repeats 
detailed %>% filter(homology_length > 1, mh_num_repeats > 0, nchar(mh_repeat_sequence) == 3) # Deletions with Microhomology and 3nt MH-repeats 

detailed %>% filter(homology_length > 1, is.na(mh_num_repeats), nchar(homology_sequence) == 2) # Deletions with Microhomology and 2nt MH sequence (No repeat) 
detailed %>% filter(homology_length > 1, is.na(mh_num_repeats), nchar(homology_sequence) == 3) # Deletions with Microhomology and 3nt MH sequence (No repeat) 
detailed %>% filter(homology_length > 1, is.na(mh_num_repeats), nchar(homology_sequence) > 3) # Deletions with Microhomology and complex MH sequence (No repeat) 

# Without Microhomology
detailed %>% filter(homology_length <= 1) # Deletions without Microhomology 

detailed %>% filter(homology_length <= 1, deletion_length == 1) # Deletions without Microhomology and 1bp length 
detailed %>% filter(homology_length <= 1, deletion_length > 1) # Deletions without Microhomology and >1bp length 

#### Compare Primary vs Recurrence
fig_s2e <- 
  d_both %>% #filter(category_simple != "1bp") %>% 
  mutate(category_simple = ifelse(category_simple == "MH", ">1bp MH", ifelse(category_simple == "Non-MH", ">1bp Non-MH", category_simple))) %>% 
  add_count(case_barcode, fraction, name = "t") %>% 
  add_count(case_barcode, fraction, category_simple, name = "n") %>% 
  mutate(prop = n/t) %>% 
  select (case_barcode, fraction, category_simple, prop, HM, received_rt, idh_codel_subtype) %>% distinct() %>% 
  spread(category_simple, prop, fill= 0) %>% 
  gather(category_simple, prop, `1bp`:`>1bp Non-MH`) %>% 
  mutate(category_simple = factor(category_simple, levels = c(">1bp MH", "1bp", ">1bp Non-MH"))) %>% 
  ggplot(aes(x=fraction, y = prop)) +
  geom_line(aes(group=case_barcode, color = category_simple), linetype="solid", size=1, alpha = 0.5) +
  geom_point(aes(group=case_barcode)) + 
  scale_color_manual(values = c("1bp" = "#5d5f60", ">1bp Non-MH" = "#2d004b", ">1bp MH" = "#f6711f")) +
  facet_grid(received_rt ~ HM + category_simple, scales = "free_y") +
  stat_compare_means(label.x = 1.4, label.y.npc= 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + theme_bw() +
  theme(strip.background =element_rect(fill="white")) + labs(x="", y = "Proportion", color = "Deletion\ncharacteristic")

fig2d <- 
  d_both %>% filter(received_rt == "RT+", HM == "Non-HM") %>% 
  mutate(category_simple = ifelse(category_simple == "MH", ">1bp MH", ifelse(category_simple == "Non-MH", ">1bp Non-MH", category_simple))) %>% 
  add_count(case_barcode, fraction, name = "t") %>% 
  add_count(case_barcode, fraction, category_simple, name = "n") %>% 
  mutate(prop = n/t) %>% 
  select (case_barcode, fraction, category_simple, prop, HM, received_rt, idh_codel_subtype) %>% distinct() %>% 
  spread(category_simple, prop, fill= 0) %>% 
  gather(category_simple, prop, `1bp`:`>1bp Non-MH`) %>% 
  mutate(category_simple = factor(category_simple, levels = c(">1bp MH", "1bp", ">1bp Non-MH"))) %>% 
  ggplot(aes(x=fraction, y = prop)) +
  geom_line(aes(group=case_barcode, color = category_simple), linetype="solid", alpha = 0.5, size = 2) +
  geom_point(aes(group=case_barcode), alpha = 0.5, size=3) + 
  scale_color_manual(values = c("1bp" = "#5d5f60", ">1bp Non-MH" = "#2d004b", ">1bp MH" = "#f6711f")) +
  facet_grid(~ category_simple, scales = "free") +
  stat_compare_means(label.x = 1.4, label.y.npc= 0.9, paired=T,
                     aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + theme_bw() + 
  theme(strip.background =element_rect(fill="white"), legend.position = "none") + 
  labs(x="", y = "Proportion of deletions", color = "Deletion\ncharacteristic")

################################################################
################################################################
######### Export figures
pdf(file = "figures/main/Fig2_GLASS_del_length.pdf", height = 5, width = 12, bg = "transparent", useDingbats = FALSE)
ggarrange(ggarrange(fig2a_left), ggarrange(fig2a_right), align="hv", widths = c(1,2))
dev.off()

pdf(file = "figures/supp/FigS2A_GLASS_reconly_del_length.pdf", height = 3, width = 3, bg = "transparent", useDingbats = FALSE)
fig_s2a
dev.off()

pdf(file = "figures/supp/FigS2E_GLASS_del_homology.pdf", height = 5, width =10 , bg = "transparent", useDingbats = FALSE)
fig_s2e
dev.off()

pdf(file = "figures/main/Fig2_GLASS_del_homology.pdf", height = 4, width =6 , bg = "transparent", useDingbats = FALSE)
fig2d
dev.off()

### END ###