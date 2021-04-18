### Validate frindings from GLASS and HMF in data from Kucab et al. 2019, Cell ####
# Data downloaded from: https://data.mendeley.com/datasets/m7r4msjb4c/2
# Extended Data Figures 1i and 1k

# Load necessary packages
library(tidyverse)
library(DBI)
library(pracma)
library(ggpubr)
library(ggrepel)

#
meta_query <- readxl::read_xlsx("data/Kucab_meta.xlsx")
indel_query <- read_tsv("data/Kucab_indels.txt")
snv_query <- read_tsv("data/Kucab_snvs.txt")

meta <- meta_query %>% select(Sample.Name = `Sample Name`, Group, Compound) %>% 
  mutate(Group = ifelse(Group == "DNA Damage Response Inhibitors", "DDR Inh.", Group)) %>% 
  mutate(Group = ifelse(Compound == "Simulated solar radiation", "UV", Group))

indel <- indel_query %>% 
  mutate(class = ifelse(indel.length == 1, "1bp", ">1bp")) %>% 
  select(Sample, Sample.Name, Type, indel.length, change, classification, class)

indel_count_simple <- indel %>% add_count(Sample, Type) %>% 
  select(Sample, Sample.Name, Type, n) %>% distinct() %>% 
  spread(Type, n, fill = 0) %>% gather(Type, n, -c("Sample", "Sample.Name"))

snv_count <- snv_query %>% add_count(Sample) %>% 
  select(Sample, Sample.Name, n) %>% distinct()

pdf(file = "figures/supp/FigS1I.pdf", height = 3, width = 4, bg = "transparent", useDingbats = FALSE)
indel_count_simple %>% left_join(meta) %>% filter(!is.na(Group)) %>% 
  mutate(Type = ifelse(Type == "Del", "Small deletions", "Small insertions")) %>% 
  filter(Group %in% c("Control", "Radiation")) %>% 
  ggplot(aes(x=Group, y=n)) + geom_boxplot(width=0.5,(aes(fill=Group))) + geom_jitter(width = 0.1, height = 0.1) + 
  facet_grid(~Type) +
  theme_bw() + stat_compare_means(label.x = 1.4, label.y.npc = 0.15,
                                  aes(label = sprintf("P = %2.1e", as.numeric(..p.format..)))) + 
  labs(x="", y="Count") + scale_fill_manual(values=c("white", "gray50")) +
  theme(strip.background =element_rect(fill="white"), legend.position = "none")
dev.off()

## Order by count
set.seed(123)

pdf(file = "figures/supp/FigS1K.pdf", height = 4, width = 8, bg = "transparent", useDingbats = FALSE)
indel_count_simple %>% left_join(meta) %>% filter(!is.na(Group)) %>% 
  filter(Type == "Del") %>% mutate(RT = ifelse(Group == "Radiation", "RT+", "RT-")) %>% 
  group_by(Group) %>% mutate(n_median=median(n), n_mean=mean(n), n_sd=sd(n)) %>% ungroup() %>% 
  ggplot(aes(x=reorder(Group, -n_mean), y=n)) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "errorbar",
               width = .2) + 
  stat_summary(fun.y = "mean", geom = "bar", alpha = .5, aes(fill = RT)) +
  stat_summary(fun.y = "mean", geom = "point",size = 1, alpha=.5) +
  stat_summary(fun.y = "median", geom = "point", color ="red", size = 1) + 
  scale_fill_manual(values = c("black", "red")) + 
  labs( x = "", y = "Small deletion count") + 
  theme_classic2() + theme(axis.text.x = element_text(angle=45,hjust=1),  legend.position = "none") + 
  EnvStats::stat_n_text(y.pos = 1, color="white", size=3.5)
dev.off()

### END ###