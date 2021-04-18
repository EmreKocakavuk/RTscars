###########################################
################ Settings #################
library(tidyverse)
library(DBI)
library(pracma)

# Create connection to database
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

# Modify function for signature fitting
fit_to_signatures = function(mut, mut_signatures)
{
  # make sure dimensions of input matrix are correct
  if (dim(mut)[1] != dim(mut_signatures)[1])
    stop(paste("Mutation matrix and signatures input have",
               "different number of mutational features"))
  
  n_features = dim(mut)[1]
  n_samples = dim(mut)[2]
  n_signatures = dim(mut_signatures)[2]
  lsq_contribution = matrix(NA, nrow=n_signatures, ncol=n_samples)
  lsq_reconstructed = matrix(NA, nrow=n_features, ncol=n_samples)
  
  # Process each sample
  for (i in 1:ncol(mut))
  {
    y = mut[,i]
    lsq = lsqnonneg(mut_signatures, y)
    lsq_contribution[,i] = lsq$x
    lsq_reconstructed[,i] = mut_signatures %*% as.matrix(lsq$x)
  }
  
  # Add row and col names
  sample_names = colnames(mut)
  signature_names = colnames(mut_signatures)
  mut_type_names = rownames(mut_signatures)
  
  colnames(lsq_contribution) = sample_names
  rownames(lsq_contribution) = signature_names
  
  colnames(lsq_reconstructed) = sample_names
  rownames(lsq_reconstructed) = mut_type_names
  
  res = list(lsq_contribution, lsq_reconstructed)
  names(res) = c("contribution", "reconstructed")
  
  return(res)
}

# Load reference signatures from Sigprofiler
ID_signatures <- read_csv("data/sigprofiler_hmf/ref/sigProfiler_ID_signatures.csv") %>% dplyr::rename(rowname = `Mutation Type`) %>% mutate(rowname = gsub("[+]", "", rowname), rowname = gsub("repeats", "R", rowname), rowname = gsub("MH", "M", rowname), rowname = gsub("DEL", "Del", rowname), rowname = gsub("INS", "Ins", rowname)) %>% 
  separate(rowname, c("A", "B", "C", "D")) %>% unite(col = `Signatures`, `C`, `A`, `B`, `D`, sep=":") %>% arrange(`Signatures`) %>% select(-1) %>% as.matrix()

SBS_signatures <- read_csv("data/submission/sigprofiler_hmf/ref/sigProfiler_SBS_signatures_2019_05_22.csv") %>% select(-1,-2) %>% as.matrix()

# Load ID signatures from HMF
ID_HMF <- read_tsv("data/sigprofiler_hmf/ID/HMF.ID83.all") %>% 
  arrange(MutationType) %>% column_to_rownames("MutationType") %>% as.matrix()

SBS_HMF <- read_tsv("data/sigprofiler_hmf/SBS/HMF.SBS96.all") %>%  
  gather(key = key, value = value, -`MutationType`) %>% 
  spread(key = key, value = value) %>% mutate(Type = substr(MutationType, 3,5)) %>% 
  arrange(Type) %>% dplyr::select(-ncol(.)) %>% column_to_rownames("MutationType") %>% as.matrix()

## Fit ID signature
fit <- fit_to_signatures(ID_HMF, ID_signatures)

fit_tmp <- t(fit$contribution)/rowSums(t(fit$contribution))

ID_sig_relative <- rownames_to_column(as.data.frame(fit_tmp), var = "sampleid")

ID_sig_absolute <- rownames_to_column(as.data.frame(t(fit$contribution)), var = "sampleid")

## Fit SBS signature
fit2 <- fit_to_signatures(SBS_HMF, SBS_signatures)

fit2_tmp <- t(fit2$contribution)/rowSums(t(fit2$contribution))

SBS_sig_relative <- rownames_to_column(as.data.frame(fit2_tmp), var = "sampleid")

SBS_sig_absolute <- rownames_to_column(as.data.frame(t(fit2$contribution)), var = "sampleid")

# Upload signatures to db

dbWriteTable(con, Id(schema="analysis",table="hmf_id_sig"), ID_sig_relative, overwrite=T)
dbWriteTable(con, Id(schema="analysis",table="hmf_id_sig_absolute"), ID_sig_absolute, overwrite=T)

dbWriteTable(con, Id(schema="analysis",table="hmf_sbs_sig"), SBS_sig_relative, overwrite=T)
dbWriteTable(con, Id(schema="analysis",table="hmf_sbs_sig_absolute"), SBS_sig_absolute, overwrite=T)

#END