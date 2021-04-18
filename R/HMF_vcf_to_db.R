## ---
##  Created: May 1, 2020
##  Updated: May 4, 2020
##  Author: FP Barthel
## ---
##  Load and parse Hartwig data and push to GLASS db
## ---

library(VariantAnnotation)
library(tidyverse)
library(DBI)
library(BSgenome.Hsapiens.UCSC.hg19)

## ---
##
##  PROCESS METADATA
##
## ---

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
meta <- read_tsv("data/hartwig-update1//metadata/metadata.tsv", locale = locale(encoding = "latin1"))

meta <- meta %>% select(sampleid = sampleId,
                        patientid = `#patientId`,
                        primary_tumor_location = primaryTumorLocation,
                        primary_tumor_subtype = cancerSubtype,
                        biopsy_date = biopsyDate,
                        biopsy_site = biopsySite,
                        biopsy_location = biopsyLocation,
                        sex = gender,
                        birth_year = birthYear,
                        death_date = deathDate,
                        systemic_pretreatment = hasSystemicPreTreatment,
                        radiotherapy_pretreatment = hasRadiotherapyPreTreatment,
                        treatment_given = treatmentGiven,
                        treatment_start_date = treatmentStartDate,
                        treatment_end_date = treatmentEndDate,
                        treatment = treatment,
                        treatment_type = treatmentType,
                        response_date = responseDate,
                        response_measured = responseMeasured,
                        first_response = firstResponse)

dbWriteTable(con, DBI::Id(schema = "public", table="hwm"), meta, overwrite = TRUE)

## ---
##
##  PROCESS VARIANTS
##
## ---

filterVCF <- function(vcf) {
  
  ## extract samples
  vcfsamples <- samples(header(vcf))
  stopifnot(length(vcfsamples)==1)
  
  ## filter FILTER = PASS
  vcf <- subset(vcf, rowRanges(vcf)$FILTER == "PASS")
  
  ## subset known chromosomes
  vcf <- subset(vcf, seqnames(vcf) %in% c(1:22,"X"))
  
  return(vcf)
}

getVcf <- function(vcf) {
  sampleid <- samples(header(vcf))
  stopifnot(length(sampleid)==1)
  
  vcf <- as(vcf, "VRanges")
  
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, sprintf("chr%s",as.character(seqnames(vcf))), start(vcf)-19, end(vcf)+20, as.character = TRUE)
  
  chrom <- case_when(as.character(seqnames(vcf)) %in% as.character(1:22) ~ as.integer(seqnames(vcf)),
                     as.character(seqnames(vcf)) == "X" ~ as.integer(23),
                     as.character(seqnames(vcf)) == "Y" ~ as.integer(24),
                     TRUE ~ NA_integer_)
  
  res <- tibble(sampleid,
                    eventid = 1:length(vcf),
                    chrom = chrom,
                    pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                    ref = ref(vcf),
                    alt = alt(vcf),
                    variant_type = case_when(nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 1 ~ "SNP",
                                             nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 2 ~ "DNP",
                                             nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 3 ~ "TNP",
                                             nchar(ref(vcf)) < nchar(alt(vcf)) & nchar(ref(vcf)) == 1 ~ "INS",
                                             nchar(ref(vcf)) < nchar(alt(vcf)) & nchar(ref(vcf)) > 1 ~ "INSDEL",
                                             nchar(ref(vcf)) > nchar(alt(vcf)) & nchar(alt(vcf)) == 1 ~ "DEL",
                                             nchar(ref(vcf)) > nchar(alt(vcf)) & nchar(alt(vcf)) > 1 ~ "DELINS",
                                             TRUE ~ NA_character_),
                    genotype = vcf$GT,
                    depth = totalDepth(vcf),
                    ref_count = refDepth(vcf),
                    alt_count = altDepth(vcf),
                    af = round(altDepth(vcf) / totalDepth(vcf),4),
                    seq = as.character(seq))
  
  write_tsv(res, path = sprintf("data/hartwig-update1/dbdump/%s.tsv", sampleid))
  #return(res)
  rm(res, vcf)
}

processVCF <- function(vcf_file) {
  message(vcf_file)
  ## read VCF
  vcf <- readVcf(vcf_file, "hg19")
  vcf <- filterVCF(vcf)
  getVcf(vcf)
  rm(vcf)
  #return(vcf)
}

vcff <- list.files(path = "data/hartwig-update1/leftaligned", recursive = TRUE, pattern = "somatic.normalized.sorted.vcf.gz$", full.names = TRUE)

mclapply(vcff, processVCF, mc.cores = 12)

#con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")
#dbWriteTable(con, DBI::Id(schema = "public", table="hwg2"), dat, overwrite = TRUE)


## ---
##
##  PROCESS DRUGS
##
## ---
## Taxanes
## Taxanes are a class of diterpenes. They were originally identified from plants of the genus Taxus (yews), and feature a taxadiene core. Paclitaxel (Taxol) and docetaxel (Taxotere) are widely used as chemotherapy agents. Cabazitaxel was FDA approved to treat hormone-refractory prostate cancer.
## ---
## pyrimidine antagonists
## The most commonly used pyrimidine antagonists are 5-fluorouracil (5-FU), gemcitabine (dFdC) and cytarabine (ara-C). Newer oral variants of 5-FU are capecitabine and tegafur
## ---
## Anthracyclines
## Anthracyclines is a class of drugs used in cancer chemotherapy that are extracted from Streptomyces bacterium. These compounds are used to treat many cancers, including leukemias, lymphomas, breast, stomach, uterine, ovarian, bladder cancer, and lung cancers.
## ---
## Alkylating
## Alkylating agents are compounds that work by adding an alkyl group to the guanine base of the DNA molecule, preventing the strands of the double helix from linking as they should. This causes breakage of the DNA strands, affecting the ability of the cancer cell to multiply.
## ---
## Platinum
## Platinum-based antineoplastic drugs (informally called platins) are chemotherapeutic agents used to treat cancer. They are coordination complexes of platinum. These drugs are used to treat almost half of people receiving chemotherapy for cancer.
## ---

drugs <- read_tsv("data/hartwig/metadata/pre_biopsy_drugs_by_patient.tsv", locale = locale(encoding = "latin1"))
table(drugs$mechanism)


drugs2 <- drugs %>% rename(patientid = `#patientId`) %>%
  group_by(patientid) %>%
  summarize(received_pyrimidine_antagonist = any(str_detect(mechanism, fixed("pyrimidine (ant)agonist", ignore_case = TRUE))),
            received_taxane = any(str_detect(mechanism, fixed("taxane", ignore_case = TRUE))),
            received_anthracyclines = any(str_detect(mechanism, fixed("anthracycline", ignore_case = TRUE))),
            received_alkylating = any(str_detect(mechanism, fixed("alkylating", ignore_case = TRUE))),
            received_platinum = any(str_detect(mechanism, fixed("platinum", ignore_case = TRUE))),
            received_hormonal_therapy = any(str_detect(type, regex("hormonal|deprivation", ignore_case = TRUE))),
            received_targeted_therapy = any(str_detect(type, fixed("targeted", ignore_case = TRUE))),
            received_immunotherapy = any(str_detect(type, fixed("immunotherapy", ignore_case = TRUE)))) %>%
  ungroup()

table(drugs2$received_pyrimidine_antagonist)
table(drugs2$received_taxane)
table(drugs2$received_anthracyclines)
table(drugs2$received_alkylating)
table(drugs2$received_platinum)
table(drugs2$received_hormonal_therapy)
table(drugs2$received_targeted_therapy)
table(drugs2$received_immunotherapy)

dbWriteTable(con, DBI::Id(schema = "public", table="hwd"), drugs2, overwrite = TRUE)

## ---
##
##  PROCESS CNVs
##
## ---

processCNV <- function(f) {
  options(scipen=999)
  dat <- read_tsv(f, col_types = cols(chromosome = col_character()))
  dat <- dat %>%
    transmute(sampleid = unlist(map(str_split(basename(f), "\\."),1)),
              chrom = ifelse(chromosome == "X", 23, ifelse(chromosome == "Y", 24, as.integer(chromosome))),
              pos = sprintf("[%s,%s]", start, end),
              copy_number = copyNumber,
              baf_count = bafCount,
              observed_baf = observedBAF,
              segment_start_support = segmentStartSupport,
              segment_end_support = segmentEndSupport,
              method,
              depth_window_count = depthWindowCount,
              gc_content = gcContent,
              min_start = minStart,
              max_start = maxStart,
              minor_allele_ploidy = minorAllelePloidy,
              major_allele_ploidy = majorAllelePloidy)
  return(dat)
}

cnvf <- list.files(path = "data/hartwig/somatics", recursive = TRUE, pattern = "cnv.somatic.tsv$", full.names = TRUE)

cnv_all <- mclapply(cnvf, processCNV, mc.cores = 16)
cnv <- data.table::rbindlist(cnv_all)

dbWriteTable(con, DBI::Id(schema = "public", table="hwc"), cnv, overwrite = TRUE)

## ---
##
##  PROCESS SVs
##
## ---

svvcff <- list.files(path = "data/hartwig/somatics", recursive = TRUE, pattern = "sv.ann.vcf.gz$", full.names = TRUE)

svFilterVCF <- function(vcf) {
  
  ## extract samples
  vcfsamples <- samples(header(vcf))
  stopifnot(length(vcfsamples)==1)
  
  ## filter FILTER = PASS
  vcf <- subset(vcf, rowRanges(vcf)$FILTER == "PASS")
  
  ## subset known chromosomes
  vcf <- subset(vcf, seqnames(vcf) %in% c(1:22,"X"))
  
  return(vcf)
}

svGetVcf <- function(vcf) {
  sampleid <- samples(header(vcf))
  stopifnot(length(sampleid)==1)
  
  vcf <- as(vcf, "VRanges")
  
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, sprintf("chr%s",as.character(seqnames(vcf))), start(vcf)-19, end(vcf)+20)
  
  chrom <- case_when(as.character(seqnames(vcf)) %in% as.character(1:22) ~ as.integer(seqnames(vcf)),
                     as.character(seqnames(vcf)) == "X" ~ as.integer(23),
                     as.character(seqnames(vcf)) == "Y" ~ as.integer(24),
                     TRUE ~ NA_integer_)
  
  res <- tibble(sampleid,
                    eventid = 1:length(vcf),
                    chrom = chrom,
                    pos = sprintf("[%s,%s]", start(vcf), end(vcf)),
                    ref = ref(vcf),
                    alt = alt(vcf),
                    variant_type = case_when(nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 1 ~ "SNP",
                                             nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 2 ~ "DNP",
                                             nchar(ref(vcf)) == nchar(alt(vcf)) & nchar(ref(vcf)) == 3 ~ "TNP",
                                             nchar(ref(vcf)) < nchar(alt(vcf) & nchar(ref(vcf)) == 1) ~ "INS",
                                             nchar(ref(vcf)) < nchar(alt(vcf) & nchar(ref(vcf)) > 1) ~ "INSDEL",
                                             nchar(ref(vcf)) > nchar(alt(vcf) & nchar(alt(vcf)) == 1) ~ "DEL",
                                             nchar(ref(vcf)) > nchar(alt(vcf) & nchar(alt(vcf)) > 1) ~ "DELINS",
                                             TRUE ~ NA_character_),
                    genotype = vcf$GT,
                    depth = totalDepth(vcf),
                    ref_count = refDepth(vcf),
                    alt_count = altDepth(vcf),
                    af = round(altDepth(vcf) / totalDepth(vcf),4),
                    seq = seq)
  
  return(res)
}

svProcessVCF <- function(vcf_file) {
  message(vcf_file)
  ## read VCF
  vcf <- readVcf(vcf_file, "hg19")
  vcf <- svFilterVCF(vcf)
  
  return(getVcf(vcf))
}
