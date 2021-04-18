###########################################
################ Settings #################
library(tidyverse)
library(DBI)

# Create connection to database and load necessary tables
con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

hmf_meta_query <- dbGetQuery(con, "SELECT 
	hwm.sampleid, 
	hwm.patientid,
	primary_tumor_location,
	primary_tumor_subtype,
	extract(year from age(biopsy_date,to_date(birth_year::varchar, 'YYYY')))::integer AS age,
	(CASE WHEN radiotherapy_pretreatment = 'Yes' THEN 1 WHEN radiotherapy_pretreatment = 'No' THEN 0 ELSE NULL END)::integer AS radiotherapy_pretreatment,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN 1 WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END)::integer AS systemic_pretreatment,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_pyrimidine_antagonist::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_pyrimidine_antagonist,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_taxane::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_taxane,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_anthracyclines::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_anthracyclines,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_alkylating::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_alkylating,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_platinum::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_platinum,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_hormonal_therapy::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_hormonal_therapy,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_targeted_therapy::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_targeted_therapy,
	(CASE WHEN systemic_pretreatment = 'Yes' THEN received_immunotherapy::integer WHEN systemic_pretreatment = 'No' THEN 0 ELSE NULL END) AS received_immunotherapy
FROM hwm 
LEFT JOIN hwd ON hwd.patientid = hwm.patientid")

# Load HMF metadata for linking
hmf_meta_link <- read_tsv("data/hmf_metadata.tsv") %>% 
  select(patientid = `#patientId`, sampleid = sampleId, sample = hmfSampleId)

# Load HRD/MSI info from HMF
hmf_hrd_msi_query <- read_tsv("data/HMF_predict_HRD_MSI.txt")
hmf_hrd_msi_df <- hmf_hrd_msi_query %>% left_join(hmf_meta_link)
hmf_hrd_msi <- hmf_hrd_msi_df %>% 
  mutate(status_hrd_msi = case_when(hr_status == "HR_deficient" ~ "HRD", 
                                    hr_status == "HR_proficient" ~ "HR_proficient",
                                    hr_status == "cannot_be_determined" & remarks_hr_status == "<50 indels" ~ "low_ID", 
                                    hr_status == "cannot_be_determined" & remarks_hr_status == "Has MSI (>14000 indel.rep)" ~ "MSI", 
                                    TRUE ~ "NA")) %>% 
  select(patientid, sampleid, status_hrd_msi)

# Load HMF RT data
hmf_rt <- readxl::read_xlsx("data/HMF_RT_20201208.xlsx")

hmf_rt_tmp <- hmf_rt %>%
  mutate(patientid = gsub(pattern = "-", replacement = "", SubjectKey)) %>% 
  filter(aim == "1=Curative")

hmf_rt_tmp2 <- hmf_rt %>% filter(!is.na(aim)) %>% 
  mutate(dose = as.numeric(dose), 
         patientid = gsub(pattern = "-", replacement = "", SubjectKey)) %>%  
  mutate(aim = ifelse(aim == "1=Curative", "cur", "pal")) %>% 
  group_by(patientid, aim) %>% 
  mutate(dose_sum = sum(dose, na.rm = T)) %>% 
  ungroup() %>% mutate(dose_sum = ifelse(dose_sum == 0, NA, dose_sum)) %>% 
  select(patientid, aim, dose_sum) %>% distinct() %>% 
  spread(key=aim, value = dose_sum) %>% 
  mutate(curative = ifelse(patientid %in% c(hmf_rt_tmp$patientid), "curative", "non_curative")) 

hmf_meta <- hmf_meta_query %>% left_join(hmf_rt_tmp2) %>% 
  mutate(radiotherapy_pretreatment = ifelse(radiotherapy_pretreatment != 1, radiotherapy_pretreatment, curative)) %>% 
  mutate(RT = ifelse(radiotherapy_pretreatment == 0, "RT-", ifelse(radiotherapy_pretreatment == "curative", "RT+ cur", "RT+ pal")), 
         RT = factor(RT, levels = c("RT-", "RT+ pal", "RT+ cur"))) %>% filter(!is.na(RT)) %>% 
  select(sampleid, patientid, primary_tumor_location, primary_tumor_subtype, systemic_pretreatment:RT, -curative, age) %>% 
  distinct() %>% mutate(primary_tumor_location = 
                          ifelse(primary_tumor_location %in% c("Bone/Soft tissue", "Breast", "Colon/Rectum", 
                                                               "Esophagus", "Head and neck", "Lung", "Nervous system", 
                                                               "Prostate", "Skin", "Urinary tract", "Uterus"), primary_tumor_location, "Others"), 
                        radiotherapy_pretreatment = ifelse(RT=="RT-", 0, 1)) %>% 
  left_join(hmf_hrd_msi)

hmf_meta %>% select(-primary_tumor_location) %>% anti_join(hmf_meta_query) #joins with query
table(hmf_meta$status_hrd_msi)

hmf_meta_final <- hmf_meta %>% filter(status_hrd_msi!= "low_ID") %>% 
  filter(!sampleid %in% c("CPCT02040173T", "CPCT02290016T", "CPCT02230024T", "CPCT02050296T", 
                         "CPCT02380037T", "CPCT02020413TII", "CPCT02020821T"))

# Upload to db
dbWriteTable(con, Id(schema="public",table="hmf_meta_test"), hmf_meta, overwrite=T)
dbWriteTable(con, Id(schema="public",table="hmf_meta"), hmf_meta_final, overwrite=T)

