library(tidyverse)
library(janitor)
library(lubridate)
library(PheWAS)


#0. get studies to include
clindata_dir <- '/path/to/clinical_data/'
pop_dir <- '/path/to/population_data/'

study_data_sarcopenia <- read_csv(paste(pop_dir, 'study_data_sarcopenia_abstract.csv', sep=''))



#PheWAS study

#1. get diagnoses for phenotyping
dx <- read_csv_chunked(paste(clindata_dir, '1/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd') %>%
  bind_rows(read_csv_chunked(paste(clindata_dir, '2/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd')) %>%
  bind_rows(read_csv_chunked(paste(clindata_dir, '3/diagnoses.csv', sep=''), DataFrameCallback$new(dx_in_final), col_types='dcdcccccccd'))

#1.1 window and aggregate code counts
#ICD10->ICD9 equivalencies
icd10_9_equiv <- read_table('/path/to/diagnosis_gems_2018/2018_I10gem.txt',col_names=c("icd10_code", "icd9_code","flags"))

dx_count_allcodes <- dx %>%
  clean_names %>%
  select(patient_id, age, icd9_code, icd10_code) %>%
  mutate(code_type=if_else(is.na(icd9_code) | str_detect(icd9_code, ','), "ICD10CM","ICD9CM"),
         code=if_else(is.na(icd9_code) | str_detect(icd9_code, ','), icd10_code, icd9_code))

#1.2 Get medical phenotypes (phecodes)

phecodes <- dx_count_allcodes %>% 
  mutate(patient_id = as.character(patient_id)) %>%
  left_join(select(study_data_sarcopenia, patient_id, age_at_scan), by=c("patient_id"="patient_id")) %>%
  filter(age < age_at_scan +3/12 & age>age_at_scan-1/12) %>%
  select(patient_id, code_type, code) %>%
  rename(vocabulary_id=code_type) %>%
  mapCodesToPhecodes(make.distinct=FALSE) %>%
  group_by(patient_id, phecode) %>% 
  summarize(n=n()) %>% 
  ungroup

custom_agg <- function(index){
  s <- sum(index)
  if (s >1){
    return(s)
  } else if (s == 0){
    return(s)
  } else {
    return(NA)
  }
}

study_data_phenotypes <- createPhenotypes(phecodes %>% mutate(vocabulary_id='phecode') %>% rename(code=phecode) %>% select(patient_id, vocabulary_id, code, n) %>% as.data.frame,
                                          id.sex=study_data_sarcopenia %>% mutate(sex=if_else(sex=='Female',"F","M")) %>% select(patient_id, sex) %>% as.data.frame,
                                          full.population.ids = study_data_sarcopenia %>% select(patient_id) %>% pull,
                                          aggregate.fun=custom_agg,
                                          translate=FALSE)
#2. get genotypes (SMI and SMD)
mean_sds <- study_data_sarcopenia %>% 
  select(sex, recent_height_cm, cross_sectional_area_muscle, hounsfield_unit_muscle) %>% #csa_muscle in mm^2
  mutate(recent_height_cm = as.numeric(recent_height_cm),
         l3smi = (cross_sectional_area_muscle/100)/((recent_height_cm/100)**2)) %>%
  group_by(sex) %>%
  summarize(avg_csi_muscle = mean(l3smi),
            sd_csi_muscle = sd(l3smi),
            avg_hu_muscle = mean(hounsfield_unit_muscle, na.rm = TRUE),
            sd_hu_muscle = sd(hounsfield_unit_muscle, na.rm = TRUE))
study_data_scores <- select(study_data_sarcopenia, patient_id, anon_id, sex, recent_height_cm) %>%
  left_join(seg_data, by="anon_id") %>%
  mutate(id = patient_id,
         recent_height_cm = as.numeric(recent_height_cm),
         l3smi = (cross_sectional_area_muscle/100)/((recent_height_cm/100)**2),
         smi_zscore = if_else(sex=='Female', (l3smi-mean_sds$avg_csi_muscle[1])/mean_sds$sd_csi_muscle[1],(l3smi-mean_sds$avg_csi_muscle[2])/mean_sds$sd_csi_muscle[2]),
         sm_hu_zscore = if_else(sex=='Female', (hounsfield_unit_muscle-mean_sds$avg_hu_muscle[1])/mean_sds$sd_hu_muscle[1],(hounsfield_unit_muscle-mean_sds$avg_hu_muscle[2])/mean_sds$sd_hu_muscle[2])) %>%
  select(id, smi_zscore, sm_hu_zscore)

#3. get covariates used for regression
study_data_covariates <- select(study_data_sarcopenia, patient_id, age_at_scan, sex) %>%
  mutate(id=patient_id,
         age = as.numeric(age_at_scan), 
         sex=if_else(sex=='Male', 1,0)) %>%
  select(id, age, sex) %>%
  as.data.frame

#4. PheWAS analysis
sm_phewas <-  phewas(phenotypes = study_data_phenotypes %>%rename(id=patient_id), 
                     genotypes=study_data_scores,
                     covariates=study_data_covariates,
                     additive.genotypes=FALSE,
                     significance.threshold = "bonferroni",
                     min.records=100, 
                     cores=4, 
                     MASS.confint.level =.95)

all_results <- sm_phewas %>%
  select(phenotype, snp, p, OR, lower, upper, n_total, n_cases, n_controls) %>%
  left_join(PheWAS::pheinfo , by=c("phenotype"="phecode"))

write_csv(all_results_beta, './sm_phewas_all_results_beta.csv') #Table E1 of manuscript

