# Purpose: Compare anthropometric readouts (WAZ, WLZ, LAZ, MUAC) at baseline
# between the post-SAM MAM and primary-MAM studies. Baseline refers to the 
# Measurements collected just prior to MDCF-2 and RUSF treatment, also 
# known as MAM-phase treatment.

library(tidyverse)
library(data.table)
library(stringr)

# Load primary MAM data and link anthro w stool sids
pM_anthro <- read.csv('pM_anthro.csv') %>% select(-X)

pM_key <- data.frame(intervention = c(0, 1,2,3,5,7,8,9,10),
                     stool =  c(1,2,3,5,6,7,8,6,12),
                     stool_code = c(rep('C', 7), rep('F', 2)))

for(pid_i in unique(pM_anthro$pid)) {
  if(pid_i == unique(pM_anthro$pid)[1]){
    pM_anthro_stool_key = mutate(pM_key, pid = pid_i)
  } else {
    pM_anthro_stool_key <- rbind(pM_anthro_stool_key, mutate(pM_key, pid = pid_i))
  }
}
pM_anthro_stool_key_final <- pM_anthro_stool_key %>% 
  mutate(stool.sid = paste0(pid, '11', str_pad(stool, width = 2, pad = '0'))) %>% 
  mutate(stool.sid = ifelse(intervention %in% 9:10, gsub('C', 'F', stool.sid), stool.sid)) %>% 
  left_join(pM_anthro) %>% 
  mutate(month = round(agedays/30.437)) #%>% 
#filter(stool.sid %in% rownames(tpm.mx))


# Load post-SAM MAM data and link anthro w stool sids
pSM_anthro <- read.csv('pSM_anthro.csv') %>% 
  mutate(arm = ifelse(group == 'Arm-1', 'MDCF2', 'RUSF'),
         gender = ifelse(gender == 1, 'Male', 'Female'))

pSM_key <- read.csv('pSM_SID_MasterKey.csv')

pSM_anthro_stool <- left_join(pSM_key, rename(pSM_anthro, anthro.tp = 'intervention') %>% 
                                select(-group, -area, -gender)) %>% 
  #filter(stool.sid %in% rownames(tpm.mx)) %>% 
  left_join(pSM_anthro %>% 
              filter(intervention == 1) %>% 
              select(pid, area, group, gender) %>% 
              unique()) %>% 
  mutate(month = round(agedays/30.437))


# Create dataframe anthropoemtric freatures of pSM and pM at baseline 
combined_anthro_baseline <- pSM_anthro_stool %>% 
  filter(study.wk == 0) %>% 
  select(waz, laz, wlz, arm , area) %>% 
  mutate(study = 'pSM') %>% 
  rbind(filter(pM_anthro_stool_key_final, month_fromEndOfTreatment == -3) %>% 
                     select(waz, laz, wlz, arm) %>% 
                     mutate(area = 'Dhaka',
                            study = 'pM')) %>% 
  mutate(study_area = paste0(study, "_", area))

combined_anthro_baseline %>% 
  pivot_longer(cols=waz:wlz,
               names_to = 'measurement',
               values_to = 'value') %>%
  ggplot(aes(x = study_area, y = value)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) +
  facet_wrap( ~ measurement)

combined_anthro_baseline_end <- pSM_anthro_stool %>% 
  filter(study.wk %in% c(0,12)) %>% 
  select(waz, laz, wlz, arm , area, muac, study.wk, gender) %>% 
  mutate(study = 'pSM') %>% 
  rbind(filter(pM_anthro_stool_key_final, month_fromEndOfTreatment %in% c(-3,0)) %>% 
          select(waz, laz, wlz, arm, muac, month_fromEndOfTreatment, gender) %>% 
          mutate(area = 'Dhaka',
                 study = 'pM',
                 study.wk = ifelse(month_fromEndOfTreatment == -3, 0, 12)) %>% 
          select(-month_fromEndOfTreatment)) %>% 
  mutate(study_area = paste0(study, "_", area)) 

# lm for each anthropometric outcome while controlling for study site
write.csv(summary(lm(wlz ~ study + area + gender, data = filter(combined_anthro_baseline_end, study.wk == 0)))$coef, 'lmWLZ_study_area_sex.csv')

write.csv(summary(lm(laz ~ study + area + gender, data = filter(combined_anthro_baseline_end, study.wk == 0)))$coef, 'lmLAZ_study_area_sex.csv')

write.csv(summary(lm(waz ~ study + area + gender, data = filter(combined_anthro_baseline_end, study.wk == 0)))$coef, 'lmWAZ_study_area_sex.csv')

write.csv(summary(lm(muac ~ study + area + gender, data = filter(combined_anthro_baseline_end, study.wk == 0)))$coef, 'lmMUAC_study_area_sex.csv')


