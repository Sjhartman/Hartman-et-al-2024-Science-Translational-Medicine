# Anthropemtric response to MDCF-2 and RUSF treatment

# Load packages
library(tidyverse)
library(lmerTest)

# Load anthropometric dataset
anthro <- read.csv('AnthropometryDataPostSAMMAM.csv') %>% 
  mutate(anthro.sid = paste0(pid, '0', time)) %>% 
  left_join(read.csv('SID_MasterKey.csv')) %>% 
  filter(dropout_status == 'non_dropout') %>% 
  mutate(group = factor(group, levels = c('RUSF', 'MDCF-2')),
         area = factor(area, levels = c('Kurigram', 'Dhaka')),
         gender = factor(gender, levels = c('1', '2')))


### Anthropometric responses in flood-unaffected participants during the MDCF-2 and RUSF treatment phase
# WLZ
summary(lmer(wlz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))

# WAZ
summary(lmer(waz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))

# LAZ
summary(lmer(laz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))

# MUAC
summary(lmer(muac ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))

# weight
summary(lmer(weight ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))

# length
summary(lmer(length ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:12,
                           flood_status == 'unaffected')))


### Anthropometric responses in flood-unaffected participants during the MDCF-2 and RUSF treatment phase 
# and 1-month follow-up
# WLZ
summary(lmer(wlz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))

# WAZ
summary(lmer(waz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))

# LAZ
summary(lmer(laz ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))

# MUAC
summary(lmer(muac ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))

# weight
summary(lmer(weight ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))

# length
summary(lmer(length ~ study.wk*group + study.wk*area + gender + (1|pid),
             data = filter(anthro, 
                           study.wk %in% 0:16,
                           flood_status == 'unaffected')))
