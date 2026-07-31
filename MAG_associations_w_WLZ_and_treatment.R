#Purpose: Identify WLZ associated MAGs in the MAM-phase

#### Libraries
library(tidyverse)
library(lmerTest)
library(progress)


######################################################
# Load datasets
######################################################
tpm.df <- read.csv("04_postSAMMAM_tpm.csv", row.names = 1)

VST_all_MAGs_fecal_info.df <- read.csv("07_VST_adjustment_counts_pSMadjusted_common_name.csv")

anthro.df <- read.csv("/AnthropometryDataPostSAMMAM_ind_bWLZ.csv")

taxa <- read.csv("Combine_GTDB_PCopriClades_Silva.csv") %>% 
  select("MAG", "genus_specie", "Clade")

MAG_IDs <- read.csv('230119_pSM_MAG_IDs_pub_common_prokka.csv')

Master_SID <- read.csv('SID_MasterKey.csv') %>% 
  filter(!is.na(study.wk))


######################################################
# Wrangle data for mixed effects model
######################################################
# Find prevalent and abunant MAGs
prev_ab_MAGs <- as.data.frame(colSums(tpm.df > 5) > 0.4*nrow(tpm.df)) %>% 
  rename(prev_ab = `colSums(tpm.df > 5) > 0.4 * nrow(tpm.df)`) %>% 
  filter(prev_ab == T)

# Filter to prevalent and abunant MAGs
VST_filt <- VST_all_MAGs_fecal_info.df %>% 
  pivot_longer(cols = S01C888.012_MAG:S99C999.v.98_sub_MAG,
               names_to = 'MAG',
               values_to = 'vst') %>% 
  rename(stool.sid = 'X') %>% 
  filter(MAG %in% rownames(prev_ab_MAGs))

# Set stool from wk 2 to wk 2.5 in master file to match up stool and anthro timepoints
for (pid in unique(Master_SID$pid)) {
  Master_SID[Master_SID$pid == pid & Master_SID$study.wk==2, 'stool.sid'] <- Master_SID[Master_SID$pid == pid & Master_SID$study.wk==2.5, 'stool.sid']
}

# Combine VST and anthro data
VST_anthro <- left_join(Master_SID %>% 
                          filter(study.wk %in% c(0,2,4,8,12)),
                        anthro.df %>% 
                            mutate(anthro.sid = paste0(pid, '0', time))) %>% 
  left_join(VST_filt)

# Set categorical vars to factors
VST_anthro$group <- factor(VST_anthro$group, levels = c("RUSF", "MDCF-2"))
VST_anthro$area <- factor(VST_anthro$area, levels = c("Kurigram", "Dhaka"))

######################################################
# Perform mixture model for WLZ-associated MAGs
######################################################
MAM_lm_res <- list()
models <- list()
MAGs <- unique(VST_anthro$MAG)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  MAG_i = MAGs[i]
  models[[MAG_i]] <- lmerTest::lmer(wlz ~ study.wk + vst + study.wk:vst + area + group + (1|pid), 
                         data = filter(VST_anthro, MAG == MAG_i))
  
  anova_i = anova(models[[MAG_i]])
  
  # save model coefficients
  MAM_lm_res[[MAG_i]] <- c(summary(models[[MAG_i]])$coefficients["study.wk",],
                            anova_i["study.wk", "Pr(>F)"],
                            summary(models[[MAG_i]])$coefficients["vst",],
                            anova_i["vst", "Pr(>F)"],
                            summary(models[[MAG_i]])$coefficients["areaDhaka",],
                            anova_i["area", "Pr(>F)"],
                           summary(models[[MAG_i]])$coefficients["groupMDCF-2",],
                           anova_i["group", "Pr(>F)"],
                           summary(models[[MAG_i]])$coefficients["study.wk:vst",],
                           anova_i["study.wk:vst", "Pr(>F)"])
}

# tabulate results
MAM_lm_res.df <- as.data.frame(MAM_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(MAM_lm_res.df) <- c("study.wk", "study.wk_sde", "study.wk_df", "study.wk_tval", "study.wk_pval", "study.wk_anova",
                             "vst", "vst_sde","vst_df","vst_tval","vst_pval","vst_anova",
                             "areaDhaka","areaDhaka_sde","areaDhaka_df","areaDhaka_tval","areaDhaka_pval","areaDhaka_anova",
                             "groupMDCF2","groupMDCF2_sde","groupMDCF2_df","groupMDCF2_tval","groupMDCF2_pval","groupMDCF2_anova",
                             "studywkXvst","studywkXvst_sde","studywkXvst_df","studywkXvst_tval","studywkXvst_pval","studywkXvst_anova")
MAM_lm_res.df <- MAM_lm_res.df %>% 
  mutate(study.wk_FDR = p.adjust(study.wk_anova, method = 'fdr'),
         vst_FDR = p.adjust(vst_anova, method = 'fdr'),
         areaDhaka_FDR = p.adjust(areaDhaka_anova, method = 'fdr'),
         groupMDCF2_FDR = p.adjust(groupMDCF2_anova, method = 'fdr'),
         studywkXvst_FDR = p.adjust(studywkXvst_anova, method = 'fdr'))

# write.csv(MAM_lm_res.df, 'MAM_phase_MAG_WLZ_associations.csv')

######################################################
# Perform mixture model for treatment associations
######################################################
MAM_lm_res <- list()
models <- list()
MAGs <- unique(VST_anthro$MAG)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  MAG_i = MAGs[i]
  models[[MAG_i]] <- lmerTest::lmer(vst ~ study.wk*group + area + (1|pid), 
                                    data = filter(VST_anthro, MAG == MAG_i))
  
  anova_i = anova(models[[MAG_i]])
  
  # save model coefficients
  MAM_lm_res[[MAG_i]] <- c(summary(models[[MAG_i]])$coefficients["study.wk",],
                           anova_i["study.wk", "Pr(>F)"],
                           summary(models[[MAG_i]])$coefficients["groupMDCF-2",],
                           anova_i["group", "Pr(>F)"],
                           summary(models[[MAG_i]])$coefficients["areaDhaka",],
                           anova_i["area", "Pr(>F)"],
                           summary(models[[MAG_i]])$coefficients["study.wk:groupMDCF-2",],
                           anova_i["study.wk:group", "Pr(>F)"])
}

# tabulate results
MAM_lm_res.df <- as.data.frame(MAM_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(MAM_lm_res.df) <- c("study.wk", "study.wk_sde", "study.wk_df", "study.wk_tval", "study.wk_pval", "study.wk_anova",
                             "groupMDCF2","groupMDCF2_sde","groupMDCF2_df","groupMDCF2_tval","groupMDCF2_pval","groupMDCF2_anova",
                             "areaDhaka","areaDhaka_sde","areaDhaka_df","areaDhaka_tval","areaDhaka_pval","areaDhaka_anova",
                             "studywkXgroupMDCF2","studywkXgroupMDCF2_sde","studywkXgroupMDCF2_df","studywkXgroupMDCF2_tval","studywkXgroupMDCF2_pval","studywkXgroupMDCF2_anova")

MAM_lm_res.df <- MAM_lm_res.df %>% 
  mutate(study.wk_FDR = p.adjust(study.wk_anova, method = 'fdr'),
         groupMDCF2_FDR = p.adjust(groupMDCF2_anova, method = 'fdr'),
         areaDhaka_FDR = p.adjust(areaDhaka_anova, method = 'fdr'),
         studywkXgroupMDCF2_FDR = p.adjust(studywkXgroupMDCF2_anova, method = 'fdr'))








