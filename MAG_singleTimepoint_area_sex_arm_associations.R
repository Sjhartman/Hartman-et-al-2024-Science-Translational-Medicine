### Create script that extracts the abundance associations
### from area, treatment arm, and sex at each timepoint

library(tidyverse)
library(data.table)
library(progress)

########################################################
# Load datasets
########################################################

# MAG IDs
MAG_IDs <- read.csv('PostSamMam/MAG_Analysis/230119_pSM_MAG_IDs_pub_common_prokka.csv') %>% 
  select(-X) %>% 
  rename(MAG='Gates_ID')

# taxa
taxa <- read.csv('PostSamMam/MAG_Analysis/Taxanomic_Assignment/GTDB_CyrusClades_SilvaCombined/Combine_GTDB_PCopriClades.csv') %>% 
  unique()

# Prev_abundant MAGs
prevAbMAGs <- read.table('PostSamMam/MAG_Analysis/Abundance_and_param_sweep/Param_sweep_passing_MAGs/MAG_passing_both_sites_allTP_0.4_5.txt',
                         header = T)

# metadata
meta <- read.csv('PostSamMam/Data/AnthropometryDataPostSAMMAM.csv') %>% 
  filter(time == 1) %>% 
  select(pid, area, group, gender) %>% 
  unique() %>% 
  left_join(read.csv('PostSamMam/Data/SID_MasterKey.csv')) %>% 
  select(pid, area, group, gender, study.wk, study.phase, stool.sid, anthro.tp) %>% 
  drop_na(stool.sid)

# Abundance
VST <- read.csv('PostSamMam/MAG_Analysis/Abundance_and_param_sweep/Abundances_w_common_names/07_VST_adjustment_counts_pSMadjusted_common_name.csv') %>% 
  rename(stool.sid = 'X')

########################################################
# Compile datasets
########################################################
compiled.dt <- data.table(MAG_IDs)[,c('publicationID', 'MAG')][as.data.table(taxa)[,c('MAG', 'genus_specie')], on = .(MAG)][as.data.table(pivot_longer(VST, cols = S01C888.012_MAG:S99C999.v.98_sub_MAG, names_to = 'MAG', values_to = 'vst')), on = .(MAG)][as.data.table(meta), on = .(stool.sid)]

compiled.dt$area <- factor(compiled.dt$area, c('Kurigram', 'Dhaka'))
compiled.dt$gender <- factor(as.character(compiled.dt$gender), c('1', '2'))
compiled.dt$group <- factor(as.character(compiled.dt$group), c('RUSF', 'MDCF-2'))

# Remove stool samples from which sequencing was not performed, and filter to prevalent and abundant MAGs
compiled_filtered.dt <- compiled.dt[!is.na(MAG), ][MAG %in% prevAbMAGs$x,]

########################################################
# Run linear models for: Enrollment (M1)
########################################################
M1_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                         data = compiled_filtered.dt[anthro.tp==1 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M1_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                           anova_i["area", "Pr(>F)"],
                           summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                           anova_i["group", "Pr(>F)"],
                           summary(models[[mag_i]])$coefficients["gender2",],
                           anova_i["gender", "Pr(>F)"])
}

# tabulate results
M1_lm_res.df <- as.data.frame(M1_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M1_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M1_lm_res.df <- mutate(M1_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')


########################################################
# Run linear models for: Baseline (M2)
########################################################
M2_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                        data = compiled_filtered.dt[study.wk==0 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M2_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                          anova_i["area", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                          anova_i["group", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["gender2",],
                          anova_i["gender", "Pr(>F)"])
}

# tabulate results
M2_lm_res.df <- as.data.frame(M2_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M2_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M2_lm_res.df <- mutate(M2_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')

########################################################
# Run linear models for: Week 4 (M3)
########################################################
M3_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                        data = compiled_filtered.dt[study.wk==4 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M3_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                          anova_i["area", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                          anova_i["group", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["gender2",],
                          anova_i["gender", "Pr(>F)"])
}

# tabulate results
M3_lm_res.df <- as.data.frame(M3_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M3_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M3_lm_res.df <- mutate(M3_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')



########################################################
# Run linear models for: Week 8 (M4)
########################################################
M4_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                        data = compiled_filtered.dt[study.wk==8 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M4_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                          anova_i["area", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                          anova_i["group", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["gender2",],
                          anova_i["gender", "Pr(>F)"])
}

# tabulate results
M4_lm_res.df <- as.data.frame(M4_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M4_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M4_lm_res.df <- mutate(M4_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')



########################################################
# Run linear models for: Week 12 (M5)
########################################################
M5_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                        data = compiled_filtered.dt[study.wk==12 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M5_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                          anova_i["area", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                          anova_i["group", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["gender2",],
                          anova_i["gender", "Pr(>F)"])
}

# tabulate results
M5_lm_res.df <- as.data.frame(M5_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M5_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M5_lm_res.df <- mutate(M5_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')



########################################################
# Run linear models for: Week 16 (M6)
########################################################
M6_lm_res <- list()
models <- list()
MAGs = unique(compiled_filtered.dt$publicationID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(MAGs))
for (i in 1:length(MAGs)) {
  pb$tick()
  
  mag_i = MAGs[i]
  models[[mag_i]] <- lm(vst ~ area + group + gender, 
                        data = compiled_filtered.dt[study.wk==16 & publicationID == mag_i,])
  anova_i = anova(models[[mag_i]])
  
  # save model coefficients
  M6_lm_res[[mag_i]] <- c(summary(models[[mag_i]])$coefficients["areaDhaka",],
                          anova_i["area", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["groupMDCF-2",],
                          anova_i["group", "Pr(>F)"],
                          summary(models[[mag_i]])$coefficients["gender2",],
                          anova_i["gender", "Pr(>F)"])
}

# tabulate results
M6_lm_res.df <- as.data.frame(M6_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(M6_lm_res.df) <- c("areaDhaka",
                            "areaDhaka_sde",
                            "areaDhaka_tvalue",
                            "areaDhaka_pval",
                            "areaDhaka_pval_anova",
                            "armMDCF2",
                            "armMDCF2_sde",
                            "armMDCF2_tvalue",
                            "armMDCF2_pval",
                            "armMDCF2_pval_anova",
                            "gender2",
                            "gender2_sde",
                            "gender2_tvalue",
                            "gender2_pval",
                            "gender2_pval_anova")

# Perform FDR
M6_lm_res.df <- mutate(M6_lm_res.df,
                       areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                       armMDCF2_anova_fdr = p.adjust(armMDCF2_pval_anova, method = 'BH'),
                       gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('publicationID')

########################################################
# Compile results/generate supplemental table
########################################################
# Combine the results
combined_results <- mutate(M1_lm_res.df, timepoint = 'enrollment') %>% 
  rbind(mutate(M2_lm_res.df, timepoint = 'baseline')) %>% 
  rbind(mutate(M3_lm_res.df, timepoint = 'week 4')) %>% 
  rbind(mutate(M4_lm_res.df, timepoint = 'week 8')) %>% 
  rbind(mutate(M5_lm_res.df, timepoint = 'week 12')) %>% 
  rbind(mutate(M6_lm_res.df, timepoint = 'week 16'))

min(combined_results$areaDhaka_anova_fdr) # 3.436831e-07
min(combined_results$armMDCF2_anova_fdr)  # 0.03091771
min(combined_results$gender2_anova_fdr)   # 0.0905779

# Filter to significant results
combined_sig_results <- combined_results %>% 
  filter(gender2_anova_fdr < 0.05 | armMDCF2_anova_fdr < 0.05 | areaDhaka_anova_fdr < 0.05)

# Format for supplemental table
supp_table <- combined_sig_results %>% 
  left_join(select(MAG_IDs, publicationID, MAG)) %>% 
  left_join(select(taxa, MAG, genus_specie)) %>% 
  select(timepoint, publicationID, genus_specie, 
         areaDhaka, areaDhaka_anova_fdr,
         armMDCF2, armMDCF2_anova_fdr,
         gender2, gender2_anova_fdr) %>% 
  rename(`MAG ID` = 'publicationID',
         `genus species` = 'genus_specie')

supp_table[,4:9] <- apply(supp_table[,4:9], 2, function(x) signif(x, digits = 3))
  
write.csv(supp_table, 'PostSamMam/MAG_Analysis/LinearMixedEffectsModels/230622_suppTable_abundance_siteEnrichment/MAG_singleTimepoint_area_sex_arm_associations.csv',
          row.names = F)


