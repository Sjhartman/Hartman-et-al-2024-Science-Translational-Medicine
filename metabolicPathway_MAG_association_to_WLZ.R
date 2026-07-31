# Identify metabolic phenotypes correlated with MAG WLZ b3 coefficients

library(tidyverse)

pSM.mcSEED.t <- read.csv("/754MAGS_consensusBPM.csv") %>% 
  select(-c(X., MAG_species, MAG_genus)) %>% 
  column_to_rownames("genome") %>% 
  t() %>% 
  as.data.frame()

# Load & subset to WLZ associated MAGs
MAG_IDs <- read.csv("230119_pSM_MAG_IDs_pub_common_prokka.csv") %>% 
  rename('MAG_ID' = Gates_ID) %>% 
  select(MAG_ID, publicationID)

WLZ.stats <- read.csv("/test5.3_wlz~wk_ab_wkXab_area_group_AllSites_TxOnly.csv") %>% 
  left_join(MAG_IDs)

WLZ.associations <- WLZ.stats %>% 
  filter(mag_abundanceXstudy_week_p_val.adj <= 0.1)

# Subset to WLZ associated MAGs
pSM.mcSEED.t_WLZ <- pSM.mcSEED.t[,WLZ.associations$MAG_ID]
# Convert values to logical for list subsetting
pSM.mcSEED.logical <- apply(pSM.mcSEED.t_WLZ, 2, function(x) as.logical(x))
rownames(pSM.mcSEED.logical) <- rownames(pSM.mcSEED.t_WLZ)


########################################################
# Find mcSEED pathways predictive of MAG WLZ association
########################################################
# Load MAG stats
WLZmcSEED.df <- WLZ.stats %>% 
  select(MAG_ID, mag_abundanceXstudy_week, mag_abundanceXstudy_week_p_val, mag_abundanceXstudy_week_p_val.adj) %>% 
  left_join(read.csv("754MAGS_consensusBPM.csv") %>% 
              select(-X., -MAG_species, -MAG_genus) %>% 
              rename('MAG_ID'=genome),
            by = 'MAG_ID')

# filter to pathways that are present in > 5% of MAGs
freq <- colSums(WLZmcSEED.df[,5:ncol(WLZmcSEED.df)]) %>% 
  as.data.frame() %>% 
  filter(. > 613*0.05)

# Subset columns
WLZmcSEED.df <- WLZmcSEED.df[,c('MAG_ID', 'mag_abundanceXstudy_week', 'mag_abundanceXstudy_week_p_val', 'mag_abundanceXstudy_week_p_val.adj', rownames(freq))]

# Convert presence/absence into factors for the LM
WLZmcSEED.df[,5:ncol(WLZmcSEED.df)] <- apply(WLZmcSEED.df[,5:ncol(WLZmcSEED.df)], 2, function(x) factor(x, levels = c("0", "1")))

# Initialize vars
WLZpredictor_stats <- list()
WLZpredictor_models <- list()

# Run lm
for (i in 5:ncol(WLZmcSEED.df)) {
  Phenotype_i = colnames(WLZmcSEED.df)[i]
  WLZpredictor_models[[i-4]] <- summary(lm(mag_abundanceXstudy_week ~ .,WLZmcSEED.df[, c("mag_abundanceXstudy_week", colnames(WLZmcSEED.df)[i])]))
  
  WLZpredictor_stats[[i-4]] <- c(Phenotype_i, WLZpredictor_models[[i-4]]$coefficients[2,], WLZpredictor_models[[i-4]]$r.squared)
}

# Convert results to data frame
WLZpredictor_stats.df <- as.data.frame(WLZpredictor_stats)
colnames(WLZpredictor_stats.df) <- WLZpredictor_stats.df[1,]
WLZpredictor_stats.df <- WLZpredictor_stats.df[-1,]  %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Phenotype') %>% 
  left_join(as.data.frame(rowSums(pSM.mcSEED.t)/ncol(pSM.mcSEED.t)) %>% 
              rownames_to_column('Phenotype') %>% 
              rename("communityFractionalRepresentation" = `rowSums(pSM.mcSEED.t)/ncol(pSM.mcSEED.t)`))
colnames(WLZpredictor_stats.df)[2:6] <- c("Estimate", "sde", "t_value", "Pr(>|t|)", "Rsq")

# reformat and apply fdr
WLZpredictor_stats.df[,2:6] <- apply(WLZpredictor_stats.df[,2:6], 2, as.numeric)
WLZpredictor_stats.df$pval_fdr <- p.adjust(WLZpredictor_stats.df$`Pr(>|t|)`, method = 'BH')














