# Plasma proteome analysis of the SAM and MAM phases

###############################################################
# Libraries
###############################################################
library(SomaDataIO)
library(tidyverse)
library(data.table)
library(progress)

###############################################################
# Datasets
###############################################################
# SomaLogic
my_adat <- read_adat("/Users/stevenhartman/Library/CloudStorage/Box-Box/000Gordon/Projects/PostSamMam/SOMALogic/WUS-216868_v4.1_EDTAPlasma_20210813/QC/QC_SS-216868_from_SS-216868_floodUnaffected_NoDropOuts/step7.SS-216868_QCd_2022-03-30.adat")

my_adat.samples <- my_adat %>% 
  filter(SampleType == "Sample") %>% 
  rownames_to_column("SampleId_SOMA") %>%
  rename("blood.sid" = SampleId)

# Create matrix of all samples in SOMAscan data
somamer.ids <- colnames(my_adat.samples)[grepl("seq.", colnames(my_adat.samples))]
SOMAscan_AllSmpl.df <- as.data.frame(my_adat.samples[,c("blood.sid", somamer.ids)])
rownames(SOMAscan_AllSmpl.df) <- NULL

SOMAscan_AllSmpl.mx <- SOMAscan_AllSmpl.df %>% 
  column_to_rownames("blood.sid") %>% 
  as.matrix()

# Log2 transform
SOMAscan_AllSmpl_transform.mx <- log2(SOMAscan_AllSmpl.mx) %>% 
  `colnames<-`(colnames(SOMAscan_AllSmpl.mx)) %>% 
  `rownames<-`(rownames(SOMAscan_AllSmpl.mx))

# Anthro df & SID key
anthro.df <- read.csv("SID_MasterKey.csv") %>% 
  dplyr::select(pid, blood.tp, anthro.sid, blood.sid, study.wk) %>% 
  full_join(read.csv("AnthropometryDataPostSAMMAM_ind_bWLZ.csv") %>% 
              mutate(anthro.sid = paste0(pid, "0", time))) %>% 
  filter(!is.na(blood.sid) | !is.na(anthro.sid))  %>% 
  mutate(study.wk = ifelse(is.na(study.wk), -1, study.wk))

# Combine SOMA transformed abundances and anthro data
soma_anthro <- left_join(as.data.frame(SOMAscan_AllSmpl_transform.mx) %>% 
                           rownames_to_column('blood.sid'),
                         anthro.df)

# Wrangle the data into a tall format before spreading to weeks in columns
soma_anthro_tall = gather(soma_anthro[,c(2:6980, 6983)],
                          key = 'soma_ID',
                          value = 'abundance',
                          -pid, 
                          -study.wk) %>% 
  unique() %>% 
  as.data.table()


soma_annotations <- attr(my_adat, "Col.Meta") %>%
  mutate(SeqId = paste0("seq.", gsub("-","\\.", SeqId))) %>%
  as.data.frame()

# Find the number of unique proteins
final_soma <- data.frame(soma = colnames(SOMAscan_AllSmpl_transform.mx))
final_soma_anno <- seq.attribute.df %>% 
  select(SeqId, EntrezGeneSymbol) %>% 
  mutate(soma = paste0('seq.', gsub('-', '.', SeqId))) %>% 
  filter(soma %in% final_soma$soma)

head(final_soma_anno)
length(unique(final_soma_anno$EntrezGeneSymbol))

###############################################################
# Calc change in plasma proteins between study phases for each participant
###############################################################
soma_anthro_tall_delta.dt <- spread(soma_anthro_tall,
                                    key = 'study.wk',
                                    value = 'abundance') %>% 
  mutate(soma_deltaSAM = `0`-`-1`,
         soma_deltaMAM_wk4 = `4` - `0`,
         soma_deltaMAM = `12` - `0`) %>% 
  as.data.table()

###############################################################
# Calc wlz rate of change per participant
###############################################################
# MAM-phase
pSM.anthro_MAM <- filter(anthro.df,
                         study.wk %in% 0:12)
wlz_slope = list()

for (pid_i in unique(pSM.anthro_MAM$pid)){
  pSM.anthro_treatment_i <- filter(pSM.anthro_MAM,
                                   pid == pid_i)
  wlz_slope[[pid_i]] <- lm(pSM.anthro_treatment_i$wlz~pSM.anthro_treatment_i$study.wk)$coefficients[2]
}
MAM_wlz_slope.df <- as.data.frame(wlz_slope) %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::rename('wlz_slope' = `pSM.anthro_treatment_i$study.wk`) %>% 
  rownames_to_column('pid')


# SAM-phase 
pSM.anthro_SAM <- filter(anthro.df,
                         study.wk %in% -1:0)

SAM_deltaWLZ.df = spread(pSM.anthro_SAM[,c('pid', 'wlz', 'study.wk')],
                         key = study.wk,
                         value = wlz) %>% 
  mutate(deltaWLZ = `0`-`-1`)



##### Combine calculations
soma_WLZ <- as.data.table(unique(filter(anthro.df[,c("pid","area","group","gender")], !is.na(area))))[as.data.table(SAM_deltaWLZ.df)[as.data.table(MAM_wlz_slope.df)[soma_anthro_tall_delta.dt, 
                                                                                                                                                                     on=.(pid)], 
                                                                                                                                     on=.(pid)], 
                                                                                                      on=.(pid)][,`i.-1`:=NULL][,`i.0`:=NULL][,`4`:=NULL][,`12`:=NULL][,`-1`:=NULL][,`0`:=NULL]
soma_WLZ$area <- factor(soma_WLZ$area, levels=c("Kurigram", "Dhaka"))
soma_WLZ$group <- factor(soma_WLZ$group, levels=c("RUSF", "MDCF-2"))
soma_WLZ$gender <- factor(soma_WLZ$gender, levels=c("1", "2"))

########################################################################
# SAM phase analysis
# Model: deltaWLZ ~ deltaSOMA + area + gender (From ALL somamers)
########################################################################
# Perform linear model
# Initialize vars
SAM_lm_res <- list()
models <- list()
SOMAmers = unique(soma_WLZ$soma_ID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(SOMAmers))
for (i in 1:length(SOMAmers)) {
  pb$tick()
  
  soma_i = SOMAmers[i]
  models[[soma_i]] <- lm(deltaWLZ ~ soma_deltaSAM + area + gender, 
                         data = filter(soma_WLZ, soma_ID == soma_i))
  
  anova_i = anova(models[[soma_i]])
  
  # save model coefficients
  SAM_lm_res[[soma_i]] <- c(summary(models[[soma_i]])$coefficients["soma_deltaSAM",],
                            anova_i["soma_deltaSAM", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["areaDhaka",],
                            anova_i["area", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["gender2",],
                            anova_i["gender", "Pr(>F)"])
}

# tabulate results
SAM_lm_res.df <- as.data.frame(SAM_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(SAM_lm_res.df) <- c("soma_deltaSAM",
                             "soma_deltaSAM_sde",
                             "soma_deltaSAM_tvalue",
                             "soma_deltaSAM_pval",
                             "soma_deltaSAM_pval_anova",
                             "areaDhaka",
                             "areaDhaka_sde",
                             "areaDhaka_tvalue",
                             "areaDhaka_pval",
                             "areaDhaka_pval_anova",
                             "gender2",
                             "gender2_sde",
                             "gender2_tvalue",
                             "gender2_pval",
                             "gender2_pval_anova")

# Perform FDR
SAM_lm_res.df <- mutate(SAM_lm_res.df,
                        soma_deltaSAM_anova_fdr = p.adjust(soma_deltaSAM_pval_anova, method = 'BH'),
                        areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                        gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH'))

sum(SAM_lm_res.df$soma_deltaSAM_anova_fdr < 0.05) # 107 total
nrow(SAM_lm_res.df[SAM_lm_res.df$soma_deltaSAM_anova_fdr < 0.05 & SAM_lm_res.df$soma_deltaSAM > 0, ]) # ) # 92 positive
nrow(SAM_lm_res.df[SAM_lm_res.df$soma_deltaSAM_anova_fdr < 0.05 & SAM_lm_res.df$soma_deltaSAM < 0, ]) # 15 negative

write.csv(SAM_lm_res.df, file = "SOMA_out/WLZ_association_SAMphase.csv")


###############################################################
# MAM phase analysis: wlz_slope ~ deltaSomamer + area + group + gender
###############################################################
# Perform linear model
# Initialize vars
MAM_lm_res <- list()
models <- list()
SOMAmers = unique(soma_WLZ$soma_ID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(SOMAmers))
for (i in 1:length(SOMAmers)) {
  pb$tick()
  
  soma_i = SOMAmers[i]
  models[[soma_i]] <- lm(wlz_slope ~ soma_deltaMAM + area + group + gender, 
                         data = filter(soma_WLZ, soma_ID == soma_i))
  
  anova_i = anova(models[[soma_i]])
  
  # save model coefficients
  MAM_lm_res[[soma_i]] <- c(summary(models[[soma_i]])$coefficients["soma_deltaMAM",],
                            anova_i["soma_deltaMAM", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["areaDhaka",],
                            anova_i["area", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["groupMDCF-2",],
                            anova_i["group", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["gender2",],
                            anova_i["gender", "Pr(>F)"])
}

# tabulate results
MAM_lm_res.df <- as.data.frame(MAM_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(MAM_lm_res.df) <- c("soma_deltaMAM",
                             "soma_deltaMAM_sde",
                             "soma_deltaMAM_tvalue",
                             "soma_deltaMAM_pval",
                             "soma_deltaMAM_pval_anova",
                             "areaDhaka",
                             "areaDhaka_sde",
                             "areaDhaka_tvalue",
                             "areaDhaka_pval",
                             "areaDhaka_pval_anova",
                             "groupMDCF2",
                             "groupMDCF2_sde",
                             "groupMDCF2_tvalue",
                             "groupMDCF2_pval",
                             "groupMDCF2_pval_anova",
                             "gender2",
                             "gender2_sde",
                             "gender2_tvalue",
                             "gender2_pval",
                             "gender2_pval_anova")

# Perform FDR
MAM_lm_res.df <- mutate(MAM_lm_res.df,
                        soma_deltaMAM_anova_fdr = p.adjust(soma_deltaMAM_pval_anova, method = 'BH'),
                        areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                        groupMDCF2_anova_fdr = p.adjust(groupMDCF2_pval_anova, method = 'BH'),
                        gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('SeqId') %>% 
  mutate(soma_annotations %>% 
           select(SeqId, SeqIdVersion, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol))


MAM_lm_res_unique.df <- MAM_lm_res.df %>% 
  arrange(EntrezGeneSymbol, soma_deltaMAM_pval_anova) %>% 
  filter(duplicated(EntrezGeneSymbol) == FALSE)

write.csv(MAM_lm_res.df, 'SOMA_out/WLZ_association_MAMphase.csv')


#######################################################################
# Combine pM, MAM-phase, and SAM-phase results into one table
#######################################################################

### Load SAM WLZ-associations
SAM <- read.csv("SOMA_out/WLZ_association_SAMphase.csv") %>% 
  rename(SeqId = 'X') %>% 
  select(SeqId, soma_deltaSAM, soma_deltaSAM_sde, soma_deltaSAM_anova_fdr)

### Load MAM WLZ-associations
MAM <- read.csv('SOMA_out/WLZ_association_MAMphase.csv') %>% 
  select(-X) %>% 
  select(SeqId, SeqIdVersion, TargetFullName, EntrezGeneSymbol, EntrezGeneID, soma_deltaMAM, soma_deltaMAM_sde, soma_deltaMAM_anova_fdr)

### Load pM WLZ-associations
# Note - this is just the WLZ assocations from Chen et al. 2012 supplementary
pM <- read.csv('ProteinAbundanceVsPonderalGrowth.csv') %>% 
  rename(EntrezGeneSymbol = 'Entrez.Gene.Symbol',
         pM_q.value = 'q.value',
         pM_pearson = 'Pearson.Rho') %>% 
  arrange(EntrezGeneSymbol, pM_q.value) %>% 
  filter(duplicated(EntrezGeneSymbol) == FALSE) %>% 
  select(EntrezGeneSymbol, pM_pearson, pM_q.value)

# Combine results
MAM_SAM_pM_compare <- full_join(MAM, SAM) %>% 
  left_join(pM) %>% 
  mutate(MAM_wlz_association = ifelse(soma_deltaMAM_anova_fdr > 0.05, 'nonsig', ifelse(soma_deltaMAM > 0, "pos", "neg")),
         SAM_wlz_association = ifelse(soma_deltaSAM_anova_fdr > 0.05, 'nonsig', ifelse(soma_deltaSAM > 0, "pos", "neg")),
         pM_wlz_association = ifelse(pM_q.value > 0.05, 'nonsig', ifelse(pM_pearson > 0, "pos", "neg"))) %>% 
  arrange(desc(soma_deltaMAM)) %>% 
  mutate(soma_deltaMAM = substr(formatC(signif(soma_deltaMAM,digits=3), digits=3,format="fg", flag="#"), 1, 7),
         soma_deltaMAM_sde = formatC(signif(soma_deltaMAM_sde,digits=3), digits=3,format="fg", flag="#"),
         soma_deltaMAM_anova_fdr = formatC(signif(soma_deltaMAM_anova_fdr,digits=3), digits=3,format="fg", flag="#"),
         soma_deltaSAM = substr(formatC(signif(soma_deltaSAM,digits=3), digits=3,format="fg", flag="#"), 1, 7),
         soma_deltaSAM_sde = formatC(signif(soma_deltaSAM_sde,digits=3), digits=3,format="fg", flag="#"),
         soma_deltaSAM_anova_fdr = formatC(signif(soma_deltaSAM_anova_fdr,digits=3), digits=3,format="fg", flag="#"),
         pM_pearson = formatC(signif(pM_pearson,digits=3), digits=3,format="fg", flag="#"),
         pM_q.value = formatC(signif(pM_q.value,digits=3), digits=3,format="fg", flag="#"))

write.csv(MAM_SAM_pM_compare, 'SOMA_out/WLZ_association_supplemental_table.csv', row.names = F)




########################################################################
# Treatment arm associations:
# Model: deltaSoma ~ group + area + gender
########################################################################
# Perform linear model
# Initialize vars
soma_WLZ_MAM <- soma_WLZ %>% 
  na.omit() %>% 
  left_join(soma_annotations %>% 
              rename('soma_ID' = SeqId) %>% 
              select(soma_ID, Target, TargetFullName))

MAM_lm_res <- list()
models <- list()
SOMAmers = unique(soma_WLZ_MAM$soma_ID)

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(SOMAmers))
for (i in 1:length(SOMAmers)) {
  pb$tick()
  
  soma_i = SOMAmers[i]
  models[[soma_i]] <- lm(soma_deltaMAM ~ area + group + gender, 
                         data = filter(soma_WLZ_MAM, soma_ID == soma_i))
  
  anova_i = anova(models[[soma_i]])
  
  # save model coefficients
  MAM_lm_res[[soma_i]] <- c(summary(models[[soma_i]])$coefficients["areaDhaka",],
                            anova_i["area", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["groupMDCF-2",],
                            anova_i["group", "Pr(>F)"],
                            summary(models[[soma_i]])$coefficients["gender2",],
                            anova_i["gender", "Pr(>F)"])
}

# tabulate results
MAM_lm_res.df <- as.data.frame(MAM_lm_res) %>% 
  t() %>% 
  as.data.frame()
colnames(MAM_lm_res.df) <- c("areaDhaka",
                             "areaDhaka_sde",
                             "areaDhaka_tvalue",
                             "areaDhaka_pval",
                             "areaDhaka_pval_anova",
                             "groupMDCF2",
                             "groupMDCF2_sde",
                             "groupMDCF2_tvalue",
                             "groupMDCF2_pval",
                             "groupMDCF2_pval_anova",
                             "gender2",
                             "gender2_sde",
                             "gender2_tvalue",
                             "gender2_pval",
                             "gender2_pval_anova")

# Perform FDR
MAM_lm_res.df <- mutate(MAM_lm_res.df,
                        areaDhaka_anova_fdr = p.adjust(areaDhaka_pval_anova, method = 'BH'),
                        groupMDCF2_anova_fdr = p.adjust(groupMDCF2_pval_anova, method = 'BH'),
                        gender2_anova_fdr = p.adjust(gender2_pval_anova, method = 'BH')) %>% 
  rownames_to_column('SeqId') %>% 
  left_join(soma_annotations %>% 
              select(SeqId, SeqIdVersion, EntrezGeneSymbol, TargetFullName)) %>% 
  select(SeqId, SeqIdVersion, EntrezGeneSymbol, TargetFullName, 
         groupMDCF2, groupMDCF2_sde, groupMDCF2_pval_anova, groupMDCF2_anova_fdr,
         areaDhaka, areaDhaka_sde, areaDhaka_pval_anova, areaDhaka_anova_fdr,
         gender2, gender2_sde, gender2_pval_anova, gender2_anova_fdr) %>% 
  mutate(groupMDCF2 = formatC(signif(groupMDCF2,digits=3), digits=3,format="fg", flag="#"),
         groupMDCF2_sde = formatC(signif(groupMDCF2_sde,digits=3), digits=3,format="fg", flag="#"),
         groupMDCF2_pval_anova = formatC(signif(groupMDCF2_pval_anova,digits=3), digits=3,format="fg", flag="#"),
         groupMDCF2_anova_fdr = formatC(signif(groupMDCF2_anova_fdr,digits=3), digits=3,format="fg", flag="#"),
         areaDhaka = formatC(signif(areaDhaka,digits=3), digits=3,format="fg", flag="#"),
         areaDhaka_sde = formatC(signif(areaDhaka_sde,digits=3), digits=3,format="fg", flag="#"),
         areaDhaka_pval_anova = formatC(signif(areaDhaka_pval_anova,digits=3), digits=3,format="fg", flag="#"),
         areaDhaka_anova_fdr = formatC(signif(areaDhaka_anova_fdr,digits=3), digits=3,format="fg", flag="#"),
         gender2 = formatC(signif(gender2,digits=3), digits=3,format="fg", flag="#"),
         gender2_sde = formatC(signif(gender2_sde,digits=3), digits=3,format="fg", flag="#"),
         gender2_pval_anova = formatC(signif(gender2_pval_anova,digits=3), digits=3,format="fg", flag="#"),
         gender2_anova_fdr = formatC(signif(gender2_anova_fdr,digits=3), digits=3,format="fg", flag="#"))

write.csv(MAM_lm_res.df, 'SOMA_out/txArm_association_MAMphase.csv')

#####################################################
# GSEA of treatment arm association
#####################################################
library(fgsea)

# Create a list of positively WLZ associated somamers during the MAM-phase
wlz_somamers <- list(pos_wlz_somamers = c(MAM_SAM_pM_compare[MAM_SAM_pM_compare$soma_deltaMAM > 0 & MAM_SAM_pM_compare$soma_deltaMAM_anova_fdr < 0.05, ]$SeqId),
                     neg_wlz_somamers = c(MAM_SAM_pM_compare[MAM_SAM_pM_compare$soma_deltaMAM < 0 & MAM_SAM_pM_compare$soma_deltaMAM_anova_fdr < 0.05, ]$SeqId))

# Create a ranked list of all modeled somamers
ranked_df <- MAM_lm_res.df %>% 
  mutate(rank = sign(groupMDCF2)*-log2(groupMDCF2_pval)) %>% 
  arrange(desc(rank))
ranked_vec = ranked_df$rank
names(ranked_vec) <- ranked_df$SeqId

# Run GSEA
fgsea_res <- fgsea(pathways = wlz_somamers,
                   stats = ranked_vec)



