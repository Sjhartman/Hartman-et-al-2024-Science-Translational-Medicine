# loading libraries -------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(readxl)
library(factoextra)
library(edgeR)
library(vegan)
library(ggbeeswarm)
library(naturalsort)
library(scales)
library(ape)
library(fgsea)

# setting up the environment ----------------------------------------------

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC_pSM"
experiment <- "DNA_PCA"

setwd("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/230109 PCA/")

cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

# loading data ------------------------------------------------------------

# metadata
map_fecal <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "fecal_map")

# anthropometry data
dat_anth <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "anthro_dat")

# MAG abundance data
dat_tpm <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/04_postSAMMAM_tpm.csv", sep = ",", header = TRUE, row.names = 1)
dat_raw <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/04_postSAMMAM_counts.csv", sep = ",", header = TRUE, row.names = 1)

# Taxonomy data
dat_tax <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/GTDB_assignments.csv", sep = ",", header = TRUE) %>%
  mutate(species = gsub(".*? ", "", genus_specie))

# loading the list of WLZ-associated MAGs determed for the MAM phase
wlz_mags <- read.table("~/Library/CloudStorage/Box-Box/post-SAM MAM shared data/MAG abundance analysis/WLZ-Associated MAGs/test5.2_wlz~wk_ab_wkXab_area_AllSites_TxOnly.csv", sep = ",", header = TRUE) %>%
  filter(mag_abundanceXstudy_week_p_val.adj < 0.1) %>%
  arrange(desc(mag_abundanceXstudy_week))

pos_wlz_mags <- wlz_mags %>%
  filter(mag_abundanceXstudy_week > 0)

# transposing the raw MAG abundance data
dat_raw <- dat_raw %>%
  rownames_to_column("SID") %>%
  pivot_longer(cols = -SID, names_to = "MAG", values_to = "count") %>%
  pivot_wider(id_cols = "MAG", names_from = "SID", values_from = "count") %>%
  column_to_rownames("MAG")

# same with TPM data
dat_tpm <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/04_postSAMMAM_tpm.csv", sep = ",", header = TRUE, row.names = 1) %>%
  rownames_to_column("SID") %>%
  pivot_longer(-SID, names_to = "bin", values_to = "ab") %>%
  pivot_wider(names_from = SID, values_from = ab) %>%
  column_to_rownames("bin")

# filtering the abundance dataset
dat_tpm.filt <- dat_tpm[rowSums(dat_tpm > 5) > (ncol(dat_tpm) * 0.40),]
sum(dat_tpm.filt)/sum(dat_tpm) # 98% of data
nrow(dat_tpm.filt) # 613 MAGs

mags <- rownames(dat_tpm.filt)

# aligning mapping data ---------------------------------------------------

map_fecal.proc1 <- map_fecal %>%
  filter(study_phase %in% c("enrollment", "acute_rehab", "baseline")) %>%
  group_by(PID) %>%
  mutate(days_baseline = age_days - max(age_days),
         study_week = days_baseline/7) %>%
  filter(study_phase != "baseline")

map_fecal.proc2 <- map_fecal %>%
  filter(study_phase %in% c("baseline", "intervention", "follow_up_1mo", "follow_up_6mo")) %>%
  group_by(PID) %>%
  mutate(days_baseline = age_days - min(age_days),
         study_age = cut(days_baseline, 
                         breaks = c(0, 5, 15, 25, 45, 75, 105, 190, 300), 
                         #                         labels = c(0, 10, 20, 30, 60, 90, 120, 270), # this is "correct" to the design of the trial
                         labels = c(0, 15, 20, 30, 60, 90, 120, 270), # this is adapted to consider the ~10 day sampling part of week 1
                         include.lowest = TRUE),
         study_week = cut(days_baseline, 
                          breaks = c(0, 5, 15, 25, 45, 75, 105, 190, 300), 
                          labels = c(0, 2, 3, 4, 8, 12, 16, 36), # this is adapted to consider the ~10 day sampling part of week 1
                          include.lowest = TRUE),
         study_age = as.numeric(as.character(study_age)),
         study_week = as.numeric(as.character(study_week)))

map_fecal.proc <- rbind(map_fecal.proc1, map_fecal.proc2) %>%
  group_by(PID) %>%
  mutate(days_enrollment = age_days - min(age_days)) %>%
  arrange(PID, SID)

# Preparing the data for PCA ------------------------------------------------------

# VST transformation prep
map_fecal.proc.sub <- map_fecal.proc %>%
  filter(SID %in% colnames(dat_raw))

dat_raw <- dat_raw %>%
  select(map_fecal.proc.sub$SID)

# create deseq dataset
dds <- DESeqDataSetFromMatrix(countData = round(dat_raw, 0),
                              colData = map_fecal.proc.sub,
                              design = ~ 1)

# estimate library size scaling factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# perform the transformation
vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")

# extract and transpose the transformed data
dat_raw.vst <- assay(vst) %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  pivot_longer(cols = -MAG, names_to = "SID", values_to = "count") %>%
  pivot_wider(id_cols = "SID", names_from = "MAG", values_from = "count") %>%
  column_to_rownames("SID")

# Performing an alternative transformation (GeTMM) ----------------------------------------------------

# Need the lengths for each MAG assembly
dat_mag_genome <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/230109 PCA/table_s9.xlsx")

# get the RPKs
dat_rpk <- dat_raw %>%
  rownames_to_column("MAG") %>%
  left_join(dat_mag_genome %>% select(MAG_ID, bases_kb), by = c("MAG" = "MAG_ID")) %>%
  column_to_rownames("MAG") %>%
  mutate(across(colnames(dat_raw), .fns = function(x) x/bases_kb)) %>%
  select(-bases_kb)

# Create a edgeR-compatible dataset
rpk.norm <- DGEList(counts = dat_rpk,
                    group = c(rep("A",ncol(dat_rpk))))

# library size normalization factors
rpk.norm <- calcNormFactors(rpk.norm)

# finish the GeTMM transformation
dat_raw.getmm <- cpm(rpk.norm)

# transpose the dataset
dat_raw.getmm <- dat_raw.getmm %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  pivot_longer(cols = -MAG, names_to = "SID", values_to = "count") %>%
  pivot_wider(id_cols = "SID", names_from = "MAG", values_from = "count") %>%
  column_to_rownames("SID")

# convert to long form
dat_raw.getmm.long <- dat_raw.getmm %>%
  rownames_to_column("SID") %>%
  pivot_longer(-SID, names_to = "MAG", values_to = "ab") %>%
  left_join(map_fecal.proc.sub, by = "SID") %>%
  mutate(log_ab = log(ab),
         group_proc = ifelse(study_phase %in% c("enrollment", "acute_rehab"), "pre_treatment", group))

# Principal Components Analysis -------------------------------------------

# get the samples of interest
map_fecal.proc.sub.sub <- map_fecal.proc.sub %>%
  filter(study_phase %in% c("baseline", "intervention", "follow_up_1mo"))

# filter the VST-transformed dataset to match the samples for analysis
dat_raw.vst.filt <- dat_raw.vst %>%
  select(all_of(mags)) %>%
  rownames_to_column("SID") %>%
  filter(SID %in% map_fecal.proc.sub.sub$SID) %>%
  column_to_rownames("SID")

# consider scaling - how does scaling affect the VST (NOTE: I did 
# not scale in the original analysis,  following the DESeq2 tutorial)
res_pca <- prcomp(dat_raw.vst, center = TRUE, scale = FALSE)

# alternatives
#res_pca <- prcomp(dat_raw.getmm, center = TRUE, scale = FALSE)
#res_pca <- prcomp(dat_raw.vst.filt, center = TRUE, scale = FALSE)

# get the variances explained by each PC
pca_res.eig <- get_eigenvalue(res_pca) %>%
  rownames_to_column("Dim") %>%
#  slice_max(variance.percent, n = 10) %>%
  mutate(Dim = factor(Dim, levels = rev(Dim)))

# random matrix PCA to determine noise threshold
# empy matrix sized to match the other dataset
elements <- ncol(dat_raw.vst)*nrow(dat_raw.vst)

# sample and fill without replacement
dat_raw.vst.rand <- matrix(sample(c(t(dat_raw.vst)), elements, replace = FALSE), ncol = ncol(dat_raw.vst)) %>%
  as.data.frame()

# set column names
colnames(dat_raw.vst.rand) <- colnames(dat_raw.vst)

# set row names
rownames(dat_raw.vst.rand) <- rownames(dat_raw.vst)

# run PCA on the randomized dataset
res_pca.rand <- prcomp(dat_raw.vst.rand, center = TRUE, scale = TRUE)

# get the variances explained by each PC for the randomized dataset
pca_res.rand.eig <- get_eigenvalue(res_pca.rand) %>%
  rownames_to_column("Dim") %>%
#  slice_max(variance.percent, n = 10) %>%
  mutate(Dim = factor(Dim, levels = rev(Dim)))

# set the 'noise threshold' as the variance explained by the top PC from
# the randomized dataset PCA
noise_threshold <- get_eigenvalue(res_pca.rand)[1,]$variance.percent

# filter the eigenspectrum of the 'test' dataset to those explaining
# variance above the noise threshold
pca_res.eig.filt <- pca_res.eig %>%
  mutate(pass_noise_filter = ifelse(variance.percent > noise_threshold, TRUE, FALSE))

# prep for plotting
pca_res.eig.plot <- pca_res.eig %>%
  select(Dim, vp = variance.percent) %>%
#  left_join(pca_res.rand.eig %>%
#              select(Dim, vp.rand = variance.percent),
#            by = "Dim") %>%
#  pivot_longer(cols = -Dim, names_to = "type", values_to = "vp") %>%
#  mutate(type = factor(type, levels = c("vp", "vp.rand"))) %>%
  filter(Dim %in% all_of(pca_res.eig.filt$Dim[pca_res.eig.filt$pass_noise_filter])) %>%
  mutate(cumulative_vp = cumsum(vp))

# scree plot
pdf(paste0(datestring, "_", study, "_", experiment, "dna_scree_plot.pdf"))
ggplot(pca_res.eig.plot, aes(y = Dim, x = vp)) +
  geom_bar(stat = "identity", position = "dodge2") +
  geom_vline(xintercept = noise_threshold, lty = 2) +
  geom_text(aes(x = vp + 1, label = round(cumulative_vp, 1))) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0), limits = c(0, round(max(pca_res.eig.plot$vp), 1)+2), breaks = seq(0,round(max(pca_res.eig.plot$vp), 1)+2,2))
dev.off()

write.table(pca_res.eig.filt, paste0(datestring, "_", study, "_", experiment, "_PCA_eigenvalues.txt"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = F)

# extract the PC projections for each sample
pca_res.ind <- get_pca_ind(res_pca)
pca_res.ind.anno <- pca_res.ind$coord %>%
  as.data.frame() %>%
  select(Dim.1, Dim.2, Dim.3) %>%
  rownames_to_column("SID") %>%
  left_join(map_fecal.proc.sub, by = "SID")

#export data
write.table(pca_res.ind.anno, paste0(datestring, "_", study, "_", experiment, "_PCA_res_scaled.txt"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = F)

# comparing PCA results to PCoA results -----------------------------------

# this was motivated by a question about choosing euclidean versus
# other distance metrics

res_pcoa <- pcoa(vegdist(dat_raw.vst, method = "bray"))

pcoa_res.ind.anno <- res_pcoa$vectors %>%
  as.data.frame() %>%
  select(Axis.1, Axis.2, Axis.3) %>%
  rownames_to_column("SID") %>%
  left_join(map_fecal.proc.sub, by = "SID")

ggplot(pcoa_res.ind.anno, aes(x = study_week, y = Axis.1, group = paste(group, study_week), fill = paste(group))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  facet_grid(. ~ area)

# contributions to axes ---------------------------------------------------

# Which features drive variance along each PC?
pca_res.var <- get_pca_var(res_pca)
pca_res.var.anno <- pca_res.var$contrib %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  left_join(dat_tax %>% select(MAG, genus, species), by = "MAG") %>%
  mutate(label = paste(MAG, genus, species, sep = " | "),
         label = factor(label, levels = label))

write.table(pca_res.var.anno %>% select(MAG, MAG, Dim.1, Dim.2, Dim.3, genus, species), paste0(datestring, "_dna_pca_contributions.txt"), sep = "\t", row.names = FALSE)

# contributions to Dim 1
pca_res.var.anno %>% 
  select(MAG, Dim.1, genus, species) %>%
  slice_max(Dim.1, n = 10)

pca_res.var.anno %>% 
  select(MAG, Dim.1, genus, species) %>%
  arrange(desc(Dim.1)) %>%
  mutate(rank = row_number()) %>%
  filter(genus == "Prevotella")

pca_res.var.anno %>% 
  select(MAG, Dim.1, genus, species) %>%
  arrange(desc(Dim.1)) %>%
  mutate(rank = row_number()) %>%
  filter(MAG %in% c("S19C888.031_MAG", "S54C999.013_MAG"))

# contributions to Dim 2
pca_res.var.anno %>% 
  select(MAG, Dim.2, genus, species) %>%
  slice_max(Dim.2, n = 10)

pca_res.var.anno %>% 
  select(MAG, Dim.2, genus, species) %>%
  arrange(desc(Dim.2)) %>%
  mutate(rank = row_number()) %>%
  filter(genus == "Prevotella")

pca_res.var.anno %>% 
  select(MAG, Dim.2, genus, species) %>%
  arrange(desc(Dim.2)) %>%
  mutate(rank = row_number()) %>%
  filter(MAG %in% c("S19C888.031_MAG", "S54C999.013_MAG"))

# verifying directionality ------------------------------------------------

# these are spot checks for up/down abundance associated with 
# contribution to each axis
pca_res.dir <- pca_res.var$coord %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  select(MAG, Dim.1, Dim.2, Dim.3) %>%
  left_join(dat_tax %>% select(MAG, genus, species), by = "MAG")

pca_res.dir %>%
  filter(MAG %in% c("S19C888.031_MAG", "S54C999.013_MAG"))

write.table(pca_res.dir %>% select(MAG, MAG, Dim.1, Dim.2, Dim.3, genus, species), paste0(datestring, "_dna_mag_pca_coordinates.txt"), sep = "\t", row.names = FALSE)

pca_res.dir %>%
  slice_max(Dim.1, n = 10)
pca_res.dir %>%
  slice_min(Dim.1, n = 10)
pca_res.dir %>%
  slice_max(Dim.2, n = 10)
pca_res.dir %>%
  slice_min(Dim.2, n = 10)

# is there an enrichment of taxa contributing to each axis?

fgsea_res.bulk <- data.frame()

sets <- split(pca_res.dir$MAG, paste(pca_res.dir$genus, pca_res.dir$species, sep = "_"))

for (i in 1:3) {
  
  ranks <- pca_res.dir[[paste0("Dim.", i)]]
  names(ranks) <- pca_res.dir$MAG
  ranks <- sort(ranks)
  
  fgsea_res <- fgsea(pathways = sets,
                     stats = ranks,
                     eps = 0) %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))
  fgsea_res
  
  fgsea_res.bulk <- fgsea_res.bulk %>%
    bind_rows(fgsea_res %>%
                mutate(dim = paste0("Dim.", i)) %>%
                rowwise() %>%
                mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>%
                select(dim, set = pathway, padj, NES, size, leadingEdge))
  
}

fgsea_res.bulk <- fgsea_res.bulk %>%
  group_by(dim) %>%
  arrange(dim, desc(NES))

write.table(fgsea_res.bulk, paste0(datestring, "_dna_mag_enrichment.txt"), sep = "\t", row.names = FALSE)

# combining direction and contribution ------------------------------------

# this is an interpretation tool - combines the contributors and directions of
# projection for each MAG along each axis
pca_res.interp <- pca_res.var$contrib %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  select(MAG, PC1_contrib = Dim.1, PC2_contrib = Dim.2, PC3_contrib = Dim.3) %>%
  left_join(pca_res.var$coord %>%
              as.data.frame() %>%
              rownames_to_column("MAG") %>%
              select(MAG, PC1_coord = Dim.1, PC2_coord = Dim.2, PC3_coord = Dim.3),
            by = "MAG") %>%
  left_join(dat_tax %>% select(MAG, genus, species), by = "MAG") %>%
  mutate(label = paste(MAG, genus, species, sep = " | "),
         label = factor(label, levels = label))

# extract the top contributors
selector <- pca_res.interp %>%
  filter(PC1_contrib > quantile(pca_res.interp$PC1_contrib, probs = seq(0,1,0.01))[100] |
           PC2_contrib > quantile(pca_res.interp$PC2_contrib, probs = seq(0,1,0.01))[100] |
           PC3_contrib > quantile(pca_res.interp$PC3_contrib, probs = seq(0,1,0.01))[100]) %>%
  mutate(genus_species = paste(genus, species),
         dir = ifelse(PC1_coord > 0, "pos", "neg")) %>%
  group_by(dir, genus_species) %>%
  add_count()

write.table(selector, file = "230628_DNA_PCA_axis_top_contributors.txt", sep = "\t", row.names = FALSE)

# What are the abundance trajectories of these top contributors?
ggplot(dat_raw.getmm.long %>% filter(MAG %in% selector$MAG), aes(x = study_week, y = log_ab, group = PID, color = study_phase)) +
  geom_line() +
  scale_fill_manual(values = cbbPalette) +
  facet_grid(MAG ~ area) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))

# preparing for plotting these
dat_raw.getmm.long.summary <- dat_raw.getmm.long %>% 
  filter(MAG %in% selector$MAG,
         study_phase %in% c("enrollment", "baseline", "intervention"),
         study_week <= 0 | study_week %in% c(2, 12)) %>%
  mutate(study_phase_proc = ifelse(study_phase %in% c("enrollment", "baseline"), study_phase, paste0(study_phase, "_", study_week))) %>%
  pivot_wider(id_cols = c(MAG, PID, area, group), names_from = study_phase_proc, values_from = log_ab) %>%
  mutate(baseline_minus_enrollment = baseline - enrollment,
         intervention_2_minus_baseline = intervention_2 - baseline,
         intervention_12_minus_intervention_2 = intervention_12 - intervention_2) %>%
  pivot_longer(c(baseline_minus_enrollment, intervention_2_minus_baseline, intervention_12_minus_intervention_2), names_to = "comp", values_to = "log_ab") %>%
  mutate(comp = factor(comp, levels = naturalsort(unique(comp))))

ggplot(dat_raw.getmm.long.summary, aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_grid(MAG ~ area) +
  theme_bw()

ggplot(dat_raw.getmm.long.summary %>% filter(MAG %in% PC1_drivers), aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_grid(MAG ~ area, scales = "free_y") +
  theme_bw()

ggplot(dat_raw.getmm.long.summary %>% filter(MAG %in% PC2_drivers), aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_grid(MAG ~ area, scales = "free_y") +
  theme_bw()

ggplot(dat_raw.getmm.long.summary %>% filter(MAG %in% PC3_drivers), aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_grid(MAG ~ area, scales = "free_y") +
  theme_bw()

ggplot(dat_raw.getmm.long %>% filter(MAG == selector$MAG[4]), aes(x = study_week, y = log_ab, group = PID, color = group)) +
  geom_line() +
  scale_fill_manual(values = cbbPalette) +
  facet_grid(area ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))

ggplot(dat_raw.getmm.long.summary %>% filter(MAG %in% intersect(pos_wlz_mags$MAG_ID, selector$MAG)), aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_beeswarm() +
  facet_grid(MAG ~ area, scales = "free_y") +
  theme_bw()

ggplot(dat_raw.getmm.long.summary %>% filter(MAG %in% intersect(pos_wlz_mags$MAG_ID, selector$MAG)), aes(x = comp, y = log_ab)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_grid(MAG ~ area, scales = "free_y") +
  theme_bw()

# plotting via the interpretation tool table

pca_res.interp.plot <- pca_res.interp %>%
  pivot_longer(c(PC1_coord, PC2_coord, PC3_coord), names_to = "axis", values_to = "coordinate") %>%
  filter(axis == "PC1_coord" & PC1_contrib > quantile(pca_res.interp$PC1_contrib, probs = seq(0,1,0.01))[100] |
           axis == "PC2_coord" & PC2_contrib > quantile(pca_res.interp$PC2_contrib, probs = seq(0,1,0.01))[100] |
           axis == "PC3_coord" & PC3_contrib > quantile(pca_res.interp$PC3_contrib, probs = seq(0,1,0.01))[100]) %>%
  mutate(genus_species = paste(genus, species)) %>%
  arrange(axis, coordinate) %>%
  mutate(genus_species = factor(genus_species, levels = rev(unique(genus_species)))) %>%
  group_by(axis, genus_species) %>%
  add_count()

write.table(pca_res.interp.plot, file = "230628_DNA_PCA_axis_top_contributors_filtered.txt", sep = "\t", row.names = FALSE)

# looking at the drivers more specifically

pca_res.interp %>%
  slice_max(PC1_contrib, n = 10)

PC1_drivers <- pca_res.interp %>%
  filter(PC1_contrib > quantile(pca_res.interp$PC1_contrib, probs = seq(0,1,0.01))[100]) %>%
  arrange(desc(PC1_contrib)) %>%
  pull(MAG)

pca_res.interp %>%
  slice_max(PC2_contrib, n = 10)

PC2_drivers <- pca_res.interp %>%
  filter(PC2_contrib > quantile(pca_res.interp$PC2_contrib, probs = seq(0,1,0.01))[100]) %>%
  arrange(desc(PC2_contrib)) %>%
  pull(MAG)

pca_res.interp %>%
  slice_max(PC3_contrib, n = 10)

PC3_drivers <- pca_res.interp %>%
  filter(PC3_contrib > quantile(pca_res.interp$PC3_contrib, probs = seq(0,1,0.01))[100]) %>%
  arrange(desc(PC3_contrib)) %>%
  pull(MAG)

pdf(paste0(datestring, "_", study, "_", experiment, "DNA_PCA_drivers2.pdf"), height = 5, width = 8)
ggplot(pca_res.interp.plot, aes(x = coordinate, y = genus_species, fill = coordinate)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = coordinate, y = genus_species, yend = genus_species), size = 0.5) +
  geom_point(size = 3, pch = 21) +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_grid(. ~ axis) +
  theme_bw()
dev.off()

ggplot(pca_res.interp.plot, aes(x = coordinate, y = genus_species, fill = coordinate)) +
  geom_vline(xintercept = 0) +
#  geom_segment(aes(x = 0, xend = coordinate, y = genus_species, yend = genus_species)) +
  geom_beeswarm(size = 2, pch = 21, groupOnX = FALSE) +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_grid(. ~ axis) +
  theme_bw()

pdf(paste0(datestring, "_", study, "_", experiment, "DNA_PCA_drivers_filtered.pdf"), height = 5, width = 8)
ggplot(pca_res.interp.plot %>% filter(n > 1), aes(x = coordinate, y = genus_species, fill = coordinate)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = coordinate, y = genus_species, yend = genus_species)) +
  geom_point(size = 2, pch = 21) +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_grid(. ~ axis) +
  theme_bw()
dev.off()

# plotting the PCA results themselves -------------------------------------

pca_res.ind.anno <- pca_res.ind.anno %>%
  mutate(group = ifelse(study_phase %in% c("enrollment", "acute_rehab"), "pre_treatment", group))

pca_res.ind.anno.filt <- pca_res.ind.anno %>%
  group_by(PID, study_phase) %>%
  slice_max(age_days)

# PC1 versus PC2
ggplot(pca_res.ind.anno.filt, aes(x = Dim.1, y = Dim.2, color = study_phase, group = PID)) +
  geom_point() +
#  geom_line() +
  scale_color_manual(values = cbbPalette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  facet_grid(. ~ area)

# plotting projections versus time

pdf(paste0(datestring, "_", study, "_", experiment, "dna_pca_boxplots.pdf"), height = 2, width = 8)

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.1, group = paste(group, study_week), fill = paste(group))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  facet_grid(. ~ area)

dev.off()

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.2, group = paste(group, study_week), fill = paste(group))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  facet_grid(. ~ area)

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.3, group = paste(group, study_week), fill = paste(group))) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)) +
  facet_grid(. ~ area)

# trajectories

pdf(paste0(datestring, "_", study, "_", experiment, "dna_pca_lines.pdf"), height = 2, width = 8)

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.1, color = group, group = PID)) +
  theme_bw() +
  geom_path() +
  facet_grid(. ~ area) +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[3], cbbPalette[4]))

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.2, color = group, group = PID)) +
  theme_bw() +
  geom_path() +
  facet_grid(. ~ area) +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[3], cbbPalette[4]))

ggplot(pca_res.ind.anno, aes(x = study_week, y = Dim.3, color = group, group = PID)) +
  theme_bw() +
  geom_path() +
  facet_grid(. ~ area) +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[3], cbbPalette[4]))

dev.off()


# PERMANOVA ---------------------------------------------------------------

# are the sample groupings on PCA significant?

# grab the subset of data for testing - baseline
pca_res.ind.sub <- pca_res.ind.anno %>%
  filter(study_phase %in% c("baseline")) %>%
  select(SID, Dim.1, Dim.2, Dim.3, group, study_week, study_phase, area) %>%
  column_to_rownames("SID")

# run permanova (adonis) - permutations not blocked
adonis_res <- adonis2(pca_res.ind.sub %>% select(Dim.1, Dim.2, Dim.3) ~ area, data = pca_res.ind.sub, method = "euclidean")
adonis_res

# grab the subset of data for testing - acute rehab and intervention
pca_res.ind.sub <- pca_res.ind.anno %>%
  filter(study_phase %in% c("acute_rehab", "intervention")) %>%
  select(SID, Dim.1, Dim.2, Dim.3, group, study_week, study_phase) %>%
  column_to_rownames("SID")

# run permanova (adonis) - permutations not blocked
adonis_res <- adonis2(pca_res.ind.sub %>% select(Dim.1, Dim.2, Dim.3) ~ study_phase, data = pca_res.ind.sub, method = "euclidean")
adonis_res

# boxplots, but aligned and simplified
pca_res.ind.anno.long <- pca_res.ind.anno %>%
  mutate(group_proc = ifelse(study_phase %in% c("enrollment", "acute_rehab"), "pre_treatment", group),
         label = ifelse(study_phase %in% c("enrollment", "acute_rehab"), study_phase, paste(study_phase, study_week, sep = "|")),
         label = factor(label, levels = c("enrollment", "acute_rehab", "baseline|0", "intervention|2", "intervention|3", "intervention|4", "intervention|8", "intervention|12", "follow_up_1mo|16"))) %>%
  pivot_longer(c(Dim.1, Dim.2, Dim.3), names_to = "axis", values_to = "projection")

# boxplots of projections over time
pdf(paste0(datestring, "_", study, "_", experiment, "_boxplots.pdf"), height = 5, width = 8)
ggplot(pca_res.ind.anno.long, aes(x = label, y = projection, fill = group_proc, group = paste(label, group_proc))) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(axis ~ area, scales = "free_y") +
  scale_fill_manual(values = c(cbbPalette[1], cbbPalette[3], cbbPalette[4])) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))
dev.off()

# plotting centroids/summaries of each point grouping in PCA space --------

# calculate the centroids
pca_res.ind.anno.cent.psm <- pca_res.ind.anno %>%
  mutate(label = ifelse(study_phase %in% c("enrollment", "acute_rehab"), study_phase, paste(study_phase, study_week, sep = "|")),
         label = factor(label, levels = c("enrollment", "acute_rehab", "baseline|0", "intervention|2", "intervention|3", "intervention|4", "intervention|8", "intervention|12", "follow_up_1mo|16"))) %>%
  group_by(area, group, label) %>%
  summarize(PC1_mean = mean(Dim.1),
            PC1_sd = sd(Dim.1),
            PC1_se = sd(Dim.1)/length(Dim.1),
            PC2_mean = mean(Dim.2),
            PC2_sd = sd(Dim.2),
            PC2_se = sd(Dim.2)/length(Dim.2),
            PC3_mean = mean(Dim.3),
            PC3_sd = sd(Dim.3),
            PC3_se = sd(Dim.3)/length(Dim.3)) %>%
  mutate(type = "centroid")

# plot the centroids
ggplot(pca_res.ind.anno, aes(x = Dim.1, y = Dim.2, shape = group, color = study_phase)) +
  geom_point(size = 2) +
  stat_ellipse() +
  coord_equal() +
  theme_bw() +
  facet_wrap(~ area) +
  scale_color_manual(values = cbbPalette) +
  xlab(paste0("PC", rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
  ylab(paste0("PC", rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance"))

# PC1 versus PC2
pdf(paste0(datestring, "_", study, "_", experiment, "dna_centroids_pc1_pc2.pdf"))
ggplot(pca_res.ind.anno.cent.psm, aes(x = PC1_mean, y = PC2_mean, color = label, shape = group, group = group)) +
#  geom_line(arrow = arrow(ends = "last", length = unit(0.1, "inches")), color = "grey") +
  geom_point(size = 4) +
    coord_equal() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
  ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance")) +
  facet_grid(group ~ area)
dev.off()

# PC1 versus PC3
pdf(paste0(datestring, "_", study, "_", experiment, "dna_centroids_pc1_pc3.pdf"), heigh = 8, width = 8)
ggplot(pca_res.ind.anno.cent.psm %>% filter(!label %in% c("intervention|3", "intervention|8")), aes(x = PC1_mean, y = PC3_mean, color = label, shape = group, group = group)) +
  geom_path(arrow = arrow(ends = "last", length = unit(0.1, "inches")), color = "grey") +
  geom_point(size = 4) +
  coord_equal() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
  ylab(paste0(rownames(pca_res.eig[3,]), " ", round(pca_res.eig[3,]$variance.percent, 2), "% variance")) +
  facet_grid(area ~ group)
dev.off()

# PC1/2/3 individually versus study phase
ggplot(pca_res.ind.anno.cent.psm, aes(x = label, y = PC3_mean, color = label, shape = group, group = group)) +
  geom_path(arrow = arrow(ends = "last", length = unit(0.1, "inches")), color = "grey") +
  geom_errorbar(aes(ymin = PC3_mean - PC3_sd, ymax = PC3_mean + PC3_sd), color = "grey") +
  geom_point(size = 4) +
  theme_classic() +
  scale_color_manual(values = cbbPalette) +
  facet_grid(area ~ group)

# breaking out these trajectories by PC
ggplot(pca_res.ind.anno.cent.psm, aes(x = label, y = PC1_mean, color = group, group = group)) +
  geom_errorbar(aes(ymin = PC1_mean - PC1_se, ymax = PC1_mean + PC1_se)) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ area) +
  theme_classic()

ggplot(pca_res.ind.anno.cent.psm, aes(x = label, y = PC2_mean, color = group, group = group)) +
  geom_errorbar(aes(ymin = PC2_mean - PC2_se, ymax = PC2_mean + PC2_se)) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ area) +
  theme_classic()

ggplot(pca_res.ind.anno.cent.psm, aes(x = label, y = PC3_mean, color = group, group = group)) +
  geom_errorbar(aes(ymin = PC3_mean - PC3_se, ymax = PC3_mean + PC3_se)) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ area) +
  theme_classic()

# summation - PC summaries by area and experimental group
pca_res.ind.anno.cent.psm.long <- pca_res.ind.anno.cent.psm %>%
  pivot_longer(cols = -c(area, group, label, type), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("axis", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(area, group, label, type, axis), names_from = statistic, values_from = value)

ggplot(pca_res.ind.anno.cent.psm.long, aes(x = label, y = mean, color = group, group = group)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_point() +
  geom_line() +
  facet_grid(axis ~ area) +
  theme_classic()
