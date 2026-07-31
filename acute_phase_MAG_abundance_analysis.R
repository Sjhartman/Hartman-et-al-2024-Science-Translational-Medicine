# load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(vegan)
#library(reshape2)
library(variancePartition)
library(edgeR)
library(fgsea)
library(DESeq2)
library(data.table)
library(lmerTest)
library(progress)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(UpSetR)

# defining questions ------------------------------------------------------

# Questions:
# 1. What are the characteristics of a SAM community at enrollment? Average abundance at enrollment - dif between sites? Commonalities? Species-level differences?
# 2. What taxa change with acute intervention? Is there an overall difference in abundance with ABX treatment/nutrition? Which specific taxa? By area? Is it different? Commonalities?
# 3. What are the characteristics of a post-SAM MAM community? (see 1) by area? commonalities?
# 4. Is there a relationship between time in treatment and later WLZ response to MDCF-2 or RUSF? Effect of diversity? Relationship between configuration and time required to acute rehab?

# time of acute treatment (baseline - enrollment)?
# starting WLZ at enrollment?
# staring WLZ at baseline?

# Thought: diversity OR PCA to look at distribution/relationships of configurations?
# Thought: look for ABX-R elements in taxa that increase with acute rehab

# setting up environment --------------------------------------------------

setwd("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/230518 wrapping up SAM phase MAG analysis/")

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC_pSM"

cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

# loading data ------------------------------------------------------------

dat_anth <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "anthro_dat")
dat_tpm <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/04_postSAMMAM_tpm.csv", sep = ",", header = TRUE, row.names = 1)
dat_tpm.pm <- read.table("/Users/hibberdm/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_tpm_full_set.txt", sep = "\t", header = TRUE, row.names = 1) %>%
  rownames_to_column("MAG_ID") %>%
  pivot_longer(-MAG_ID, names_to = "SID", values_to = "tpm") %>%
  pivot_wider(id_cols = SID, names_from = MAG_ID, values_from = tpm) %>%
  column_to_rownames("SID")
dat_raw <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/04_postSAMMAM_counts.csv", sep = ",", header = TRUE, row.names = 1)
tax <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/210914_MAG_abundance_data/GTDB_assignments.csv", sep = ",", header = TRUE)

mag_name_key <- read.table("~/Library/CloudStorage/Box-Box/post-SAM MAM shared data/Reports/Manuscript/versions/230522/230523 pSM MAG name key.txt", header = TRUE, sep = "\t")

map_fecal <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "fecal_map")
map_fecal.pm <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")

anth_key <- data.frame(time = c(1:9),
                       study_age = c(0, 0, 15, 30, 45, 60, 75, 90, 120),
                       study_week = c(0, 0, 2, 4, 6, 8, 10, 12, 16),
                       study_phase = c("enrollment", "baseline", rep("intervention", 6), "follow_up_1mo"))

dat_anth.anno <- dat_anth %>%
  left_join(anth_key, by = "time") %>%
  filter(dropout_status == "non_dropout")

dat_anth.anno.acute <- dat_anth.anno %>%
  filter(study_phase %in% c("enrollment", "baseline") & flood_status_jan == "unaffected") %>%
  select(-flood_status_original, -flood_status_jan, -flood_status_may, -dropout_status, -study_week) %>%
  pivot_wider(id_cols = c(pid, area, group, gender), names_from = study_phase, values_from = c(wlz, waz, laz, muac, agedays)) %>%
  mutate(acute_rehab_days = agedays_baseline - agedays_enrollment,
         acute_rehab_bin = cut(acute_rehab_days, breaks = seq(6, 20, 2)))

# these data were calculated within the "anthropometry analysis.R" notebook
acute_recovery_metrics <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220107 anthropometry re-analysis/220201_1244_MDCF_POC_pSM_enrollment_anthro_vs_rate_of_change.txt", sep = "\t", header = TRUE) %>%
  left_join(dat_anth.anno.acute %>%
              select(-area, -group, -gender),
            by = "pid")

# formatting/aligning the mapping file
# this is necessary to convert the sampling timepoint labels, which are sequential, to their actual timing (in days/weeks)
map_fecal.proc1 <- map_fecal %>%
  filter(study_phase %in% c("enrollment", "acute_rehab")) %>%
  group_by(PID)

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
  mutate(days_enrollment = age_days - min(age_days)) %>%
  arrange(PID, SID)

all_participants <- unique(map_fecal.proc$PID)

# Exploring filtering criteria between the current and former 
# studies of MAG abundance relationships with anthropometry

# post-SAM MAM (current study)
dat_tpm.filt <- dat_tpm
dat_tpm.filt[dat_tpm.filt < 5] <- 0

dat_tpm.filt2 <- dat_tpm[colSums(dat_tpm > 5) > (nrow(dat_tpm) * 0.40)]
sum(dat_tpm.filt2)/sum(dat_tpm)
ncol(dat_tpm.filt2) #613 MAGs

# primary MAM (Hibberd, Webber, et al)
dat_tpm.pm.filt <- dat_tpm.pm
dat_tpm.pm.filt[dat_tpm.pm.filt < 5] <- 0

dat_tpm.pm.filt <- dat_tpm.pm[colSums(dat_tpm.pm > 5) > (nrow(dat_tpm.pm) * 0.40)]

# Microbial diversity analyses --------------------------------------------

# calculate beta diversity metrics using vegan
bray_df <- dat_tpm.filt %>%
  vegdist(method = "bray", diag = TRUE) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("sample_1") %>%
  pivot_longer(-sample_1, names_to = "sample_2", values_to = "bray") %>%
  mutate(pid_1 = substr(sample_1, 0, 7),
         pid_2 = substr(sample_2, 0, 7)) %>%
  filter(pid_1 == pid_2) %>%
  arrange(sample_1) %>%
  group_by(pid_1) %>%
  filter(pid_1 == first(pid_1))

# annotate the table with metadata
bray_df.anno <- bray_df %>%
  left_join(map_fecal.proc %>% select(SID, days_enrollment), by = c("sample_1" = "SID")) %>%
  filter(days_enrollment == 0) %>%
  select(-days_enrollment, -PID) %>%
  left_join(map_fecal.proc, by = c("sample_2" = "SID")) %>%
  mutate(bray_dist = bray,
         study_phase = factor(study_phase, levels = c("enrollment", "acute_rehab", "baseline", "intervention", "follow_up_1mo"))) %>%
  arrange(PID, days_enrollment) %>%
  distinct() %>%
  left_join(acute_recovery_metrics %>%
              select(pid, wlz, wlz_severity, wlz_severity_q, acute_rehab_bin, acute_rehab_days), by = c("PID" = "pid")) %>%
  distinct()

write.table(bray_df.anno, paste(datestring, study, "bray_div_data.txt", sep = "_"), sep = "\t", row.names = FALSE)

# get the participants for which we have microbiota data
participants <- bray_df.anno %>%
  select(PID, area, group) %>%
  arrange(area, group, PID) %>%
  distinct() %>%
  pull(PID)

# filter the dataset to the acute rehabilitation phase
bray_df.anno.acute <- bray_df.anno %>%
  filter(study_phase %in% c("enrollment", "acute_rehab"))

# plotting beta diversity trajectories
pdf(paste(datestring, study, "bray_div_acute_phase_trajectories.pdf", sep = "_"), width = 11, height = 8.5)
p <- ggplot(bray_df.anno.acute, aes(x = days_enrollment, y = bray_dist)) +
  geom_line(aes(group = PID)) +
  geom_point(aes(color = study_phase), size = 2) +
  geom_smooth(method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
  #  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p)
p <- ggplot(bray_df.anno.acute, aes(x = days_enrollment, y = bray_dist)) +
  geom_line(aes(group = PID)) +
  geom_point(aes(color = study_phase), size = 2) +
  geom_smooth(method = "loess") +
  scale_y_continuous(limits = c(0, 1)) +
#  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(acute_rehab_bin ~ area)
print(p)

p2 <- ggplot(bray_df.anno.acute, aes(x = days_enrollment, y = bray_dist, group = days_enrollment)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1)) +
#  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p2)
dev.off()

pdf(paste(datestring, study, "bray_div_acute_phase_enrollment_vs_baseline.pdf", sep = "_"), width = 11, height = 8.5)
p <- ggplot(bray_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = bray_dist)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1)) +
  #  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p)

p2 <- ggplot(bray_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = bray_dist, group = PID)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  #  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p2)

p3 <- ggplot(bray_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = bray_dist, group = PID)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  #  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(wlz_severity_q ~ area)
print(p3)
dev.off()

# alpha diversity calculations

shan_df.anno <- diversity(dat_tpm.filt, index = "shannon", MARGIN = 1, base = exp(1)) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("SID") %>%
  transmute(SID = SID, shannon = V1) %>%
  left_join(map_fecal.proc %>% select(SID, days_enrollment, study_phase, group, area), by = "SID") %>%
  mutate(study_phase = factor(study_phase, levels = c("enrollment", "acute_rehab", "baseline", "intervention", "follow_up_1mo"))) %>%
  arrange(PID, days_enrollment) %>%
  left_join(acute_recovery_metrics %>%
              select(pid, wlz, wlz_severity, wlz_severity_q), by = c("PID" = "pid")) %>%
  distinct()

write.table(shan_df.anno, paste(datestring, study, "shannon_div_data.txt", sep = "_"), sep = "\t", row.names = FALSE)

# filter to acute rehabilitation phase
shan_df.anno.acute <- shan_df.anno %>%
  filter(study_phase %in% c("enrollment", "acute_rehab"))

pdf(paste(datestring, study, "shannon_div_acute_phase_enrollment_vs_baseline.pdf", sep = "_"), width = 11, height = 8.5)
p3 <- ggplot(shan_df.anno.acute, aes(x = days_enrollment, y = shannon, group = PID)) +
  geom_line() +
  geom_point(aes(color = study_phase), size = 2) +
  #  scale_y_continuous(limits = c(0, 1)) +
#  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p3)

p4 <- ggplot(shan_df.anno %>% filter(study_phase %in% c("enrollment", "acute_rehab", "baseline")), aes(x = days_enrollment, y = shannon, group = study_phase, fill = study_phase)) +
  geom_boxplot() +
  #  scale_y_continuous(limits = c(0, 1)) +
  #  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  facet_grid(. ~ area)
print(p4)

p <- ggplot(shan_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = shannon, group = PID)) +
  geom_line() +
  theme_classic() +
  facet_grid(. ~ area)
print(p)

p <- ggplot(shan_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = shannon, group = study_phase, fill = study_phase)) +
  geom_boxplot() +
  theme_classic() +
  facet_grid(. ~ area)
print(p)

p <- ggplot(shan_df.anno %>% filter(study_phase %in% c("enrollment", "baseline")), aes(x = study_phase, y = shannon, group = study_phase, fill = study_phase)) +
  geom_boxplot() +
  theme_classic() +
  facet_grid(area ~ wlz_severity_q)
print(p)

dev.off()

# testing relationships between diversity and study phase in the acute rehabilitation phase

# alpha
test_dat <- shan_df.anno %>% 
  filter(study_phase %in% c("enrollment", "acute_rehab")) %>%
  group_by(PID) %>%
  filter(row_number() %in% c(1, n()))

mod2 <- lmer(shannon ~ study_phase + area + days_enrollment + (1|PID), data = test_dat)
summary(mod2)

t.test(shannon ~ study_phase, data = test_dat)

# beta
test_dat <- bray_df.anno %>% 
  filter(study_phase %in% c("enrollment", "acute_rehab")) %>%
  group_by(PID)

mod2 <- lmer(bray ~ study_phase + area + days_enrollment + (1|PID), data = test_dat)
mod2 <- lm(bray ~ study_phase + area, data = test_dat)
summary(mod2)

# Identifying MAGs that increase/decrease significantly during  --------
# acute rehabilitation

# get metadata for the acute rehabiliation phase
metadata <- map_fecal.proc %>%
  filter(SID %in% rownames(dat_raw),
         study_phase %in% c("enrollment", "acute_rehab", "baseline")) %>%
  column_to_rownames("SID") %>%
  mutate(group = factor(group, levels = c("RUSF", "MDCF-2")))

# align the abundance dataset with the selected metadata
dat_raw.t.filt <- dat_raw %>%
  rownames_to_column("SID") %>%
  pivot_longer(-SID, names_to = "MAG") %>% 
  pivot_wider(names_from=SID, values_from=value) %>%
  select(c(MAG, rownames(metadata))) %>%
  filter(MAG %in% colnames(dat_tpm.filt2)) %>%
  column_to_rownames("MAG")

# prepare a limma/voom/edgeR/dream-compatible data structure
# edgeR doesn't handle mixed effects - dream does but the overall approach is different
dge <- DGEList(counts = dat_raw.t.filt,
               samples = metadata)
# calculate library size normalization factors
dge <- calcNormFactors(dge)

# define the model
form <- ~ area + days_enrollment + (1|PID)

# estimate precision weights (mean-variance relationships) for each MAG. This approach builds on the 
# 'traditional' limma/voom approach to allow mixed effects. 
vobjDream <- voomWithDreamWeights(counts = dge, 
                                  formula = form, 
                                  data = metadata,
                                  span = 0.8, 
                                  plot = TRUE,
                                  save.plot = TRUE)

# Fit the models for each MAG, then perform hypothesis tests as defined in the contrast matrix (default or custom)
# degrees of freedom are estimated using the Satterthwaite approximation by default (Kenward-Roger is more accurate, but slower - recommended for < 10 samples)
fitmm <- dream(exprObj = vobjDream,
               formula = form, 
               data = metadata)

# exploring the modeling objects and contrasts - this can help you identify/specify your calls for results
#L = getContrast(vobjDream, form, metadata, c('study_armX'))
#plotContrasts(L)
# Examine design matrix
head(fitmm$design, 3)

# Extract results of DE testing for coefficients of interest
res.days_enrollment <- topTable(fitmm, coef="days_enrollment", adjust.method = "BH", n = Inf) %>%
  arrange(z.std)

res.days_enrollment.anno <- res.days_enrollment %>%
  filter(adj.P.Val < 0.05) %>%
  rownames_to_column("MAG") %>%
  left_join(tax %>%
              select(MAG, genus, genus_species = genus_specie),
            by = "MAG")

# There are a lot of MAGs with significant changes here. I need to find a way to break this down to show/summarize

# what are the most common genera in this list of MAGs
top_genera <- res.days_enrollment.anno %>%
  group_by(genus, sign(logFC)) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  filter(n > 10)

# plot the top genera results
ggplot(res.days_enrollment.anno, aes(x = logFC, y = -1*log10(adj.P.Val))) +
  geom_point(data = res.days_enrollment.anno %>% filter(!genus %in% top_genera$genus)) +
  geom_point(aes(color = genus), data = res.days_enrollment.anno %>% filter(genus %in% top_genera$genus), size = 3) +
  scale_color_manual(values = cbbPalette) +
  theme_bw()

# top species
top_species <- res.days_enrollment.anno %>%
  filter(genus_species != "") %>%
  group_by(genus_species, sign(logFC)) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  filter(n > 5)

# top positive species
top_species.up <- res.days_enrollment.anno %>%
  filter(genus_species != "" & logFC > 0) %>%
  group_by(genus_species) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# top mags
top_mags <- res.days_enrollment.anno %>%
  filter(abs(logFC) > quantile(abs(logFC), probs = seq(0,1,0.05))[20])

# more plotting
ggplot(res.days_enrollment.anno, aes(x = logFC, y = -1*log10(adj.P.Val))) +
  geom_point(data = res.days_enrollment.anno %>% filter(!MAG %in% top_mags$MAG)) +
  geom_point(aes(color = genus_species), data = res.days_enrollment.anno %>% filter(MAG %in% top_mags$MAG), size = 3) +
  geom_label_repel(aes(label = genus_species), data = res.days_enrollment.anno %>% filter(MAG %in% top_mags$MAG)) +
  scale_color_manual(values = cbbPalette) +
  theme_bw()

ggplot(res.days_enrollment.anno %>% filter(genus_species %in% top_species$genus_species), aes(x = logFC, y = genus_species)) + 
  geom_density_ridges() +
  geom_vline(xintercept = 0, color = "blue", lty = 2)

ggplot(res.days_enrollment.anno %>% filter(genus %in% top_genera$genus), aes(x = logFC, y = genus)) + 
  geom_density_ridges() +
  geom_vline(xintercept = 0, color = "blue", lty = 2)

# top positive MAGs
top_mags.up <- res.days_enrollment.anno %>%
  filter(logFC > 0 & logFC > quantile(logFC, probs = seq(0,1,0.05))[20])

# Looking at ABX-R - does ABX resistance drive positive responses?  --------

# this data comes from an aggregation of AMRFinderPlus results, run on all MAGs
dat_amr <- read.table("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/230221_amr_profiling/230317_1225_MDCF_POC_pSM_MAG_AMR_data.txt", header = TRUE, sep = "\t", quote = "") %>%
  mutate(MAG = gsub("$", "_MAG", MAG_ID))

# looking at subsets of the bulk AMR dataset - here, top positive mags
dat_amr.filt <- dat_amr %>%
  filter(MAG %in% top_mags.up$MAG)

# looking at each top mag individually
dat_amr.filt %>%
  filter(MAG == top_mags.up$MAG[1] & type == "AMR")

# shifting to looking at the whole spectrum of responses via enrichment analyses
dat_amr.amr <- dat_amr %>%
  filter(type == "AMR") %>%
  select(MAG, class) %>%
  distinct()

# prepare MAG sets
amr_sets <- split(dat_amr.amr$MAG, dat_amr.amr$class)

# rankings for GSEA are 
ranks <- res.days_enrollment$logFC
names(ranks) <- rownames(res.days_enrollment)

# run gsea - are sets of similar AMR markers enriched in
# organisms whose abundance increases during acute rehabiliation
fgsea_res.amr <- fgsea(amr_sets,
                           minSize = 5,
                           maxSize = length(ranks)/2,
                           stats = ranks,
                           eps = 0) %>%
  arrange(desc(NES))

# filter results
fgsea_res.amr.sig <- fgsea_res.amr %>%
  filter(padj < 0.05 & pathway != "")

# convert to a form easier to write to file
fgsea_res.amr.sig.format <- fgsea_res.amr.sig %>%
  mutate(leadingEdge = NA)
for (i in 1:nrow(fgsea_res.amr.sig)) {
  LE <- data.frame(MAG_ID = fgsea_res.amr.sig$leadingEdge[i][[1]]) %>%
    left_join(mag_name_key, by = "MAG_ID")
  fgsea_res.amr.sig.format$leadingEdge[i] <- paste(sort(LE$MAG), collapse = ";")
}

write.table(fgsea_res.amr.sig.format, file = "fgsea_res_sig_amr.txt", sep = "\t", row.names = FALSE)

# are there common MAGs that drive these enrichments?
amr_sig_leading_edges <- fgsea_res.amr.sig$leadingEdge
names(amr_sig_leading_edges) <- fgsea_res.amr.sig$pathway

dat_amr.upset <- data.frame(MAG = character(0))
for (i in 1:length(amr_sig_leading_edges)) {
  name <- names(amr_sig_leading_edges)[i]
  dat_amr.upset.sub <- data.frame(MAG = amr_sig_leading_edges[i][[1]], 
                                  path = 1)
  colnames(dat_amr.upset.sub)[2] <- name
  dat_amr.upset <- dat_amr.upset %>%
    full_join(dat_amr.upset.sub, by = "MAG")
}
dat_amr.upset[is.na(dat_amr.upset)] <- 0

# upset plot helps show overlapping MAG drivers of difference AMR associations
upset(dat_amr.upset)

dat_amr.upset.anno <- dat_amr.upset %>%
  column_to_rownames("MAG") %>%
  mutate(rowsum = rowSums(.)) %>%
  rownames_to_column("MAG") %>%
  left_join(res.days_enrollment.anno, by = "MAG") %>%
  arrange(desc(logFC))

# What taxa are abundant at enrollment between the study sites ------------

# Generally, the question here is whether this is commonality between the
# microbial configurations of participants presenting with SAM. Obvious
# enteropathogens or dominant organisms?

# use the TPM-normalized dataset from kallisto
dat_tpm.filt.long <- dat_tpm.filt %>%
  rownames_to_column("SID") %>%
  left_join(map_fecal.proc %>% select(SID, area, study_phase), by = "SID") %>%
  filter(study_phase == "enrollment") %>%
  pivot_longer(colnames(dat_tpm.filt), names_to = "MAG", values_to = "abund_tpm") %>%
  left_join(tax %>% select(MAG, family, genus, genus_specie), by = "MAG")

dat_tpm.filt.summary <- dat_tpm.filt.long %>%
  group_by(area, MAG) %>%
  summarize(mean = mean(abund_tpm), 
            sd = sd(abund_tpm), 
            sem = sd(abund_tpm)/sqrt(length(abund_tpm))) %>%
  mutate(mean_frac = mean / 1E6)

dat_tpm.filt.summary.plot <- dat_tpm.filt.long %>%
  filter(MAG %in% (dat_tpm.filt.summary %>%
                     slice_max(mean, n = 10) %>%
                     pull(MAG))) %>%
  arrange(genus, desc(abund_tpm)) %>%
  mutate(MAG = factor(MAG, levels = unique(MAG)))

pdf(paste(datestring, study, "abundant_mags_at_enrollment.pdf", sep = "_"), height = 8.5, width = 11)
ggplot(dat_tpm.filt.summary.plot, aes(x = MAG, y = abund_tpm, fill = genus)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(area ~ .) +
  scale_fill_manual(values = cbbPalette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))
dev.off()

# Looking at whether aspects of anthropometric recovery are relate --------

# testing the relationship of rate of wlz response/time in acute treatment to enrollment MAG abundance

# 1 - does enrollment abundance of individual MAGs predict rate of response?

metadata.enrollment <- map_fecal.proc %>%
  filter(SID %in% rownames(dat_tpm.filt),
         study_phase %in% c("enrollment")) %>%
  left_join(dat_anth.anno %>% select(pid, study_phase, weight, length, waz, laz, wlz, muac), by = c("PID" = "pid", "study_phase")) %>%
  column_to_rownames("SID")

mags <- rownames(dat_vst.filt)

enrollment_config_summary <- data.frame()

dat_vst.filt.t.enrollment <- dat_vst.filt %>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG, names_to = "SID") %>% 
  pivot_wider(names_from=MAG, values_from=value) %>%
  inner_join(metadata.enrollment %>% 
               rownames_to_column("SID") %>% 
               select(-wlz, -area, -gender), 
             by = "SID") %>%
  left_join(acute_recovery_metrics %>% 
              filter(metric == "wlz_rate") %>%
              mutate(wlz_rate = value), 
            by = c("PID" = "pid"))

for (i in 1:length(mags)) {
  mag_abundance <- dat_vst.filt.t.enrollment[[mags[i]]]
#  mod <- lm(wlz ~ area + gender + mag_abundance, data = dat_vst.filt.t.enrollment)
  mod <- lm(wlz_rate ~ area + gender + mag_abundance, data = dat_vst.filt.t.enrollment)
  ano <- anova(mod)
  
  enrollment_config_summary <- enrollment_config_summary %>%
    bind_rows(data.frame(MAG = mags[i],
                         response = "wlz_rate",
                         beta = get_coef(mod, "mag_abundance"),
                         pval = get_anova_p(ano, "mag_abundance")))
}

enrollment_config_summary <- enrollment_config_summary %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(padj) %>%
  left_join(tax %>% select(MAG, genus, genus_specie), by = "MAG")

ggplot(enrollment_config_summary, aes(y = -1*log10(padj), x = beta)) +
  geom_point() +
  geom_hline(yintercept = -1*log10(0.1))

# 2 - does rate of response relate to rate of change of MAG abundance?

metadata.acute_rehab <- map_fecal.proc %>%
  filter(SID %in% colnames(dat_vst.filt),
         study_phase %in% c("enrollment", "acute_rehab", "baseline")) %>%
  left_join(dat_anth.anno %>% select(pid, study_phase, weight, length, waz, laz, wlz, muac), by = c("PID" = "pid", "study_phase")) %>%
  column_to_rownames("SID")

dat_vst.filt.t.acute_rehab <- dat_vst.filt %>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG, names_to = "SID") %>% 
  pivot_wider(names_from=MAG, values_from=value) %>%
  inner_join(metadata.acute_rehab %>% 
               rownames_to_column("SID") %>% 
               select(-wlz, -area, -gender), 
             by = "SID") %>%
  left_join(acute_recovery_metrics %>% 
              filter(metric == "wlz_rate") %>%
              mutate(wlz_rate = value), 
            by = c("PID" = "pid")) %>%
  mutate(weeks_enrollment = days_enrollment/7,
         gender = factor(as.character(gender)),
         area = factor(area),
         PID = factor(PID))


# WLZ-associated MAGs in the 'SAM phase' ----------------------------------

# I have rates of WLZ change for each participant. I want to relate those rates to rates of change of various MAGs. So, I need the rate of change
# of MAG abundance for each MAG in each participant.

# This is necessary since the fecal and anthropometry sampling timepoints don't
# align in the SAM phase

mags_vs_time_per_pid <- data.frame()

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(mags)*length(participants))
for (i in 1:length(mags)) {
  for (j in 1:length(participants)) {
    pb$tick()
    dat_vst.filt.t.acute_rehab.sub <- dat_vst.filt.t.acute_rehab %>%
      filter(PID %in% participants[j]) %>%
      select(mags[i], days_enrollment)
    mag_abundance <- dat_vst.filt.t.acute_rehab.sub[[mags[i]]]
    mod <- lm(mag_abundance ~ days_enrollment, data = dat_vst.filt.t.acute_rehab.sub)
    mags_vs_time_per_pid <- mags_vs_time_per_pid %>%
      bind_rows(data.frame(pid = participants[j],
                           mag = mags[i],
                           beta_days_enrollment = get_coef(mod, "days_enrollment")))
  }
}

acute_recovery_metrics.wide <- acute_recovery_metrics %>%
  pivot_wider(names_from = metric, values_from = value)

mags_vs_time_per_pid.wide <- mags_vs_time_per_pid %>%
  pivot_wider(id_cols = "pid", names_from = "mag", values_from = "beta_days_enrollment") %>%
  left_join(acute_recovery_metrics.wide, by = "pid")

rate_vs_rate_summary <- data.frame()

for (i in 1:length(mags)) {
  mag_betas <- mags_vs_time_per_pid.wide[[mags[i]]]
  mod <- lm(wlz_rate ~ mag_betas + area, data = mags_vs_time_per_pid.wide)
  ano <- anova(mod)
  rate_vs_rate_summary <- rate_vs_rate_summary %>%
    bind_rows(data.frame(MAG = mags[i],
                         beta_rate_rate = get_coef(mod, "mag_betas"),
                         pval_rate_rate = get_anova_p(ano, "mag_betas"),
                         beta_area = get_coef(mod, "areaKurigram"),
                         pval_area = get_anova_p(ano, "area")))
}

rate_vs_rate_summary <- rate_vs_rate_summary %>%
  mutate(padj_rate_rate = p.adjust(pval_rate_rate, method = "fdr"),
         padj_area = p.adjust(pval_area, method = "fdr")) %>%
  arrange(padj_rate_rate) %>%
  left_join(tax %>% select(MAG, genus, genus_specie), by = "MAG")

mags_vs_time_per_pid.summary1 <- mags_vs_time_per_pid %>%
  group_by(mag) %>%
  summarize(mean_beta_days_enrollment = mean(beta_days_enrollment),
            median_beta_days_enrollment = median(beta_days_enrollment),
            sd_beta_days_enrollment = sd(beta_days_enrollment),
            se_beta_days_enrollment = sd(beta_days_enrollment)/sqrt(length(beta_days_enrollment))) %>%
  arrange(mag)

mags_vs_time_per_pid.summary <- mags_vs_time_per_pid %>%
  left_join(acute_recovery_metrics.wide %>%
              select(pid, area) %>%
              distinct(),
            by = "pid") %>%
  group_by(area, mag) %>%
  summarize(mean_beta_days_enrollment = mean(beta_days_enrollment),
            median_beta_days_enrollment = median(beta_days_enrollment),
            sd_beta_days_enrollment = sd(beta_days_enrollment),
            se_beta_days_enrollment = sd(beta_days_enrollment)/sqrt(length(beta_days_enrollment))) %>%
  arrange(mag, area) %>%
  pivot_wider(id_cols = mag, names_from = area, values_from = c(mean_beta_days_enrollment, median_beta_days_enrollment, sd_beta_days_enrollment, se_beta_days_enrollment)) %>%
  left_join(mags_vs_time_per_pid.summary1, by = "mag")

acute_recovery_metrics.wide.summary <- acute_recovery_metrics.wide %>%
  summarize(area = "all",
            n = length(pid),
            mean_starting_wlz = mean(wlz),
            sd_starting_wlz = sd(wlz),
            mean_wlz_rate = mean(wlz_rate),
            sd_wlz_rate = sd(wlz_rate),
            mean_waz_rate = mean(waz_rate),
            sd_waz_rate = sd(waz_rate),
            mean_laz_rate = mean(laz_rate),
            sd_laz_rate = sd(laz_rate),
            mean_muac_rate = mean(muac_rate),
            sd_muac_rate = sd(muac_rate))

acute_recovery_metrics.wide.summary.area <- acute_recovery_metrics.wide %>%
  group_by(area) %>%
  summarize(n = length(pid),
            mean_starting_wlz = mean(wlz),
            sd_starting_wlz = sd(wlz),
            mean_wlz_rate = mean(wlz_rate),
            sd_wlz_rate = sd(wlz_rate),
            mean_waz_rate = mean(waz_rate),
            sd_waz_rate = sd(waz_rate),
            mean_laz_rate = mean(laz_rate),
            sd_laz_rate = sd(laz_rate),
            mean_muac_rate = mean(muac_rate),
            sd_muac_rate = sd(muac_rate))

acute_recovery_metrics.wide.summary <- acute_recovery_metrics.wide.summary %>%
  bind_rows(acute_recovery_metrics.wide.summary.area) %>%
  pivot_wider(names_from = area, values_from = -area)

rate_vs_rate_summary.anno <- rate_vs_rate_summary %>%
  left_join(mags_vs_time_per_pid.summary, by = c("MAG" = "mag")) %>%
  mutate(signed_log10_padj = -1*log10(padj_rate_rate)*sign(beta_rate_rate)) %>%
  arrange(signed_log10_padj)

write.table(rate_vs_rate_summary.anno, file = paste0(datestring, "_MDCF_POC_pSM_rate_versus_rate.txt"), sep = "\t", row.names = FALSE)

rate_vs_rate_summary.anno.sig <- rate_vs_rate_summary.anno %>%
  filter(padj_rate_rate < 0.05) %>%
  #  arrange(beta_rate_rate) %>%
  mutate(MAG = factor(MAG, levels = MAG)) %>%
  bind_cols(acute_recovery_metrics.wide.summary)

sum(rate_vs_rate_summary.anno.sig$beta_rate_rate < 0)

# plotting to visualize results

wlz_mags_acute <- rate_vs_rate_summary.anno.sig %>%
  arrange(beta_rate_rate) %>%
  pull(MAG)

pdf(paste(datestring, study, "abundance_rate_vs_wlz_rate.pdf", sep = "_"))
for (i in 1:length(wlz_mags_acute)) {
  sub <- rate_vs_rate_summary.anno.sig %>%
    filter(MAG == wlz_mags_acute[i])
  
  p <- ggplot(mags_vs_time_per_pid.wide, aes(x = .data[[as.character(wlz_mags_acute[i])]], y = wlz_rate)) +
    geom_point(aes(color = area)) +
    geom_smooth(method = "lm") +
    ggtitle(paste0("beta_rr = ", round(sub$beta_rate_rate, 3), " | padj_rr = ", round(sub$padj_rate_rate, 3)),
            paste0("beta_area = ", round(sub$beta_area, 3), " | padj_rr = ", round(sub$padj_area, 3))) +
    xlab(paste(wlz_mags_acute[i], sub$genus_specie, sep = "|")) +
    theme_classic()
  print(p)
  
}
dev.off()

# functions ---------------------------------------------------------------

get_anova_p <- function(anova, term) {
  anova_df <- as.data.frame(anova) %>%
    rownames_to_column(var = "term")
  pval <- anova_df[anova_df$term == term, "Pr(>F)"]
  return(pval)
}

get_confint <- function(model, term, level = 0.95) {
  if (class(model) == "lm") {
    confint_df <- as.data.frame(confint(model, level = level)) %>%
      rownames_to_column("term")
  } else if (class(model) == "lmerModLmerTest") {
    confint_df <- as.data.frame(confint(profile(model), level = level)) %>%
      rownames_to_column("term")
  } else {
    stop("ERROR: model class not recognized!")
  }
  lci <- confint_df[confint_df$term == term,]$`2.5 %`
  uci <- confint_df[confint_df$term == term,]$`97.5 %`
  if (length(lci) > 0) {
    return(c(lci, uci))
  } else {
    return(paste(NA, NA, sep = ";"))
  }
  
}

get_coef <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef)
  }
}
