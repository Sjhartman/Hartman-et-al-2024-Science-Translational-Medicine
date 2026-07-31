# load packages -----------------------------------------------------------

library(readxl)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(FSA)

# set up environment ------------------------------------------------------

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC_pSM"

setwd("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220523_enteropathogen_data/")

# load enteropathogen abundances from ddPCR -------------------------------

dat_entero_1 <- read.table("Gordon_TummyBug_05-19-2022/Final_Results_allSamples.xls", sep = "\t", header = TRUE) %>%
  mutate(batch = "A")
dat_entero_2 <- read.table("Gordon_TummyBug_05-20-2022/Final_Results_allSamples.xls", sep = "\t", header = TRUE) %>%
  mutate(batch = "B")

dat_entero <- dat_entero_1 %>%
  bind_rows(dat_entero_2) %>%
  mutate(SID = gsub("_P|_Q", "", Sample))

# load metadata -----------------------------------------------------------

map_fecal <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "fecal_map")

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

# join the enteropathogen abundances to the metadata
dat_entero.anno <- dat_entero %>%
  left_join(map_fecal.proc, by = "SID") %>%
  mutate(study_phase = factor(study_phase, levels = c("enrollment", "baseline", "intervention")))

# starting to look at abundances ------------------------------------------

# expecting enrollment, baseline, end of intervention samples
dat_entero.anno %>%
  filter(PID %in% "S01C888") %>%
  select(study_phase, days_enrollment) %>%
  distinct()

# some summary statistics
dat_entero.anno.stat <- dat_entero.anno %>%
  group_by(area, group, study_phase, Bug) %>%
  summarize(mean = mean(Pred.Conc),
            sd = sd(Pred.Conc),
            n = n()) %>%
  mutate(study_phase = factor(study_phase, levels = c("enrollment", "baseline", "intervention")))


# plotting to look at enteropathogen abundances ---------------------------

enteropathogens <- unique(dat_entero.anno.stat$Bug)[unique(dat_entero.anno.stat$Bug) != "16S"]

pdf("pSM_enteropathogens.pdf")
for (i in 1:length(enteropathogens)) {
  p <- ggplot(dat_entero.anno %>% filter(Bug == enteropathogens[i]), aes(x = study_phase, y = Pred.Conc)) +
    geom_boxplot() +
    ggtitle(enteropathogens[i]) +
    facet_grid(group ~ area) +
    theme_bw()
  print(p)
}
dev.off()

# enrollment comparisons --------------------------------------------------

# for each enteropathogen, run a kruskal wallis test and calculate summary
# statistics, then store the results
res_enrollment <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase == "enrollment")
  
  res <- kruskal.test(Pred.Conc ~ area, data = dat_entero.anno.sub)
  res_enrollment <- res_enrollment %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka = dat_entero.anno.sub %>% filter(area == "Dhaka") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka = dat_entero.anno.sub %>% filter(area == "Dhaka") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram = dat_entero.anno.sub %>% filter(area == "Kurigram") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram = dat_entero.anno.sub %>% filter(area == "Kurigram") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

# run p val adjustements
res_enrollment$padj <- p.adjust(res_enrollment$pval, method = "fdr")
res_enrollment %>%
  arrange(padj)

write.table(res_enrollment, file = paste0(datestring, "_", study, "_enrollment_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# baseline comparisons ----------------------------------------------------

# as above, at baseline
res_baseline <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase == "baseline")
  
  res <- kruskal.test(Pred.Conc ~ area, data = dat_entero.anno.sub)
  res_baseline <- res_baseline %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka = dat_entero.anno.sub %>% filter(area == "Dhaka") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka = dat_entero.anno.sub %>% filter(area == "Dhaka") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram = dat_entero.anno.sub %>% filter(area == "Kurigram") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram = dat_entero.anno.sub %>% filter(area == "Kurigram") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

res_baseline <- res_baseline %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         ratio = mean_kurigram/mean_dhaka) %>%
  arrange(padj)

write.table(res_baseline, file = paste0(datestring, "_", study, "_baseline_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# baseline comparisons by group -------------------------------------------

# as above, using a dummy variable combining area and group at baseline
res_baseline <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase == "baseline") %>%
    mutate(area_group = paste(area, group, sep = "_"))
  
  res <- kruskal.test(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  dunnTest(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res_baseline <- res_baseline %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka_MDCF2 = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "MDCF-2") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_MDCF2 = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "MDCF-2") %>% pull(Pred.Conc) %>% sd,
                         mean_dhaka_RUSF = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "RUSF") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_RUSF = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "RUSF") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_MDCF2 = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "MDCF-2") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_MDCF2 = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "MDCF-2") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_RUSF = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "RUSF") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_RUSF = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "RUSF") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

res_baseline <- res_baseline %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(padj)

write.table(res_baseline, file = paste0(datestring, "_", study, "_baseline_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# treatment comparisons by group -------------------------------------------

# treatment phase
res_trt <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase == "intervention") %>%
    mutate(area_group = paste(area, group, sep = "_"))
  
  res <- kruskal.test(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
#  dunnTest(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res_trt <- res_trt %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka_MDCF2 = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "MDCF-2") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_MDCF2 = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "MDCF-2") %>% pull(Pred.Conc) %>% sd,
                         mean_dhaka_RUSF = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "RUSF") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_RUSF = dat_entero.anno.sub %>% filter(area == "Dhaka" & group == "RUSF") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_MDCF2 = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "MDCF-2") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_MDCF2 = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "MDCF-2") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_RUSF = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "RUSF") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_RUSF = dat_entero.anno.sub %>% filter(area == "Kurigram" & group == "RUSF") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

res_trt <- res_trt %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(padj)

write.table(res_trt, file = paste0(datestring, "_", study, "_treatment_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# SAM vs post-SAM MAM -------------------------------------------

# comparing enrollment to baseline, using paired nonparametric tests
res_sam_vs_psm <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase %in% c("enrollment", "baseline")) %>%
    mutate(area_group = paste(area, group, sep = "_")) %>%
    add_count(PID) %>%
    filter(n == 2) %>%
    arrange(PID)
  
#  res <- kruskal.test(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res <- wilcox.test(Pred.Conc ~ study_phase, data = dat_entero.anno.sub, paired = TRUE)
  #  dunnTest(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res_sam_vs_psm <- res_sam_vs_psm %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka_sam = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "enrollment") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_sam = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "enrollment") %>% pull(Pred.Conc) %>% sd,
                         mean_dhaka_psm = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "baseline") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_psm = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "baseline") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_sam = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "enrollment") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_sam = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "enrollment") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_psm = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "baseline") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_psm = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "baseline") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

res_sam_vs_psm <- res_sam_vs_psm %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(padj)

write.table(res_sam_vs_psm, file = paste0(datestring, "_", study, "_sam_vs_psm_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# post-SAM MAM vs treatment -------------------------------------------

# comparing baseline to treatment, using paired nonparametric tests
res_psm_vs_int <- data.frame()
for (enteropathogen in enteropathogens) {
  dat_entero.anno.sub <- dat_entero.anno %>%
    filter(Bug == enteropathogen & study_phase %in% c("baseline", "intervention")) %>%
    mutate(area_group = paste(area, group, sep = "_")) %>%
    add_count(PID) %>%
    filter(n == 2) %>%
    arrange(PID)
  
  #  res <- kruskal.test(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res <- wilcox.test(Pred.Conc ~ study_phase, data = dat_entero.anno.sub, paired = TRUE)
  #  dunnTest(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
  res_psm_vs_int <- res_psm_vs_int %>%
    bind_rows(data.frame(pathogen = enteropathogen,
                         mean_dhaka_psm = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "baseline") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_psm = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "baseline") %>% pull(Pred.Conc) %>% sd,
                         mean_dhaka_int = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "intervention") %>% pull(Pred.Conc) %>% mean,
                         sd_dhaka_int = dat_entero.anno.sub %>% filter(area == "Dhaka" & study_phase == "intervention") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_psm = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "baseline") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_psm = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "baseline") %>% pull(Pred.Conc) %>% sd,
                         mean_kurigram_int = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "intervention") %>% pull(Pred.Conc) %>% mean,
                         sd_kurigram_int = dat_entero.anno.sub %>% filter(area == "Kurigram" & study_phase == "intervention") %>% pull(Pred.Conc) %>% sd,
                         pval = res$p.value))
}

res_psm_vs_int <- res_psm_vs_int %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(padj)

write.table(res_psm_vs_int, file = paste0(datestring, "_", study, "_psm_vs_int_enteropathogens.txt"), row.names = FALSE, sep = "\t")

# post-SAM MAM vs treatment, by group  -------------------------------------------

# baseline versus treatment, using paired test and dummy variable combining group and phase
res_psm_vs_int2 <- data.frame()
for (a in c("Dhaka", "Kurigram")) {
  for (g in c("MDCF-2", "RUSF")) {
    for (enteropathogen in enteropathogens) {
      dat_entero.anno.sub <- dat_entero.anno %>%
        filter(Bug == enteropathogen & study_phase %in% c("baseline", "intervention") & area == a & group == g) %>%
        #    mutate(group_phase = paste(group, study_phase, sep = "_")) %>%
        add_count(PID) %>%
        filter(n == 2) %>%
        arrange(PID)
      
      #  res <- kruskal.test(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
      res <- wilcox.test(Pred.Conc ~ study_phase, data = dat_entero.anno.sub, paired = TRUE)
      #  dunnTest(Pred.Conc ~ area_group, data = dat_entero.anno.sub)
      res_psm_vs_int2 <- res_psm_vs_int2 %>%
        bind_rows(data.frame(area = a,
                             group = g,
                             pathogen = enteropathogen,
                             mean_psm_MDCF_2 = dat_entero.anno.sub %>% filter(study_phase == "baseline") %>% pull(Pred.Conc) %>% mean,
                             sd_psm_MDCF_2 = dat_entero.anno.sub %>% filter(study_phase == "baseline") %>% pull(Pred.Conc) %>% sd,
                             mean_int_MDCF_2 = dat_entero.anno.sub %>% filter(study_phase == "intervention") %>% pull(Pred.Conc) %>% mean,
                             sd_int_MDCF_2 = dat_entero.anno.sub %>% filter(study_phase == "intervention") %>% pull(Pred.Conc) %>% sd,
                             pval = res$p.value))
    }}}

res_psm_vs_int2 <- res_psm_vs_int2 %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(area, group, padj)

write.table(res_psm_vs_int2, file = paste0(datestring, "_", study, "_psm_vs_int_enteropathogens_by_group.txt"), row.names = FALSE, sep = "\t")

# relating enteropathogen burden to WLZ -----------------------------------

# load the anthropometry data
dat_anth <- read_excel("~/Documents/Projects/human_studies/MDCF POC Post-SAM MAM/220112 abundance re-analysis/220112_pSM_mapping_file.xlsx", sheet = "anthro_dat")
anth_key <- data.frame(time = c(1:9),
                       study_age = c(-1, 0, 15, 30, 45, 60, 75, 90, 120),
                       study_week = c(-1, 0, 2, 4, 6, 8, 10, 12, 16),
                       study_phase = c("enrollment", "baseline", rep("intervention", 6), "follow_up_1mo"))
dat_anth.anno <- dat_anth %>%
  left_join(anth_key, by = "time")

# annotate the enteropathogen data with anthropometry data
dat_entero.anno.anth <- dat_entero %>%
  left_join(map_fecal.proc, by = "SID") %>%
  mutate(study_week = ifelse(is.na(study_week), -1, study_week)) %>%
  inner_join(dat_anth.anno %>% 
               mutate(age_days = agedays) %>%
               select(pid, study_phase, study_week, weight, length, waz, laz, wlz, muac, age_days), by = c("PID" = "pid", "study_phase", "study_week")) %>%
  mutate(group = factor(group, levels = c("RUSF", "MDCF-2")))

# are there enteropathogen abundances associated with WLZ at enrollment (cross-section)
# using linear models, controlling for area (study site)
dat_entero.anno.anth.sam <- dat_entero.anno.anth %>%
  filter(study_phase == "enrollment")

res <- data.frame()

for (i in 1:length(enteropathogens)) {
  
  dat_entero.anno.anth.sam.sub <- dat_entero.anno.anth.sam %>%
    filter(Bug == enteropathogens[i])
  
  mod <- lm(wlz ~ Pred.Conc + area, data = dat_entero.anno.anth.sam.sub)
  
  res <- res %>%
    bind_rows(get_coef_p(mod, "Pred.Conc") %>%
                mutate(enteropathogen = enteropathogens[i],
                       metric = "wlz"))
  
}

res <- res %>%
  mutate(padj = p.adjust(pval, "fdr")) %>%
  arrange(padj)

res %>%
  filter(padj < 0.05)

# does the number of enteropathogens confidently detected relate to WLZ?
dat_entero.anno.anth.sam.bin <- dat_entero.anno.anth.sam %>%
  group_by(SID, PID, area, group, study_phase, wlz) %>%
  summarize(n_enteropathogens = sum(Result_Status == "Positive")) %>%
  mutate(group = factor(group, levels = c("RUSF", "MDCF-2")))

mod <- lm(wlz ~ n_enteropathogens + area, data = dat_entero.anno.anth.sam.bin)
summary(mod)

# functions ---------------------------------------------------------------

get_coef_p <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  pval <- coef_df[term, "Pr(>|t|)"]
  return(data.frame(term = term, coef = coef, pval = pval))
}

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

get_std_coef <- function(model, term, response, predictor) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  coef_std <- coef*(sd(predictor)/sd(response))
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef_std)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef_std)
  }
}

get_intercept <- function(model) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df["(Intercept)", "Estimate"]
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef)) {
    stop(paste0("ERROR: term ", term, " does not exist in model!"))
  } else {
    return(coef)
  }
}

sem <- function(x, na.rm) sd(x, na.rm = na.rm)/sqrt(length(x[!is.na(x)]))

get_model_se <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Std. Error"]
  if (is.na(coef)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef)
  }
}
