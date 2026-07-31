# loading libraries -------------------------------------------------------

library(tidyverse)
library(scales)
library(lmerTest)
library(readxl)
library(ggbeeswarm)

# setting up the environment ----------------------------------------------

# setwd("~/Library/CloudStorage/Box-Box/Projects/human_studies/MDCF POC Post-SAM MAM/221006 Anthropometry for tables/")

# for file naming
datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC_pSM"

# custom color palette
cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

# loading data ------------------------------------------------------------

# anthropometry
dat_anthropometry <- read_excel("/Users/stevenhartman/Downloads/AnthropometryDataPostSAMMAM.xlsx", sheet = "AnthropometryDataPostSAMMAM")
anthropometry_week_key <- read_excel("/Users/stevenhartman/Downloads/AnthropometryDataPostSAMMAM.xlsx", sheet = "key")
dat_anthropometry <- dat_anthropometry %>%
  left_join(anthropometry_week_key, by = c("time" = "timepoint")) %>%
  mutate(group = factor(group, levels = c("RUSF", "MDCF-2")))

# how many participants in each group/site, regardless of trial compliance?
length(unique(dat_anthropometry$pid))
dat_anthropometry %>%
  group_by(area, group) %>%
  summarize(length(unique(pid)))

# remove participants who did not complete the trial per-protocol
dat_anthropometry.no_dropout <- dat_anthropometry %>%
  filter(dropout_status == "non_dropout")

# how many participants after dropout removal?
dat_anthropometry.no_dropout %>%
  group_by(area, group) %>%
  summarize(length(unique(pid)))

dat_anthropometry.no_dropout.trt <- dat_anthropometry.no_dropout %>%
  filter(!study_phase %in% c("enrollment", "followup"))

dat_acute_rehab <- dat_anthropometry.no_dropout %>%
  select(pid, agedays, study_phase) %>%
  filter(study_phase %in% c("enrollment", "baseline")) %>%
  pivot_wider(id_cols = pid, names_from = study_phase, values_from = agedays) %>%
  mutate(acute_rehab_time = baseline - enrollment)

# first look at WLZ in participants prior to flood analysis ---------------

# lme
mod <- lmer(wlz ~ gender + area + group + study_week + group*study_week + (1|pid), data = dat_anthropometry.no_dropout.trt)
anova <- stats::anova(mod, ddf="Kenward-Roger") 

# interpret the coefficients of interest
summary(mod)
get_coef(mod, "groupMDCF-2:study_week")
get_anova_p(anova, "group:study_week") # not significant

# plot the effect
ggplot(dat_anthropometry.no_dropout.trt, aes(x = study_week, y = wlz, color = group)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
#  scale_color_manual(values = c("red", "black")) +
  theme_classic()

# flood timing analysis ---------------------------------------------------
# We were informed of the flood in Kurigram affecting a large proportion of our trial
# participants. Here we explore the effect of this flood on participant outcomes.

# load a dataset describing which participant were affected by the flood, and begin
# to organize the data to understand precisely when the flood occurred
dat_flood <- read_excel("/Users/stevenhartman/Library/CloudStorage/Box-Box/000Gordon/Projects/PostSamMam/Data/FLOOD (01-Mar-21).xlsx") %>%
  left_join(dat_acute_rehab %>% select(pid, acute_rehab_time), by = "pid") %>%
  filter(area == "Kurigram" & pid %in% dat_anthropometry.no_dropout$pid) %>%
  arrange(Feeding_start_date) %>%
  mutate(feeding_start_day_of_year = as.numeric(strftime(Feeding_start_date, format = "%j")),
         feeding_start_month = as.numeric(strftime(Feeding_start_date, format = "%m")),
         feeding_start_season = cut(feeding_start_month, breaks = c(1, 3, 6, 11, 12), labels = c("cool_and_dry", "hot", "monsoon", "cool_and_dry"), right = FALSE, include.lowest = TRUE),
         feeding_start_season = factor(feeding_start_season, levels = c("cool_and_dry", "hot", "monsoon")),
         flood_time = ifelse(pid %in% c("S50C999", "S49C999", "S47C999", "S48C999"), "affected_acute_rehab", flood_time),
         flood_time = ifelse(pid %in% c("S03C999", "S02C999", "S05C999", "S04C999", "S01C999"), "affected_followup", flood_time),
         flood_time = factor(flood_time, levels = c("affected_acute_rehab", 1, 2, 3, 4, "affected_followup", 0), 
                             labels = c("affected_acute_rehab", "affected_first_quartile", "affected_second_quartile", "affected_third_quartile", "affected_fourth_quartile", "affected_followup", "not_affected")),
         pid = factor(pid, levels = pid),
         group = ifelse(group == "Arm-1", "MDCF-2", "RUSF"))

# how many children were affected by the flood?
sum(dat_flood$flood_time == "not_affected")
sum(dat_flood$flood_time != "not_affected")

# March-May: pre-monsoon hot season
# June-October: rainy
# November-February: cool and dry

# plot the flood timing
# pdf(paste0(datestring, "_MDCF_POC_pSM_flood_timing_in_Kurigram.pdf"), height = 11, width = 8.5)
ggplot(dat_flood, aes(x = Feeding_start_date, y = pid)) +
  geom_errorbarh(aes(xmin = Feeding_start_date - acute_rehab_time*24*60*60, xmax = Feeding_end_date + 30*24*60*60), color = "black") +
  geom_segment(aes(x = Feeding_start_date, xend = Feeding_end_date, y = pid, yend = pid, color = flood_time), size = 2) +
  geom_segment(aes(x = flood_start_date, xend = flood_end_date, y = pid, yend = pid), size = 2, color = "red") +
  scale_color_manual(values = cbbPalette) +
  theme_classic() +
  facet_wrap(~ group, scales = "free_y")
# dev.off()

# version 1 - January list (most conservative) - excludes participants who experienced flooding in acute rehab or follow-up in addition to treatment phase
dat_anthropometry.no_dropout.no_flood <- dat_anthropometry %>%
  filter(dropout_status == "non_dropout" & flood_status_jan == "unaffected") %>%
  left_join(dat_flood %>% select(pid, Feeding_start_date, feeding_start_day_of_year, feeding_start_month, feeding_start_season), by = "pid")
length(unique(dat_anthropometry.no_dropout.no_flood$pid))

# check study-wide numbers after excluding flood-affected participants
dat_anthropometry.no_dropout.no_flood %>%
  group_by(area, group) %>%
  summarize(length(unique(pid)))

# anthropometric effects of flooding in Kurigram --------------------------

anova(lm(wlz ~ area + flood_status_jan, dat_anthropometry.no_dropout %>% filter(study_phase == "enrollment")))

anova(lm(wlz ~ flood_status_jan, dat_anthropometry.no_dropout %>% filter(study_phase == "enrollment" & area == "Kurigram")))

# did the anthropometry of children affected by flooding differ at enrollment?
x <- dat_anthropometry.no_dropout %>% 
  filter(study_phase == "enrollment" & area == "Kurigram" & flood_status_jan == "unaffected")
y <- dat_anthropometry.no_dropout %>% 
  filter(study_phase == "enrollment" & area == "Kurigram" & flood_status_jan == "flood_affected")

t.test(x$wlz, y$wlz)
t.test(x$laz, y$laz)
t.test(x$muac, y$muac)
t.test(x$waz, y$waz)

anova(lm(wlz ~ flood_status_jan, dat_anthropometry.no_dropout %>% filter(study_phase == "baseline" & area == "Kurigram")))

# did the anthropometry of children affected by flooding differ at baseline?
x <- dat_anthropometry.no_dropout %>% 
  filter(study_phase == "baseline" & area == "Kurigram" & flood_status_jan == "unaffected")
y <- dat_anthropometry.no_dropout %>% 
  filter(study_phase == "baseline" & area == "Kurigram" & flood_status_jan == "flood_affected")

t.test(x$wlz, y$wlz)
t.test(x$laz, y$laz)
t.test(x$muac, y$muac)
t.test(x$waz, y$waz)

# visualize these effects
ggplot(dat_anthropometry.no_dropout %>% filter(study_phase == "enrollment"), aes(x = paste(area, flood_status_jan), y = wlz)) +
  geom_beeswarm()

ggplot(dat_anthropometry.no_dropout %>% filter(study_phase == "enrollment"), aes(x = paste(area, flood_status_jan), y = wlz)) +
  geom_boxplot()

ggplot(dat_anthropometry.no_dropout %>% filter(study_phase == "baseline"), aes(x = paste(area, flood_status_jan), y = wlz)) +
  geom_boxplot()

dat_anthropometry.no_dropout.kurigram <- dat_anthropometry.no_dropout %>%
  filter(area == "Kurigram" & !study_phase %in% c("enrollment"))

# looking specifically at the treatment phase
dat_anthropometry.no_dropout.kurigram.trt <- dat_anthropometry.no_dropout %>%
  filter(area == "Kurigram" & study_phase %in% c("baseline", "treatment"))

ggplot(dat_anthropometry.no_dropout.kurigram.trt, aes(x = study_week, y = wlz, color = flood_status_jan)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  scale_color_manual(values = c("red", "black")) +
  theme_classic()

# statistical analysis
dat <- dat_anthropometry.no_dropout.kurigram
treatment_groups <- c("MDCF-2", "RUSF", "all")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("flood_status_janunaffected", "flood_status_janunaffected:study_week")
anova_terms <- c("flood_status_jan", "flood_status_jan:study_week")

stat_summary.flood <- data.frame()
dat_tx <- filter(dat, study_week %in% 0:12)

for (treatment_group in treatment_groups) {
  if (treatment_group == "all") {
    dat_sub <- dat_tx
  } else {
  dat_sub <- dat_tx %>%
    filter(group == treatment_group)
  }
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    if (treatment_group == "all") {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + flood_status_jan + flood_status_jan:study_week + group + study_week + group:study_week + (1|pid), data = dat_sub)
    } else {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + flood_status_jan + flood_status_jan:study_week + study_week + (1|pid), data = dat_sub)
    }
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary.flood <- stat_summary.flood %>%
      bind_rows(data.frame(data = "treatment_phase_only",
                           group = treatment_group,
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}
}

stat_summary.flood

# write results to file
# write.table(stat_summary.flood, file = paste0(datestring, "_", study, "_anthropometry_analysis_summary_kurigram_effect_of_flood.txt"), row.names = FALSE, sep = "\t")

# enrollment anthropometry (all participants) --------------------------------------------------

# filter dataset to just enrollment timepoint
dat_anthropometry.no_dropout.enrollment <- dat_anthropometry %>%
  filter(study_phase %in% c("enrollment") & dropout_status == "non_dropout")

# descriptive summaries
dat_anthropometry.no_dropout.enrollment %>%
  group_by(area) %>%
  summarize(n = length(pid),
          age_mean = mean(agedays),
          age_sd = sd(agedays),
          female_n = sum(gender == 2),
          male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.enrollment %>%
  group_by(area) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -area, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(area, metric), names_from = statistic, values_from = value) %>%
  transmute(area = area, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = area, values_from = statistic)

# write to file
# write.table(summary, "PSM_enrollment_summary_flood_included.txt", sep = "\t")

# run tests on metrics of interest
fisher.test(dat_anthropometry.no_dropout.enrollment %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.enrollment %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.316
t.test(wlz ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.0106
t.test(waz ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.0031
t.test(laz ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.0113
t.test(weight ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.008066
t.test(length ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.01369
t.test(muac ~ area, data = dat_anthropometry.no_dropout.enrollment) # P = 0.03529

# enrollment anthropometry (all participants, Kurigram) --------------------------------------------------

# enrollment analysis per site
dat_anthropometry.no_dropout.enrollment.kurigram <- dat_anthropometry %>%
  filter(study_phase %in% c("enrollment") & dropout_status == "non_dropout" & area == "Kurigram")

dat_anthropometry.no_dropout.enrollment.kurigram %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.enrollment.kurigram %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_enrollment_summary_flood_included_kurigram.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.enrollment %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.enrollment %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.3274
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.4962
t.test(waz ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.8152
t.test(laz ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.7001
t.test(weight ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.3811
t.test(length ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.2842
t.test(muac ~ group, data = dat_anthropometry.no_dropout.enrollment.kurigram) # P = 0.4113

# enrollment anthropometry (all participants, Dhaka) --------------------------------------------------

dat_anthropometry.no_dropout.enrollment.dhaka <- dat_anthropometry %>%
  filter(study_phase %in% c("enrollment") & dropout_status == "non_dropout" & area == "Dhaka")

dat_anthropometry.no_dropout.enrollment.dhaka %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.enrollment.dhaka %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_enrollment_summary_flood_included_dhaka.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.enrollment %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.enrollment %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.584
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.1228
t.test(waz ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.7025
t.test(laz ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.875
t.test(weight ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.4967
t.test(length ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.4865
t.test(muac ~ group, data = dat_anthropometry.no_dropout.enrollment.dhaka) # P = 0.1308

# enrollment anthropometry (flood unaffected participants) --------------------------------------------------

# similar analysis as above, but looking at only flood-unaffected participants
dat_anthropometry.no_dropout.no_flood.enrollment <- dat_anthropometry %>%
  filter(study_phase %in% c("enrollment") & dropout_status == "non_dropout" & flood_status_jan == "unaffected")
length(unique(dat_anthropometry.no_dropout.no_flood.enrollment$pid))

dat_anthropometry.no_dropout.no_flood.enrollment %>%
  group_by(area) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.no_flood.enrollment %>%
  group_by(area) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -area, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(area, metric), names_from = statistic, values_from = value) %>%
  transmute(area = area, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = area, values_from = statistic)

# write.table(summary, "PSM_enrollment_summary_flood_excluded.txt", sep = "\t")

fisher.test(dat_anthropometry.no_dropout.no_flood.enrollment %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.no_flood.enrollment %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.603
t.test(wlz ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.002754
t.test(waz ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.001903
t.test(laz ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.007565
t.test(weight ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.0223
t.test(length ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.03466
t.test(muac ~ area, data = dat_anthropometry.no_dropout.no_flood.enrollment) # P = 0.01954

# enrollment anthropometry (flood unaffected participants, Kurigram) --------------------------------------------------

# flood-unaffected participants, restricted to the Kurigram site
dat_anthropometry.no_dropout.no_flood.enrollment.kurigram <- dat_anthropometry %>%
  filter(study_phase %in% c("enrollment") & dropout_status == "non_dropout" & flood_status_jan == "unaffected" & area == "Kurigram")
length(unique(dat_anthropometry.no_dropout.no_flood.enrollment.kurigram$pid))

dat_anthropometry.no_dropout.no_flood.enrollment.kurigram %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.no_flood.enrollment.kurigram %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_enrollment_summary_flood_excluded_kurigram.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.no_flood.enrollment %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.no_flood.enrollment %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.7843
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.2853
t.test(waz ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.6854
t.test(laz ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.9221
t.test(weight ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.9937
t.test(length ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.8883
t.test(muac ~ group, data = dat_anthropometry.no_dropout.no_flood.enrollment.kurigram) # P = 0.6468

# baseline anthropometry (all participants) --------------------------------------------------

# all participants, including flood-affected, that completed the trial
dat_anthropometry.no_dropout.baseline <- dat_anthropometry %>%
  filter(study_phase %in% c("baseline") & dropout_status == "non_dropout")
length(unique(dat_anthropometry.no_dropout.baseline$pid))

dat_anthropometry.no_dropout.baseline %>%
  group_by(area) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.baseline %>%
  group_by(area) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -area, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(area, metric), names_from = statistic, values_from = value) %>%
  transmute(area = area, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = area, values_from = statistic)

# write.table(summary, "PSM_baseline_summary_flood_included.txt", sep = "\t")

fisher.test(dat_anthropometry.no_dropout.baseline %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.baseline %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.3833
t.test(wlz ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.001069
t.test(waz ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.0003833
t.test(laz ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.01509
t.test(weight ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.002444
t.test(length ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.02191
t.test(muac ~ area, data = dat_anthropometry.no_dropout.baseline) # P = 0.001897

# baseline anthropometry (all participants, Kurigram) --------------------------------------------------

# similar analysis blocks as above, applied to the baseline (after SAM acute rehab)
dat_anthropometry.no_dropout.baseline.kurigram <- dat_anthropometry %>%
  filter(study_phase %in% c("baseline") & dropout_status == "non_dropout" & area == "Kurigram")
length(unique(dat_anthropometry.no_dropout.baseline.kurigram$pid))

dat_anthropometry.no_dropout.baseline.kurigram %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.baseline.kurigram %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_baseline_summary_flood_included_kurigram.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.baseline %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.baseline %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.3226
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.7432
t.test(waz ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.5826
t.test(laz ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.703
t.test(weight ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.2718
t.test(length ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.2856
t.test(muac ~ group, data = dat_anthropometry.no_dropout.baseline.kurigram) # P = 0.1573

# baseline anthropometry (all participants, Dhaka) --------------------------------------------------

dat_anthropometry.no_dropout.baseline.dhaka <- dat_anthropometry %>%
  filter(study_phase %in% c("baseline") & dropout_status == "non_dropout" & area == "Dhaka")

dat_anthropometry.no_dropout.baseline.dhaka %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.baseline.dhaka %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_baseline_summary_flood_included_dhaka.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.baseline %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.baseline %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.5265
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.2202
t.test(waz ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.5569
t.test(laz ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.8276
t.test(weight ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.5189
t.test(length ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.4828
t.test(muac ~ group, data = dat_anthropometry.no_dropout.baseline.dhaka) # P = 0.2033

# baseline anthropometry (flood-unaffected participants, Kurigram) --------------------------------------------------

dat_anthropometry.no_dropout.no_flood.baseline.kurigram <- dat_anthropometry %>%
  filter(study_phase %in% c("baseline") & dropout_status == "non_dropout" & area == "Kurigram" & flood_status_jan == "unaffected")

dat_anthropometry.no_dropout.no_flood.baseline.kurigram %>%
  group_by(group) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.no_flood.baseline.kurigram %>%
  group_by(group) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(group, metric), names_from = statistic, values_from = value) %>%
  transmute(group = group, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = group, values_from = statistic)

# write.table(summary, "PSM_baseline_summary_flood_excluded_kurigram.txt", sep = "\t")

#fisher.test(dat_anthropometry.no_dropout.baseline %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.baseline %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.792
t.test(wlz ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.6727
t.test(waz ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.7729
t.test(laz ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.9271
t.test(weight ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.9635
t.test(length ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.8883
t.test(muac ~ group, data = dat_anthropometry.no_dropout.no_flood.baseline.kurigram) # P = 0.6166

# baseline anthropometry (flood-unaffected participants) --------------------------------------------------

dat_anthropometry.no_dropout.no_flood.baseline <- dat_anthropometry %>%
  filter(study_phase %in% c("baseline") & dropout_status == "non_dropout" & flood_status_jan == "unaffected")
length(unique(dat_anthropometry.no_dropout.no_flood.baseline$pid))

dat_anthropometry.no_dropout.no_flood.baseline %>%
  group_by(area) %>%
  summarize(n = length(pid),
            age_mean = mean(agedays),
            age_sd = sd(agedays),
            female_n = sum(gender == 2),
            male_n = sum(gender == 1))

summary <- dat_anthropometry.no_dropout.no_flood.baseline %>%
  group_by(area) %>%
  summarize(wlz_mean = mean(wlz),
            wlz_sd = sd(wlz),
            waz_mean = mean(waz),
            waz_sd = sd(waz),
            laz_mean = mean(laz),
            laz_sd = sd(laz),
            weight_mean = mean(weight),
            weight_sd = sd(weight),
            length_mean = mean(length),
            length_sd = sd(length),
            muac_mean = mean(muac),
            muac_sd = sd(muac)) %>%
  pivot_longer(cols = -area, names_to = "statistic", values_to = "value") %>%
  separate(statistic, into = c("metric", "statistic"), sep = "_") %>%
  pivot_wider(id_cols = c(area, metric), names_from = statistic, values_from = value) %>%
  transmute(area = area, metric = metric, statistic = paste(round(mean, 2), round(sd, 2), sep = " ± ")) %>%
  pivot_wider(id_cols = metric, names_from = area, values_from = statistic)

# write.table(summary, "PSM_baseline_summary_flood_excluded.txt", sep = "\t")

fisher.test(dat_anthropometry.no_dropout.no_flood.baseline %>% filter(area == "Dhaka") %>% pull(gender), dat_anthropometry.no_dropout.no_flood.baseline %>% filter(area == "Kurigram") %>% pull(gender))
t.test(agedays ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.6893
t.test(wlz ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.0003257
t.test(waz ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.0002239
t.test(laz ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.01091
t.test(weight ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.008208
t.test(length ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.05551
t.test(muac ~ area, data = dat_anthropometry.no_dropout.no_flood.baseline) # P = 0.007891

# trial anthropometry responses -------------------------------------------
# these are the primary outcome analysis/analyses for the study

# Following our analysis of flooding effects on anthropometry, we subset the data to
# participants unaffected by flooding and who also adhered to the trial protocol
# Our first analyses concern the MAM treatment phase only
dat_anthropometry.no_dropout.no_flood.trt <- dat_anthropometry.no_dropout.no_flood %>%
  filter(study_phase %in% c("baseline", "treatment")) %>%
  mutate(treatment_phase = ifelse(study_week <= 2, "early", "late"))

# summarize the numbers of participants per treatment group/arm
dat_anthropometry.no_dropout.no_flood.trt %>%
  select(area, group, pid) %>%
  distinct() %>%
  group_by(area, group) %>%
  tally()

# generalizable looped modeling framework is used, so set the dataset, variables, etc. here
dat <- dat_anthropometry.no_dropout.no_flood.trt
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

# run the models and collect the results
stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
#  mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + treatment_phase + area:study_week + group*study_week + (1|pid), data = dat)
  anova <- stats::anova(mod, ddf="Kenward-Roger")
  for (i in 1:length(model_terms)) {
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_phase_only",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

stat_summary %>%
  arrange(term)

# ranova - test for importance of random effects
# ls_means - calculate and compare least squares means
# step - stepwise simplification

# trial plus follow-up ----------------------------------------------------

# Our next set of analyses concern the MAM treatment phase plus one month of follow-up
dat_anthropometry.no_dropout.with_followup <- dat_anthropometry.no_dropout %>%
  filter(!study_phase %in% c("enrollment"))

dat <- dat_anthropometry.no_dropout.with_followup
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

#stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_plus_followup",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

# trial anthropometry responses (all participants) -------------------------------------------

# Similar strategy applied to participants who completed the trial regardless to flood-affected status
dat_anthropometry.no_dropout.trt %>%
  select(area, group, pid) %>%
  distinct() %>%
  group_by(area, group) %>%
  tally()
unique(dat_anthropometry.no_dropout.trt$study_phase)

dat <- dat_anthropometry.no_dropout.trt
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_phase_only",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

stat_summary %>%
  arrange(term)

# trial plus follow-up (all participants) ----------------------------------------------------

# as above, including follow-up
dat_anthropometry.no_dropout.trt_pl_followup <- dat_anthropometry.no_dropout %>%
  filter(!study_phase %in% c("enrollment"))

dat <- dat_anthropometry.no_dropout.trt_pl_followup
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

#stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_plus_followup",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

# write.table(stat_summary, file = paste0(datestring, "_", study, "_anthropometry_analysis_summary_flood_included.txt"), row.names = FALSE, sep = "\t")

# trial plus follow-up (all participants, ITT) ----------------------------------------------------

# as above - this is a 'true' intent-to-treat for the MAM treatment phase and follow up
dat_anthropometry.trt_pl_followup <- dat_anthropometry %>%
  filter(!study_phase %in% c("enrollment"))

dat <- dat_anthropometry.trt_pl_followup
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_plus_followup",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

stat_summary %>%
  arrange(term)


# trial anthropometry responses (ITT) -------------------------------------------

# as above - this is a 'true' intent-to-treat for the MAM treatment phase only
dat_anthropometry.trt <- dat_anthropometry %>%
  filter(!study_phase %in% c("enrollment", "followup"))

dat_anthropometry.trt %>%
  select(area, group, pid) %>%
  distinct() %>%
  group_by(area, group) %>%
  tally()
unique(dat_anthropometry.no_dropout.trt$study_phase)

dat <- dat_anthropometry.trt
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week", "areaKurigram:study_week")
anova_terms <- c("group:study_week", "area:study_week")

stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + area + group + study_week + area:study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary <- stat_summary %>%
      bind_rows(data.frame(data = "treatment_phase_only",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

stat_summary %>%
  arrange(term)

# ranova - test for importance of random effects
# ls_means - calculate and compare least squares means
# step - stepwise simplification

# trial anthropometry responses by treatment group -------------------------------------------

# Above analyses included results form children across both arms of the trial. Here, we break down
# the analyses to look at results within each treatment group
dat_anthropometry.no_dropout.no_flood.trt <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment", "followup"))

dat <- dat_anthropometry.no_dropout.no_flood.trt
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

stat_summary <- data.frame()

for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_phase_only",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}
  
stat_summary

# ranova - test for importance of random effects
# ls_means - calculate and compare least squares means
# step - stepwise simplification

# trial anthropometry responses by treatment group (ITT) -------------------------------------------

dat <- dat_anthropometry.trt
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

stat_summary <- data.frame()

for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_phase_only",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}

stat_summary

# trial plus follow-up by treatment group ----------------------------------------------------

dat_anthropometry.no_dropout.no_flood <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment"))

dat <- dat_anthropometry.no_dropout.no_flood
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_plus_followup",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}

# trial plus follow-up by treatment group (ITT) ----------------------------------------------------

dat <- dat_anthropometry.trt_pl_followup
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

stat_summary <- data.frame()
for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_plus_followup",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}

stat_summary %>%
  arrange(term)

# trial anthropometry responses by treatment group (all participants) -------------------------------------------

dat <- dat_anthropometry.no_dropout.trt
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

stat_summary <- data.frame()

for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_phase_only",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}

stat_summary

# trial plus follow-up by treatment group (all participants) ----------------------------------------------------

dat <- dat_anthropometry.no_dropout.trt_pl_followup
treatment_groups <- c("MDCF-2", "RUSF")
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("study_week")
anova_terms <- c("study_week")

for (treatment_group in treatment_groups) {
  dat_sub <- dat %>%
    filter(group == treatment_group)
  for (dep_var in anthro_metrics) {
    for (i in 1:length(model_terms)) {
      mod <- lmer(dat_sub[[dep_var]] ~ gender + area + study_week + area:study_week + (1|pid), data = dat_sub)
      anova <- stats::anova(mod, ddf="Kenward-Roger")
      stat_summary <- stat_summary %>%
        bind_rows(data.frame(data = "treatment_plus_followup",
                             group = treatment_group,
                             metric = dep_var,
                             term = model_terms[i],
                             coef = get_coef(mod, model_terms[i]),
                             lci = get_confint(mod, model_terms[i])[1],
                             uci = get_confint(mod, model_terms[i])[2],
                             pval = get_anova_p(anova, anova_terms[i])))
    }
  }
}

# write.table(stat_summary, file = paste0(datestring, "_", study, "_anthropometry_analysis_summary_by_group_flood_included.txt"), row.names = FALSE, sep = "\t")

# Relating enrollment anthropometry to length of acute treatment ----------

# Changing gears, here, to look at anthropometry responses within the acute rehabilitation phase
# of the study, prior to MAM treatment with MDCF-2/RUSF

# subset to the correct dataset between enrollment and MAM treatment
dat_anthropometry.no_dropout.acute <- dat_anthropometry.no_dropout %>%
  filter(study_phase %in% c("enrollment", "baseline") & flood_status_jan == "unaffected") %>%
  select(-flood_status_original, -flood_status_jan, -flood_status_may, -dropout_status, -study_week) %>%
  pivot_wider(id_cols = c(pid, area, group, gender), names_from = study_phase, values_from = c(wlz, waz, laz, muac, agedays)) %>%
  mutate(acute_rehab_days = agedays_baseline - agedays_enrollment)

# summary of length of time spent in acute rehab
dat_anthropometry.no_dropout.acute %>%
  group_by(area) %>%
  summarize(range(acute_rehab_days))

# density of enrollment and baseline WLZ scores
ggplot(dat_anthropometry.no_dropout.acute) +
  geom_density(aes(x = wlz_enrollment), fill = "red", alpha = 0.5) +
  geom_density(aes(x = wlz_baseline), fill = "blue", alpha = 0.5) +
  theme_classic()

# There's a lot of continuous variation of anthropometry scores. This exploratory 
# analysis seeks to break down the severity of wasting further into categories or quartiles

# The research question here is whether the severity of wasting affects the trajectory of recovery
# for each participant
severity_quartiles <- dat_anthropometry.no_dropout %>%
  filter(study_phase %in% "enrollment" & flood_status_jan == "unaffected") %>%
  select(pid, wlz, area) %>%
  mutate(wlz_severity = cut(wlz, breaks = c(-Inf, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5)),
         wlz_severity_q = cut(wlz, breaks = quantile(wlz), include.lowest = TRUE)) %>%
  group_by(area) %>%
  mutate(wlz_severity_area = cut(wlz, breaks = c(-Inf, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5)),
         wlz_severity_q_area = cut(wlz, breaks = quantile(wlz), include.lowest = TRUE)) %>%
  ungroup() %>%
  select(-area)

# preparing to do some plotting via this strategy
dat_anthropometry.no_dropout.acute.plot <- dat_anthropometry.no_dropout %>%
  filter(study_phase %in% c("enrollment", "baseline") & flood_status_jan == "unaffected") %>%
  select(-flood_status_original, -flood_status_jan, -flood_status_may, -dropout_status, -study_week) %>%
  group_by(pid) %>%
  mutate(days_acute_treatment = agedays - min(agedays)) %>%
  left_join(severity_quartiles %>% select(-wlz), by = "pid")
  
ggplot(dat_anthropometry.no_dropout.acute.plot, aes(x = agedays, y = wlz, group = pid)) +
  geom_line()

ggplot(dat_anthropometry.no_dropout.acute.plot, aes(x = days_acute_treatment, y = wlz, group = pid)) +
  geom_line()

ggplot(dat_anthropometry.no_dropout.acute.plot, aes(x = days_acute_treatment, y = wlz, group = pid)) +
  geom_line() +
  facet_grid(. ~ wlz_severity)

ggplot(dat_anthropometry.no_dropout.acute.plot %>% filter(study_phase == "baseline"), aes(y = days_acute_treatment, x = wlz_severity, fill = wlz_severity)) +
  geom_boxplot()

# summaries of the categorizations by predefined categories
dat_anthropometry.no_dropout.acute.plot %>%
  filter(study_phase == "baseline") %>%
  group_by(wlz_severity) %>%
  summarize(mean = mean(days_acute_treatment),
            median = median(days_acute_treatment),
            n = length(days_acute_treatment))

# baseline plots
ggplot(dat_anthropometry.no_dropout.acute.plot %>% filter(study_phase == "baseline"), aes(y = days_acute_treatment, x = wlz_severity_q, fill = wlz_severity_q)) +
  geom_boxplot()

# summaries of the categorizations by quartiles
dat_anthropometry.no_dropout.acute.plot %>%
  filter(study_phase == "baseline") %>%
  group_by(wlz_severity_q) %>%
  summarize(mean = mean(days_acute_treatment),
            median = median(days_acute_treatment),
            n = length(days_acute_treatment))

# preparing for plotting
dat_anthropometry.no_dropout.acute.plot2 <- dat_anthropometry.no_dropout.acute.plot %>%
  mutate(days_acute_treatment2 = max(agedays) - min(agedays))

ggplot(dat_anthropometry.no_dropout.acute.plot2 %>% filter(study_phase == "enrollment"), aes(x = days_acute_treatment2, y = wlz)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(dat_anthropometry.no_dropout.acute.plot2 %>% filter(study_phase == "enrollment"), aes(x = days_acute_treatment2, y = wlz, group = area)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(. ~ area, scales = "free")

dat_anthropometry.no_dropout.no_flood.enrollment.cortest <- dat_anthropometry.no_dropout.no_flood.enrollment %>%
  left_join(dat_anthropometry.no_dropout.acute.plot %>% select(pid, days_acute_treatment, wlz_severity, wlz_severity_q, wlz_severity_area, wlz_severity_q_area) %>% distinct(), by = "pid")

# preparing for a statistical test, here. How are the data distributed?
hist(dat_anthropometry.no_dropout.no_flood.enrollment.cortest$wlz)
hist(dat_anthropometry.no_dropout.no_flood.enrollment.cortest$days_acute_treatment)

# Is there a correlation between enrollment WLZ and days required to exit acute rehabilitation?

# across the whole trial?
cor.test(dat_anthropometry.no_dropout.no_flood.enrollment.cortest[dat_anthropometry.no_dropout.no_flood.enrollment.cortest$days_acute_treatment != 0,]$wlz,dat_anthropometry.no_dropout.no_flood.enrollment.cortest[dat_anthropometry.no_dropout.no_flood.enrollment.cortest$days_acute_treatment != 0,]$days_acute_treatment, method = 'spearman')
         

# within the Dhaka site?
cor.test(dat_anthropometry.no_dropout.no_flood.enrollment.cortest %>% filter(area == "Dhaka") %>% pull(wlz), 
         dat_anthropometry.no_dropout.no_flood.enrollment.cortest %>% filter(area == "Dhaka") %>% pull(days_acute_treatment), method = "spearman")

# within the Kurigram site?
cor.test(dat_anthropometry.no_dropout.no_flood.enrollment.cortest %>% filter(area == "Kurigram") %>% pull(wlz), 
         dat_anthropometry.no_dropout.no_flood.enrollment.cortest %>% filter(area == "Kurigram") %>% pull(days_acute_treatment), method = "spearman")

# using an alternative statistical strategy that allows covariates
# here we test the relationship between time in acute rehab and wlz, severity, or severity quartile

mod.enrollment <- lm(days_acute_treatment ~ wlz + gender + area, data = dat_anthropometry.no_dropout.no_flood.enrollment.cortest)
anova(mod.enrollment)

mod.enrollment2 <- lm(days_acute_treatment ~ wlz_severity + gender + area, data = dat_anthropometry.no_dropout.no_flood.enrollment.cortest)
anova(mod.enrollment2)

mod.enrollment3 <- lm(days_acute_treatment ~ wlz_severity_q + gender + area, data = dat_anthropometry.no_dropout.no_flood.enrollment.cortest)
anova(mod.enrollment3)

# Summary of SAM recovery trajectories ------------------------------------

# Changing gears again, here. I know the fecal samples were collected at different timepoints than
# the anthropometry data. To relate the two, I need some summary of trajectories for each.

# calculating simple rates/velocities to use with the microbiota data
dat_anthropometry.no_dropout.acute.rates <- dat_anthropometry.no_dropout.acute %>%
  mutate(wlz_rate = (wlz_baseline - wlz_enrollment) / acute_rehab_days,
         waz_rate = (waz_baseline - waz_enrollment) / acute_rehab_days,
         laz_rate = (laz_baseline - laz_enrollment) / acute_rehab_days,
         muac_rate = (muac_baseline - muac_enrollment) / acute_rehab_days) %>%
  select(pid, area, group, gender, wlz_rate, waz_rate, laz_rate, muac_rate) %>% 
  pivot_longer(cols = c(-pid, -area, -group, -gender), names_to = "metric") %>%
  left_join(severity_quartiles, by = "pid")

# plotting the results
ggplot(dat_anthropometry.no_dropout.acute.rates, aes(x = paste(group, area), y = value, fill = metric)) +
  geom_hline(yintercept = 0, lty = 2, color = "blue") +
  geom_boxplot() +
  facet_wrap( ~ metric, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))

ggplot(dat_anthropometry.no_dropout.acute.rates %>% filter(metric == "wlz_rate"), aes(x = wlz, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

# looking at enrollment anthropometry versus time required for acute rehab
dat_anthropometry.no_dropout.acute.rates.cortest <- dat_anthropometry.no_dropout.acute.rates %>% 
  filter(metric == "wlz_rate")

mod.enrollment4 <- lm(value ~ wlz + gender + area, data = dat_anthropometry.no_dropout.acute.rates.cortest)
anova(mod.enrollment4)

get_coef(mod.enrollment4, "wlz")
get_coef(mod.enrollment4, "area")

# write.table(dat_anthropometry.no_dropout.acute.rates, paste0(datestring, "_", study, "_enrollment_anthro_vs_rate_of_change.txt"), sep = "\t", row.names = FALSE)

# based on the above analysis, I see a trend toward a relationship between starting WLZ and time required for acute rehabilitation but it's not significant
## even when controlled for gender, area, etc. So, I suppose this leaves the door open to microbial configurations that are associated with or predictors of
## rapidity of response. It also gets to thinking about enrollment versus baseline configurations and how different they are along the spectrum of 
## SAM severity.

# Dhaka anthropometry responses -------------------------------------------

# Returning to analyses similar to those earlier in this document, but now looking at recovery
# within each study site
dat_anthropometry.no_dropout.no_flood.trt.dhaka <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment", "followup") & area == "Dhaka")
length(unique(dat_anthropometry.no_dropout.no_flood.trt.dhaka$pid))

dat <- dat_anthropometry.no_dropout.no_flood.trt.dhaka
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week")
anova_terms <- c("group:study_week")

stat_summary.site <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + group + study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary.site <- stat_summary.site %>%
      bind_rows(data.frame(data = "Dhaka_treatment",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

dat_anthropometry.no_dropout.no_flood.dhaka <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment") & area == "Dhaka")
length(unique(dat_anthropometry.no_dropout.no_flood.dhaka$pid))

dat <- dat_anthropometry.no_dropout.no_flood.dhaka
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week")
anova_terms <- c("group:study_week")

#stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + group + study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary.site <- stat_summary.site %>%
      bind_rows(data.frame(data = "Dhaka_treatment_plus_followup",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

# Kurigram anthropometry responses -------------------------------------------

dat_anthropometry.no_dropout.no_flood.trt.kurigram <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment", "followup") & area == "Kurigram")
length(unique(dat_anthropometry.no_dropout.no_flood.trt.kurigram$pid))

dat <- dat_anthropometry.no_dropout.no_flood.trt.kurigram
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
#model_terms <- c("gender", "groupMDCF-2", "study_week", "groupMDCF-2:study_week")
#anova_terms <- c("gender", "group", "study_week", "group:study_week")
model_terms <- c("groupMDCF-2:study_week")
anova_terms <- c("group:study_week")

#stat_summary.site <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
#    mod <- lmer(dat[[dep_var]] ~ gender + group + group*feeding_start_month + study_week + group*study_week + (1|pid), data = dat)
    mod <- lmer(dat[[dep_var]] ~ gender + group + study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary.site <- stat_summary.site %>%
      bind_rows(data.frame(data = "Kurigram_treatment",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

dat_anthropometry.no_dropout.no_flood.kurigram <- dat_anthropometry.no_dropout.no_flood %>%
  filter(!study_phase %in% c("enrollment") & area == "Kurigram")
length(unique(dat_anthropometry.no_dropout.no_flood.kurigram$pid))

dat <- dat_anthropometry.no_dropout.no_flood.kurigram
anthro_metrics <- c("wlz", "waz", "laz", "muac", "weight", "length")
model_terms <- c("groupMDCF-2:study_week")
anova_terms <- c("group:study_week")

#stat_summary <- data.frame()
for (dep_var in anthro_metrics) {
  for (i in 1:length(model_terms)) {
    mod <- lmer(dat[[dep_var]] ~ gender + group + study_week + group*study_week + (1|pid), data = dat)
    anova <- stats::anova(mod, ddf="Kenward-Roger")
    stat_summary.site <- stat_summary.site %>%
      bind_rows(data.frame(data = "Kurigram_treatment_plus_followup",
                           metric = dep_var,
                           term = model_terms[i],
                           coef = get_coef(mod, model_terms[i]),
                           lci = get_confint(mod, model_terms[i])[1],
                           uci = get_confint(mod, model_terms[i])[2],
                           pval = get_anova_p(anova, anova_terms[i])))
  }
}

# write.table(stat_summary, file = paste0(datestring, "_", study, "_anthropometry_analysis_summary.txt"), row.names = FALSE, sep = "\t")
# write.table(stat_summary.site, file = paste0(datestring, "_", study, "_anthropometry_analysis_summary_by_site.txt"), row.names = FALSE, sep = "\t")

# Summary plots for results of anthropometry analyses ---------------------

stat_summary.agg <- stat_summary %>%
  filter(term == "groupMDCF-2:study_week") %>%
  bind_rows(stat_summary.site)

# pdf(paste0(datestring, "_", study, "_anthropometry_effects.pdf"))
ggplot(stat_summary.agg %>% filter(term == "groupMDCF-2:study_week" & !metric %in% c("weight", "length")), aes(x = coef, y = metric)) +
  geom_vline(xintercept = 0, alpha = 0.5, color = "blue", lty = 2) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.25) +
  geom_point(aes(fill = term), size = 2, pch = 21) +
  facet_grid( ~ data) +
  theme_bw()
# dev.off()

##

ggplot(dat_anthropometry.no_dropout.no_flood.trt, aes(x = study_week, y = wlz, group = group, color = group)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic()

ggplot(dat_anthropometry.no_dropout.flood.trt, aes(x = study_week, y = wlz, group = group, color = group)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic()

ggplot(dat_anthropometry.no_dropout.no_flood.trt, aes(x = study_week, y = wlz, group = paste(group, study_week, sep = "_"), fill = group)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic()

ggplot(dat_anthropometry.no_dropout.no_flood.trt, aes(x = study_week, y = wlz, group = group, color = group)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic() +
  facet_grid(~ area)

ggplot(dat_anthropometry.no_dropout.trt, aes(x = study_week, y = wlz, group = group, color = group)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic() +
  facet_grid(~ area)

ggplot(dat_anthropometry.no_dropout.no_flood.trt, aes(x = study_week, y = wlz, group = paste(group, study_week, sep = "_"), fill = group)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic() +
  facet_grid(~ area)

ggplot(dat_anthropometry.no_dropout, aes(x = study_week, y = wlz, group = paste(group, study_week, sep = "_"), fill = group)) +
  geom_boxplot() +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  theme_classic() +
  facet_grid(~ area)

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

