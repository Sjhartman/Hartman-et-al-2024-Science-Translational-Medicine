# Purpose: QC SOMAmers by the following criteria:
#   1. Collect summary stats of plates and controls
# 2. Exclude flagged samples and document noted samples
# 3. Exclude non-human aptamers
# 4. Exclude non-protein targeted aptamers
# 5. Exclude SOMAmers that are flagged on more than one plate
# 6. Exclude low-activity aptamers
# 7. Save adat


### Libraries
library(SomaDataIO)
library(tidyverse)


prjDir <- ""
adatDir <- "OriginalAdatFiles"
runDir <- ""
studyID <- "SS-216868"
my_adat <- read_adat(paste0(prjDir, "/", adatDir, "/SS-216868_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"))
my_adat <- read_adat('/Users/stevenhartman/Library/CloudStorage/Box-Box/000Gordon/Projects/PostSamMam/SOMALogic/WUS-216868_v4.1_EDTAPlasma_20210813/OriginalAdatFiles/SS-216868_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat')


### Filter flood unaffected and non-dropout
# Get patient IDs in format used for filtering MAGs
PatientId.df <-  filter(my_adat, SampleType == "Sample") %>% 
  select(SampleId) %>% 
  mutate(PatiendId = str_extract(SampleId, "^[:alnum:]{7}"))
PatientId.df$PatiendId <- gsub("F", "C", PatientId.df$PatiendId)
# Add patient IDs to adat
my_adat <- left_join(my_adat, PatientId.df, by = "SampleId")
# Load unaffected list
floodUnaffected.NoDropout <- read.table("/Users/stevenhartman/Library/CloudStorage/Box-Box/000Gordon/Projects/PostSamMam/Data/PID_Flood_Unaffected_No_Dropout.txt")$V1
# Use unaffected list to create affected list
floodAffected_Dropout <- setdiff(PatientId.df$PatiendId, floodUnaffected.NoDropout)
# Filter adat
my_adat <- my_adat %>% 
  filter(!PatiendId %in% floodAffected_Dropout)


#################################################
###   Step 1: Collect summary stats from the run
################################################
# - Total number of samples
# - Number of plates (96 well format)
# - Number of plasma samples, calibrators, QCs and buffer samples / plate (df)


print(paste0("Number of patients: ", length(unique(my_adat$PatiendId))-1))

paste0("Total number of control samples: ", nrow(filter(my_adat,
                                                        SampleType != "Sample")))
paste0("Total number of samples (sample only): ", nrow(filter(my_adat,
                                                              SampleType == "Sample")))
paste0("Number of 96-well plates: ", length(unique(my_adat$PlateId)))
# Number of plasma samples, calibrators, QCs and buffer samples / plate (df)
plate.summary.df <- my_adat %>% count(PlateId, SampleType) %>% 
  spread(SampleType, n) %>% 
  mutate(Total = Buffer + Calibrator + QC + Sample)

#################################################
###   Step 2: Exclude flagged samples and document noted samples
################################################
my_adat$PID <-substr(my_adat$SampleId, 1, 7)
my_adat_flood_unaffected <- filter(my_adat, PID %in% floodUnaffected.NoDropout | substr(PID, 1, 1) != "S" )

# Extract flagged samples & determine reason for flag
flaggedTagged <- my_adat_flood_unaffected  %>% 
  select("SampleId", "SampleType", "AssayNotes", "HybControlNormScale", "NormScale_20", "NormScale_0_5", "NormScale_0_005", "RowCheck") %>% 
  filter(RowCheck == "FLAG" | !is.na(AssayNotes)) %>% 
  mutate(Action = ifelse(NormScale_20 > 0.4 & NormScale_20 < 2.5, 
                         ifelse(NormScale_0_5 > 0.4 & NormScale_0_5 < 2.5,
                                ifelse(NormScale_0_005 > 0.4 & NormScale_0_005 < 2.5, 
                                       ifelse(HybControlNormScale > 0.4 & HybControlNormScale < 2.5, "Retain", 
                                              "HybControlNormScale_OutOfRange"), 
                                       "NormScale_0_005_OutOfRange"), 
                                "NormScale_0_5_OutOfRange"), 
                         "NormScale_20_OutOfRange")
  ) %>% 
  arrange(RowCheck, desc(Action))
flaggedTagged

# Save the flagged samples as csv
# write.csv(flaggedTagged, paste0(prjDir, "/", runDir, "/step2.", studyID, "_FlaggedTaggedSmpls.csv"))

# Summary stats
paste0("Num of flagged samples: ", sum(my_adat_flood_unaffected$RowCheck == "FLAG"))
paste0("Num of flagged cntrl samples: ", 
       sum(my_adat_flood_unaffected$RowCheck == "FLAG" & my_adat_flood_unaffected$SampleType != "Sample"))
paste0("Num of noted samples: ", 
       sum(!is.na(my_adat_flood_unaffected$AssayNotes) & my_adat_flood_unaffected$SampleType == "Sample"))
paste0("Num of noted cntrls: ", 
       sum(!is.na(my_adat_flood_unaffected$AssayNotes) & my_adat_flood_unaffected$SampleType != "Sample"))
paste0("Num of samples both flagged and noted: ", sum(my_adat_flood_unaffected$RowCheck == "FLAG" & 
                                                        !is.na(my_adat_flood_unaffected$AssayNotes) &
                                                        my_adat_flood_unaffected$SampleType == "Sample"))
paste0("Num of cntrls both flagged and noted: ", sum(my_adat_flood_unaffected$RowCheck == "FLAG" & 
                                                       !is.na(my_adat_flood_unaffected$AssayNotes) &
                                                       my_adat_flood_unaffected$SampleType != "Sample"))
paste0("Passing samples: ", 
       sum(my_adat_flood_unaffected$RowCheck == "PASS" & my_adat_flood_unaffected$SampleType == "Sample"))
paste0("Passing cntrls: ", 
       sum(my_adat_flood_unaffected$RowCheck == "PASS" & my_adat_flood_unaffected$SampleType != "Sample"))

paste0("Percent passing samples: ", sum(my_adat_flood_unaffected$RowCheck == "FLAG")/nrow(my_adat_flood_unaffected) * 100)

# Filter out the flagged samples
my_adat_flood_unaffected_normFilt <- my_adat_flood_unaffected %>% 
  filter(!SampleId %in% flaggedTagged$SampleId[flaggedTagged$RowCheck == "FLAG"])


#################################################
###   Step 3 Remove non-human aptamers, create summary table of non-human samples
################################################
seq.attribute.df <- attr(my_adat_flood_unaffected_normFilt, "Col.Meta")
organsimIds <- unique(seq.attribute.df$Organism)

# df with SOMAmers per organism
organism.summary.df <- count(seq.attribute.df, Organism) %>% 
  arrange(desc(n))
# write.csv(organism.summary.df, 
#           paste0(prjDir, "/", runDir, "/step3.", studyID, ".SOMAmersPerOrganism.csv"))

# df with SOMAmer counts per organism and target
organism.target.summary.df <- count(seq.attribute.df, Organism, Target) %>% 
  arrange(desc(n)) %>% 
  arrange(Organism, desc(n))

# write.csv(organism.target.summary.df, 
#           paste0(prjDir, "/", runDir, "/step3.", studyID, ".SOMAmersPerOrganismAndTarget.csv"))

# Create df with only human-targeting SOMAmers
nonHumanseq_attr.df <- seq.attribute.df %>% 
  filter(Organism != "Human")
Humanseq_attr.df <- seq.attribute.df %>% 
  filter(Organism == "Human")
my_adat_flood_unaffected_normHumanFilt <- my_adat_flood_unaffected_normFilt %>% 
  select(!seqid2apt(nonHumanseq_attr.df$SeqId))
# write.table(nonHumanseq_attr.df$SeqId, paste0(prjDir, "/", runDir, "/step3.", studyID, ".nonHumanTrgtAptamer.txt"),
#             row.names = F,
#             col.names = F,
#             quote = F)

# Summary stats
paste0("Total num SOMAmers: ", nrow(seq.attribute.df))
paste0("Total human SOMAmers: ", nrow(Humanseq_attr.df))
paste0("Total non-human SOMAmers: ", nrow(nonHumanseq_attr.df))


#################################################
### Step 4 Summarize and exclude non-protein targetting aptamers
################################################
seq.attribute.Human.only.df <- attr(my_adat_flood_unaffected_normHumanFilt, "Col.Meta")
# df of SOMAmer counts for different target types 
TargetType.summary.df <- count(seq.attribute.Human.only.df, Type) %>% 
  arrange(desc(n))
# write.csv(TargetType.summary.df, 
#           paste0(prjDir, "/", runDir, "/step4.", studyID, ".Human.SOMAmers.TargetTypes.csv"))

# df with SOMAmer counts per target type and target
TargetType.and.target.non_protein.summary.df <- count(filter(seq.attribute.Human.only.df,
                                                             Type != "Protein"), 
                                                      Type, Target, TargetFullName)
# write.csv(TargetType.and.target.non_protein.summary.df, 
#           paste0(prjDir, "/", runDir, "/step4.", studyID, ".Human.nonProtein.SOMAmers.TargetTypesAndTarget.csv"))

# filter out non-protein targeting aptamers
seq.attribute.HumanAndProtein.only.df <- filter(seq.attribute.Human.only.df, Type == "Protein")
seq.attribute.HumanAndnonProtein.only.df <- filter(seq.attribute.Human.only.df, Type != "Protein")
my_adat_flood_unaffected_normHuman.protein.Filt <- my_adat_flood_unaffected_normHumanFilt %>% 
  select(!seqid2apt(seq.attribute.HumanAndnonProtein.only.df$SeqId))

# Summary stats
paste0("Total num SOMAmers: ", nrow(seq.attribute.Human.only.df))
paste0("non-protien, human-targetting SOMAmers: ", sum(seq.attribute.Human.only.df$Type != "Protein"))
paste0("Remaining SOMAmers: ", nrow(seq.attribute.HumanAndProtein.only.df))

#################################################
### Step 5 Exclude SOMAmers that are flagged on more than one plate
################################################
# Assign max and min CalQcRatios
minRatio = 0.8
maxRatio = 1.2

# Create df of only CalQcRatios
seq.attribute.QC.df <- seq.attribute.HumanAndProtein.only.df %>% 
  select(SeqId, starts_with("CalQcRatio")) %>% 
  column_to_rownames("SeqId") 

# Create df of TRUE/FALSE for CalQcRatios within expected range
pass.plate.df =  seq.attribute.QC.df > minRatio & seq.attribute.QC.df < maxRatio

# Count fails per SOMAmer in each plate
FailCounts = as.data.frame(apply(pass.plate.df, 1, function(x) sum(x == FALSE)))
colnames(FailCounts) = "FailCounts"

# Extract SOMAmers that fail on more than 1 plate
FailedSOMAmers <- rownames(filter(FailCounts, FailCounts > 1))
# write.table(FailedSOMAmers, paste0(prjDir, "/", runDir, "/step5.", studyID, ".SOMAmersFlaggedMultiplePLates.txt"),
#             row.names = F,
#             col.names = F,
#             quote = F)

# Reduce seq.attribute.HumanAndProtein.only.df
seq.attribute.HumanAndProtein.only.passing.df <- filter(seq.attribute.HumanAndProtein.only.df,
                                                        !SeqId %in% FailedSOMAmers)
seq.attribute.HumanAndProtein.only.fail.df <- filter(seq.attribute.HumanAndProtein.only.df,
                                                     SeqId %in% FailedSOMAmers)

my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.Filt <- my_adat_flood_unaffected_normHuman.protein.Filt %>% 
  select(!seqid2apt(seq.attribute.HumanAndProtein.only.fail.df$SeqId))

# Summary stats
paste0("Number of flagged SOMAmers: ", sum(seq.attribute.HumanAndProtein.only.df$ColCheck == "FLAG"))
paste0("Number of failed SOMAmers (Flagged in more than 1 plate): ", length(FailedSOMAmers))
paste0("Number of remaining SOMAmers: ", nrow(seq.attribute.HumanAndProtein.only.passing.df))

#################################################
### Step 6 Remove aptamers with mean < buffer mean + 3 buffer SDs
################################################
# Create buffer only df
cntrl.buffer.df <- my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.Filt %>% 
  filter(SampleType == "Buffer") %>% 
  select(starts_with("seq"))

# Creatre sample only df
cntrl.sample.df <- my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.Filt %>% 
  filter(SampleType == "Sample") %>% 
  select(starts_with("seq"))

# Create data summary required low response aptamer identification
cntrl.buffer.smpl.summarys.df <- data.frame(buff.mean = sapply(cntrl.buffer.df, mean),
                                            buff.sd = sapply(cntrl.buffer.df, sd)) %>% 
  mutate("buff.mean_three.sd" = buff.mean + 3*buff.sd) %>% 
  mutate("Sample.mean" = sapply(cntrl.sample.df, mean)) %>% 
  mutate("SampleMean>Thresh" = Sample.mean > buff.mean_three.sd)

# Create summary stats
failedSOMAMers.cnt =  sum(cntrl.buffer.smpl.summarys.df$`SampleMean>Thresh` == FALSE)
passedSOMAMers.cnt = sum(cntrl.buffer.smpl.summarys.df$`SampleMean>Thresh` == TRUE)
lowAbSOMAmers = rownames(cntrl.buffer.smpl.summarys.df)[cntrl.buffer.smpl.summarys.df$`SampleMean>Thresh` == FALSE]
# write.table(lowAbSOMAmers, paste0(prjDir, "/", runDir, "/step6.", studyID, "LowOutputAptamers.txt"),
#             row.names = F,
#             quote = F,
#             col.names = F)

# Summary
paste0("Failed SOMAmers: ", failedSOMAMers.cnt)
paste0("Passing SOMAmers: ", passedSOMAMers.cnt)

# Subset and save adat
my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.lowAb.Filt <- my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.Filt %>% 
  select(!all_of(lowAbSOMAmers))
write_adat(my_adat_flood_unaffected_normHuman.protein.SOMAmerPass.lowAb.Filt,
           paste0('SOMA_out/filtered_', Sys.Date(),  ".adat"))

# Find the number of unique proteins
final_soma <- data.frame(soma = colnames(as.data.frame(final_ab))[34:ncol(as.data.frame(final_ab))])
final_soma_anno <- seq.attribute.df %>% 
  select(SeqId, EntrezGeneSymbol) %>% 
  mutate(soma = paste0('seq.', gsub('-', '.', SeqId))) %>% 
  filter(soma %in% final_soma$soma)


