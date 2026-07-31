library(CMA)
library(dplyr)
library(tidyr)
library(ggplot2)

#########################################################################
# Generate input files. Each file is from a unique study and needs:
# - One column called 'markname', which will have aptamer IDs)
# - Another column called 'pval' which has unadjusted p-values
#########################################################################
treatment_association <- read.csv('treatment_assocation_primaryMAM/supplemental_table.csv')

treatment_association %>% 
  filter(pM_armMDCF2_pval != "  NA") %>% 
  select(SeqId, pM_armMDCF2_pval) %>% 
  rename(markname = 'SeqId',
         pval = 'pM_armMDCF2_pval') %>% 
  write.csv('raw_data/pM_tx_res.csv')

treatment_association %>% 
  filter(pM_armMDCF2_pval != "  NA") %>% 
  select(SeqId, pSM_armMDCF2_pval) %>% 
  rename(markname = 'SeqId',
         pval = 'pSM_armMDCF2_pval') %>% 
  write.csv('raw_data/pSM_MAM_tx_res.csv')

#########################################################################
# Vaha's CMA code (untouched besides input/output dir update)
#########################################################################

# Specify the directory containing your CSV files
input_dir <- "raw_data/"

# List all CSV files in the directory
input_files <- list.files(input_dir, pattern = ".csv$", full.names = TRUE)

output_dir <- '001_CMA_res/'

# Extract the variable names from the file names
input_vars <- lapply(input_files, function(x) strsplit(basename(x), "\\.")[[1]][1]) %>% unlist

# Read only relevant columns from files "markname,pval"
df_input_list <- lapply(input_files, function(x) select(read.csv(x), markname, pval) )

# Create a vector that counts how many files each markname appears in
markname_counts <- table(unlist(lapply(df_input_list, function(x) x$markname)))

# Filter dataframes to only include marknames that appear in more than one file
df_input_list <- lapply(df_input_list, function(x) filter(x, markname %in% names(markname_counts)[markname_counts > 1]))

# Merge by markname
df_input <- Reduce(function(x,y) merge(x,y, by='markname', all=TRUE), df_input_list)

# Remove duplicates
df_input <- df_input[!duplicated(df_input),]

colnames(df_input)[-1] = input_vars
write.csv(df_input, paste0(output_dir, "merged_input.csv"), row.names = FALSE)


# ------------- CMA --------------- #
# calculate tetrachoric correlations
result_tetracorr <- CMA::tetracorr(df_input, input_vars)

# write tetrachoric correlations to disk
# create output directory if it doesn't exist
dir.create(output_dir, showWarnings=FALSE)
write.table(result_tetracorr$sigma, paste0(output_dir, "tetrachor_sigma.txt"), quote=FALSE, sep="\t", row.names=FALSE, na="NA")


# calcaulte fisher's p-value
result_fisher <- CMA::fishp(df_input, input_vars, result_tetracorr$sigma, result_tetracorr$sum_sigma)
#print(paste0("# fisher result nrows: ", nrow(result_fisher)))
#print(paste0("# fisher result nrows (unique): ", nrow(result_fisher[!duplicated(result_fisher),])))
#print(head(result_fisher))

# merge final output file
result_fisher$corr <- (result_fisher$sum_sigma_var - 2) / 2

# calculate nlog10p for input columns
result_fisher <- result_fisher %>%
    mutate(across(all_of(input_vars), ~ -log10(.), .names = "{.col}_nlogp"))

# subset final output file columns
CMA_output <- result_fisher %>%
    select(c("markname","meta_nlog10p","meta_p"), contains(input_vars, ignore.case=TRUE), "corr") %>%
    arrange(meta_p)
print(paste0("# CMA_output result nrows: ", nrow(CMA_output)))
print(head(CMA_output))

# write CMA result to disk
write.csv(CMA_output, paste0(output_dir, "CMA_meta.csv"), quote=FALSE, row.names=FALSE, na="NA")

print("DONE")


#########################################################################
# SJH: Annotate the output, add FDR (Benjamini Hochberg), and original values compared
#########################################################################
# load each file of interest
CMA_output <- read.csv('001_CMA_res/CMA_meta.csv')
pSM_MAM_pM <- read.csv('treatment_assocation_primaryMAM/supplemental_table.csv')

#####  Create p-val dist
pSM_MAM_pM %>% 
  filter(pM_armMDCF2_pval != "  NA") %>% 
  ggplot(aes(x = as.numeric(pM_armMDCF2_pval))) + 
  geom_histogram() + 
  ggtitle(paste0('primary-MAM p-value distribution\nnum of markers: ', nrow(pSM_MAM_pM %>% 
                                                                              filter(pM_armMDCF2_pval != "  NA")))) + 
  theme_bw()
ggsave(paste0(output_dir, '/histogram_pM_pval_dist.png'))

pSM_MAM_pM %>% 
  ggplot(aes(x = as.numeric(pSM_armMDCF2_pval))) + 
  geom_histogram() + 
  ggtitle(paste0('post-SAM MAM p-value distribution\nnum of markers: ', nrow(pSM_MAM_pM))) + 
  theme_bw()
ggsave(paste0(output_dir, '/histogram_pSM_MAM_pval_dist.png'))

# Find the marks consistently positive across the datasets
pos_marks <- intersect(pSM_MAM_pM[pSM_MAM_pM$pSM_armMDCF2 > 0, ]$SeqId, pSM_MAM_pM[pSM_MAM_pM$pM_armMDCF2 > 0, ]$SeqId)

# Find the marks consistently negative across the datasets
neg_marks <- intersect(pSM_MAM_pM[pSM_MAM_pM$pSM_armMDCF2 < 0, ]$SeqId, pSM_MAM_pM[pSM_MAM_pM$pM_armMDCF2 < 0, ]$SeqId)


soma_anno <- read.csv('CMA/soma_anno.csv')

soma_anno_full <- read.csv('soma_anno.csv')

### Analysis
# - Apply FDR (BH and Bonferonni for comparison)
# - Join soma annotations
# - Remove duplicate targets (Retain the most significant)
# - Join the original fit results for each protein
annotated <- CMA_output %>% 
  mutate(meta_pval_BH = p.adjust(meta_p, method = 'BH'),
         meta_pval_Bonferroni = p.adjust(meta_p, method = 'bonferroni')) %>% 
  left_join(pSM_MAM_pM %>% 
              dplyr::rename(markname = 'SeqId',
                     pSM_treatment_arm_coeff = 'pSM_armMDCF2',
                     pSM_treatment_arm_pval = 'pSM_armMDCF2_pval',
                     pM_treatment_arm_coeff = 'pM_armMDCF2',
                     pM_treatment_arm_pval = 'pM_armMDCF2_pval') %>% 
              dplyr::select(markname, SeqIdVersion, EntrezGeneSymbol, TargetFullName, pSM_treatment_arm_coeff, pSM_treatment_arm_pval, pM_treatment_arm_coeff, pM_treatment_arm_pval)) %>%  
  mutate(pM_direction = ifelse(pM_treatment_arm_coeff > 0, 'pos', 'neg'),
         pSM_direction = ifelse(pSM_treatment_arm_coeff > 0, 'pos', 'neg')) %>% 
  mutate(direction = paste0(pSM_direction,'/', pM_direction)) %>% 
  left_join(soma_anno_full %>% dplyr::select(SeqId, EntrezGeneID) %>% dplyr::rename(markname = 'SeqId'))

write.csv(annotated, '001_CMA_res/cma_fdr_res.csv')

#########################################################################
# GSEA of WLZ associated proteins groups in the MDCF-2 or RUSF treatment arms
#########################################################################
library(fgsea)
library(data.table)

# Load the WLZ-associated CMA results
wlz_cma <- read.csv('003_module/cma_res_supplemental.csv') %>%
  filter(direction_MAMphase_pM %in% c('pos/pos', 'neg/neg')) %>% 
  mutate(meta_pval_BH = as.numeric(meta_pval_BH)) %>% 
  filter(meta_pval_BH < 0.05)

treatment_arm_CMA <- read.csv('001_CMA_res/CMA_meta.csv')
pSM_MAM_pM <- read.csv('supplemental_table.csv')

soma_anno_full <- read.csv('soma_anno.csv')

treatment_CMA_direction <- treatment_arm_CMA %>% 
  mutate(meta_pval_BH = p.adjust(meta_p, method = 'BH'),
         meta_pval_Bonferroni = p.adjust(meta_p, method = 'bonferroni')) %>% 
  left_join(pSM_MAM_pM %>% 
              dplyr::rename(markname = 'SeqId',
                            pSM_treatment_arm_coeff = 'pSM_armMDCF2',
                            pSM_treatment_arm_pval = 'pSM_armMDCF2_pval',
                            pM_treatment_arm_coeff = 'pM_armMDCF2',
                            pM_treatment_arm_pval = 'pM_armMDCF2_pval') %>% 
              dplyr::select(markname, SeqIdVersion, EntrezGeneSymbol, TargetFullName, pSM_treatment_arm_coeff, pSM_treatment_arm_pval, pM_treatment_arm_coeff, pM_treatment_arm_pval)) %>%  
  mutate(pM_direction = ifelse(pM_treatment_arm_coeff > 0, 'pos', 'neg'),
         pSM_direction = ifelse(pSM_treatment_arm_coeff > 0, 'pos', 'neg')) %>% 
  mutate(direction = paste0(pSM_direction,'/', pM_direction)) %>% 
  left_join(soma_anno_full %>% dplyr::select(SeqId, EntrezGeneID) %>% dplyr::rename(markname = 'SeqId'))%>% 
  filter(direction %in% c('pos/pos', 'neg/neg')) %>% 
  mutate(rank = ifelse(direction == 'pos/pos', meta_nlog10p, -1*meta_nlog10p))

soma_pathways <- list(pos_soma = wlz_cma[wlz_cma$direction_MAMphase_pM == 'pos/pos',]$markname,
                 neg_soma = wlz_cma[wlz_cma$direction_MAMphase_pM == 'neg/neg',]$markname)

soma_ranks <- treatment_CMA_direction$rank
names(soma_ranks) <- treatment_CMA_direction$markname
soma_ranks <- sort(soma_ranks)

gsea_res <- fgsea(pathways = soma_pathways,
      stats = soma_ranks) %>% 
  mutate(leadingEdge = sapply(leadingEdge, toString))

write.csv(gsea_res, '002_CMA_treatment_GSEA/WLZ_proteins_ranked_by_treamtent_assocation.csv')

plotEnrichment(soma_pathways$pos_soma, soma_ranks)
plotEnrichment(soma_pathways$neg_soma, soma_ranks)

library(stringr)
x <- as.data.frame(gsea_res)

leading_edge_interp <- data.frame(SeqId = c(str_split(x$leadingEdge[[2]], pattern = ', ')[[1]], str_split(x$leadingEdge[[1]], pattern = ', ')[[1]])) %>% 
  left_join(soma_anno_full %>% 
              select(SeqId, SeqIdVersion, EntrezGeneSymbol, TargetFullName) %>% 
              unique()) 

write.csv(leading_edge_interp, '002_CMA_treatment_GSEA/leading_edges.csv')

library(enrichplot)
library(clusterProfiler)
test <- GSEA(geneList = sort(soma_ranks, decreasing = T),
             TERM2GENE = select(wlz_cma, direction_MAMphase_pM, markname),
             pvalueCutoff = 1)


gseaplot2(test, geneSetID = 1:2)




