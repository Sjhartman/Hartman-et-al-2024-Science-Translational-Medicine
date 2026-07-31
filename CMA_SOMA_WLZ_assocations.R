library(CMA)
library(dplyr)
library(ggplot2)

#########################################################################
# Generate input files for CMA
#########################################################################
pM_res <- read.csv('SOMApublished_results_wAptamerIDs/20.04.03 pearcor betaWLZ vs deltaPlasmaP.csv') %>% 
  # mutate(markname = gsub('_', '.', gsub(".[0-9]$", '', gsub('SL.', 'seq.', PlasmaSOMAmer)))) %>% 
  mutate(markname = gsub('_', '.', gsub('SL.', 'seq.', PlasmaSOMAmer))) %>% 
  rename(pval = 'pVal') %>% 
  select(markname, pval, pAdj, pearsonRho, EntrezGeneSymbol)

pSM_res <- read.csv('MAM_lm_res_wlz_soma_area_group_gender.csv') %>% 
  rename(markname = 'SeqId',
         pval = 'soma_deltaMAM_pval_anova') %>% 
  mutate(markname = paste0(markname, '.', SeqIdVersion)) %>% 
  select(markname, soma_deltaMAM, pval, soma_deltaMAM_anova_fdr, EntrezGeneSymbol) %>% 
  filter(markname %in% pM_res$markname)

pM_res <- pM_res %>% 
  filter(markname %in% pSM_res$markname)

dim(pM_res) # 4777
dim(pM_res) # 4777. Note, filtered to all aptamer IDs and versions match

write.csv(pM_res, 'raw_data/pM.csv')

write.csv(pSM_res, 'raw_data/pSM_MAM.csv')

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
library(tidyverse)
library(ggplot2)

# load each file of interest
CMA_output <- read.csv('001_CMA_res/CMA_meta.csv')
pM <- read.csv('raw_data/pM.csv')
pSM_MAM <- read.csv('raw_data/pSM_MAM.csv')

#####  Create p-val dist
ggplot(pM, aes(x = pval)) + 
  geom_histogram() + 
  ggtitle(paste0('primary-MAM p-value distribution\nnum of markers: ', nrow(pM))) + 
  theme_bw()
ggsave(paste0(output_dir, '/histogram_pM_pval_dist.png'))

ggplot(pSM_MAM, aes(x = pval)) + 
  geom_histogram() + 
  ggtitle(paste0('post-SAM MAM p-value distribution\nnum of markers: ', nrow(pSM_MAM))) + 
  theme_bw()
ggsave(paste0(output_dir, '/histogram_pSM_MAM_pval_dist.png'))

# Find the marks consistently positive across the datasets
pos_marks <- intersect(pM[pM$pearson > 0, ]$markname, pSM_MAM[pSM_MAM$soma_deltaMAM > 0, ]$markname)

# Find the marks consistently negative across the datasets
neg_marks <- intersect(pM[pM$pearson < 0, ]$markname, pSM_MAM[pSM_MAM$soma_deltaMAM < 0, ]$markname)

soma_anno <- read.csv('CMA/soma_anno.csv')

soma_anno_full <- read.csv('soma_anno.csv')

### Analysis
# - Filter the outputs to just pos and neg proteins across all conditions
# - Join soma annotations
# - Remove duplicate targets (Retain the most significant)
# - Join the original fit results for each protein
# Apply FDR (BH and Bonferonni for comparison)
annotated <- CMA_output %>% 
  mutate(meta_pval_BH = p.adjust(meta_p, method = 'BH'),
         meta_pval_Bonferroni = p.adjust(meta_p, method = 'bonferroni')) %>% 
  filter(markname %in% c(pos_marks, neg_marks)) %>% 
  mutate(direction = ifelse(markname %in% pos_marks, 'pos', 'neg')) %>% 
  left_join(soma_anno %>% 
              rename(markname = 'SeqId') %>% 
              select(TargetFullName, EntrezGeneSymbol, markname) %>% 
              unique()) %>% 
  left_join(pM %>% 
              rename(pM_coeff = 'pearsonRho',
                     pM_pAdj = 'pAdj',
                     EntrezGeneSymbol_pM = 'EntrezGeneSymbol') %>% 
              select(-pval)) %>% 
  left_join(pSM_MAM %>% 
              rename(pSM_MAM_coeff = 'soma_deltaMAM') %>% 
              select(markname, pSM_MAM_coeff)) %>% 
  relocate(pM, pM_nlogp, .after = pM_coeff) %>% 
  relocate(pSM_MAM, pSM_MAM_nlogp, .after = pSM_MAM_coeff) %>% 
  relocate(meta_pval_BH, meta_pval_Bonferroni, .after = meta_p)

write.csv(annotated, '002_Filtered_res/CMA_meta_annotated.csv')

#### Check the number of pos & neg hits by either FDR method:
# Bonferonni
nrow(annotated[annotated$direction == 'pos' & annotated$meta_pval_Bonferroni < 0.05, ]) # 89
nrow(annotated[annotated$direction == 'neg' & annotated$meta_pval_Bonferroni < 0.05, ]) # 3

view(annotated[annotated$direction == 'pos' & annotated$meta_pval_Bonferroni < 0.05, ])

# Hochberg
nrow(annotated[annotated$direction == 'pos' & annotated$meta_pval_BH < 0.05, ]) # 222
nrow(annotated[annotated$direction == 'neg' & annotated$meta_pval_BH < 0.05, ]) # 44

### Check if the overlaping positives are in the 29 proteins...
# Load overlap res:
MAM_SAM_elisa_pM_compare <- read.csv("MAM_SAM_pM_compare_0.05fdr.csv") %>% 
  filter(soma_deltaMAM_anova_fdr < 0.05 & pM_fdr < 0.05) %>% 
  select(EntrezGeneSymbol) %>% 
  unique()
sum(annotated[annotated$meta_pval_BH < 0.05,]$EntrezGeneSymbol %in% MAM_SAM_elisa_pM_compare$EntrezGeneSymbol) # 29, all are present

#######################################################################################
### Generate figure:
#######################################################################################
cma_res <- read.csv('002_Filtered_res/CMA_meta_annotated.csv') %>% 
  select(-X) %>% 
  mutate(dir_log10_BHpval = ifelse(direction == 'pos', -log10(meta_pval_BH), log10(meta_pval_BH)))

cma_res_ordered <- cma_res %>% 
  arrange(dir_log10_BHpval) %>% 
  mutate(num = 1:nrow(cma_res))

cma_res_ordered_pos <- cma_res_ordered %>% 
  filter(meta_pval_BH < 0.05 & direction == 'pos')

cma_res_ordered_nonsig <- cma_res_ordered %>% 
  filter(meta_pval_BH >= 0.05)

cma_res_ordered_neg <- cma_res_ordered %>% 
  filter(meta_pval_BH < 0.05 & direction == 'neg')
#######################################################################################
### over-representation analysis results:
#######################################################################################
library(clusterProfiler)
library(wordcloud)
library(tidyverse)
library(msigdbr)
library(enrichplot)

# Load Msigdb hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# Load human genes
organism = 'org.Hs.eg.db'
# BiocManager::install(organism, character.only = TRUE, force = TRUE) # do not update packages when offered during installation!
library(organism, character.only = TRUE)

# Create list of the genes to screen:
cma_res_anno <- read.csv('002_Filtered_res/CMA_meta_01_extended_annotations.csv')
pos_genes <- unique(cma_res_anno[cma_res_anno$meta_pval_BH < 0.05 & cma_res_anno$direction == 'pos',]$EntrezGeneSymbol)

# Load the entrez IDs 
entrezIDs_pos <- read.csv('MAM_lm_res_wlz_soma_area_group_gender.csv') %>% 
  filter(EntrezGeneSymbol %in% pos_genes) %>% 
  dplyr::select(EntrezGeneID) %>% 
  unique()

go_enrich_pos <- enrichGO(gene = entrezIDs_pos$EntrezGeneID,
                      # universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05)

go_enrich_pos.df <- as.data.frame(go_enrich_pos)

go_enrich_pos.df %>% 
  filter(grepl('LEP', geneID)) %>% 
  view()

leading_edges <- sort(str_split(go_enrich_pos.df[go_enrich_pos.df$Description== 'ossification',]$geneID, "/")[[1]])

length(unique(leading_edges)) # 25

write.csv(go_enrich_pos.df, '003_module/GO_results_pos_proteins.csv')

# Create list of the genes to screen:
neg_genes <- unique(cma_res_anno[cma_res_anno$meta_pval_BH < 0.05 & cma_res_anno$direction == 'neg',]$EntrezGeneSymbol)

# Load the entrez IDs 
entrezIDs_neg <- read.csv('MAM_lm_res_wlz_soma_area_group_gender.csv') %>% 
  filter(EntrezGeneSymbol %in% neg_genes) %>% 
  dplyr::select(EntrezGeneID) %>% 
  unique()

go_enrich_neg <- enrichGO(gene = entrezIDs_neg$EntrezGeneID,
                          # universe = names(gene_list),
                          OrgDb = organism, 
                          keyType = 'ENTREZID',
                          readable = T,
                          ont = "BP",
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1)

go_enrich_neg.df <- as.data.frame(go_enrich_neg)


write.csv(go_enrich_neg.df, '003_module/GO_results_neg_proteins.csv')

somaDf <- cma_res_anno %>% 
  left_join(soma_anno_full %>% 
              dplyr::select(EntrezGeneID, EntrezGeneSymbol, TargetFullName) %>% 
              dplyr::arrange(EntrezGeneID) %>% 
              dplyr::filter(!duplicated(EntrezGeneID))) %>% 
  mutate(dir_pval = ifelse(direction == 'pos', -log10(meta_p), log10(meta_p))) %>% 
  dplyr::select(EntrezGeneID, EntrezGeneSymbol, dir_pval) %>% 
  na.omit()

### Study the responses of the individual markers
gsea_res.df <- as.data.frame(gsea_res) %>% 
  separate_rows(core_enrichment, sep = '/') %>% 
  left_join(soma_anno_full %>% 
              rename(core_enrichment = 'EntrezGeneID') %>% 
              dplyr::select(core_enrichment, EntrezGeneSymbol, TargetFullName) %>% 
              dplyr::arrange(core_enrichment) %>% 
              dplyr::filter(!duplicated(core_enrichment))) %>% 
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge) %>% 
  dplyr::summarise(core_enrichment = paste(core_enrichment, collapse = ","),
                   EntrezGeneSymbol = paste(EntrezGeneSymbol, collapse = ","),
                   TargetFullName = paste(TargetFullName, collapse = ",")) %>% 
  arrange(qvalue)

gsea_res.df[gsea_res.df$ID=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]$EntrezGeneSymbol
gsea_res.df[gsea_res.df$ID=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',]$TargetFullName

gsea_res.df[gsea_res.df$ID=='HALLMARK_INTERFERON_GAMMA_RESPONSE',]$EntrezGeneSymbol
gsea_res.df[gsea_res.df$ID=='HALLMARK_INTERFERON_GAMMA_RESPONSE',]$TargetFullName

gsea_res.df[gsea_res.df$ID=='HALLMARK_INTERFERON_ALPHA_RESPONSE',]$EntrezGeneSymbol
gsea_res.df[gsea_res.df$ID=='HALLMARK_INTERFERON_ALPHA_RESPONSE',]$TargetFullName


detach("package:org.Hs.eg.db", unload=TRUE)
detach("package:clusterProfiler", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)

#######################################################################################
### Recreate figure colored by the ORA results:
#######################################################################################

go_enrich_pos.df <- read.csv('003_module/GO_results_pos_proteins.csv')
go_enrich_neg.df <- read.csv('003_module/GO_results_neg_proteins.csv')

cellularStress <- c('GDF15', 'HSPA1A', 'HSP90AB1', 'HSP90AA1', 'HSPH1')
immuneResponse <- c('C9', 'C2', 'C4A|C4B', 'C6', 'GNLY')

axon_dev = cma_res$EntrezGeneSymbol[cma_res$EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='axon development',]$geneID, split = '/')[[1]]]
axon_dev = cma_res$EntrezGeneSymbol[cma_res$EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='axon development',]$geneID, split = '/')[[1]]]

cma_res %>% 
  mutate(axon_dev = ifelse(EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='axon development',]$geneID, split = '/')[[1]] & as.numeric(meta_pval_BH) < 0.05, 1, 0),
         ECM = ifelse(EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='extracellular matrix organization',]$geneID, split = '/')[[1]] & as.numeric(meta_pval_BH) < 0.05, 1, 0),
         nutrient = ifelse(EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='response to nutrient levels',]$geneID, split = '/')[[1]]& as.numeric(meta_pval_BH) < 0.05, 1, 0),
         ossification = ifelse(EntrezGeneSymbol %in% strsplit(go_enrich_pos.df[go_enrich_pos.df$Description=='ossification',]$geneID, split = '/')[[1]]& as.numeric(meta_pval_BH) < 0.05, 1, 0),
         cellularStress = ifelse(EntrezGeneSymbol %in% cellularStress & as.numeric(meta_pval_BH) < 0.05, 1, 0),
         immuneResponse = ifelse(EntrezGeneSymbol %in% immuneResponse & as.numeric(meta_pval_BH) < 0.05, 1, 0)) %>% 
  mutate(BP = ifelse(axon_dev==1, 'axon_dev',
                     ifelse(ECM==1, 'ECM',
                            ifelse(nutrient==1, 'nutrient',
                                   ifelse(ossification==1,'ossification',
                                          ifelse(cellularStress==1,'cellularStress',
                                                 ifelse(immuneResponse==1, 'immuneResponse', 'other'))))))) %>% 
  # filter(BP!='other') %>% view()
  arrange(dir_log10_BHpval) %>% 
  mutate(num = 1:nrow(cma_res)) %>% 
  ggplot(aes(x = num, y = dir_log10_BHpval, color = BP)) + 
  geom_point(shape = 18, size = 2) + 
  # geom_segment(aes(x=num, xend = num, y = 0, yend = dir_log10_BHpval, color = color), alpha = 0.05) + 
  scale_color_manual(values = c('blue3', 'yellow3', 'purple3', 'red3', 'green3', 'orange3', 'grey')) +
  coord_flip() + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(-1.3, 1.3), linetype = 'dashed') + 
  coord_flip() + 
  xlab('Ranked plasma protein') + 
  ylab('log10(q-value)') + 
  theme_classic() + 
  scale_y_continuous(breaks = c(-6,-1.3, 0, 1.3, 3, 6, 9))+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  ggtitle('q-value distribution of covarying proteins')

ggsave('/Users/stevenhartman/Library/CloudStorage/Box-Box/000Gordon/Projects/Misc/CMA/SOMA_CMAaptamerID1_filter2_posNeg_pSMMAM_pM/002_Filtered_res/qval_dist_proteinsSameDirAcrossStudies_002.pdf', width = 5, height = 3)

#####################################################################################
# Create the figure table (top 30)
#####################################################################################
# Load annotated table, add in the CMA p-value, and add the presence or absence
# in biological processes

# Annotated table:
axonDev = strsplit(go_enrich_pos.df[go_enrich_pos.df$Description == 'axon development',]$geneID, split = '/')

BP_genes <- strsplit(go_enrich_pos.df$geneID, split = '/')
names(BP_genes) <- go_enrich_pos.df$Description


pub_soma_anno <- read.csv('003_module/publication_annotations.csv') %>% 
  separate(col = protein, into = c('protein', 'EntrezGeneID'), sep = ' \\(') %>% 
  mutate(EntrezGeneID = gsub('\\)', '', EntrezGeneID)) %>% 
  mutate(axonDev = ifelse(EntrezGeneID %in% BP_genes$`axon development`, 1, 0),
         ECM = ifelse(EntrezGeneID %in% BP_genes$`extracellular matrix organization`, 1, 0),
         nutrient = ifelse(EntrezGeneID %in% BP_genes$`response to nutrient levels`, 1, 0),
         ossification = ifelse(EntrezGeneID %in% BP_genes$ossification, 1, 0)) %>% 
  left_join(cma_res %>% 
              dplyr::rename(EntrezGeneID = 'markname') %>% 
              select(EntrezGeneID, dir_log10_BHpval))

cma_res_top_hits[cma_res_top_hits$EntrezGeneID %in% c('CNTN5', 'CRABP2', 'ROBO2', 'RET'), 'manualTheme'] <- 'CNS/Neuronal function'
cma_res_top_hits[cma_res_top_hits$EntrezGeneID %in% c('COMP', 'COL9A1', 'CD248'), 'manualTheme'] <- 'Connective tissue/musculoskeletal and organ development'
cma_res_top_hits[cma_res_top_hits$EntrezGeneID %in% c('DLK1'), 'manualTheme'] <- 'Other'

cma_res_top_hits <- cma_res_top_hits %>% 
  mutate(axonDev = ifelse(EntrezGeneID %in% BP_genes$`axon development`, 1, 0),
         ECM = ifelse(EntrezGeneID %in% BP_genes$`extracellular matrix organization`, 1, 0),
         nutrient = ifelse(EntrezGeneID %in% BP_genes$`response to nutrient levels`, 1, 0),
         ossification = ifelse(EntrezGeneID %in% BP_genes$ossification, 1, 0)) %>% 
  mutate(final_annotation = ifelse(!is.na(pSM_man_function), pSM_man_function, soma100_panel_man_function)) %>% 
  arrange(manualTheme, desc(dir_log10_BHpval))

cma_res_top_hits %>% 
  select(dir_log10_BHpval, EntrezGeneID, axonDev, ECM, nutrient, ossification) %>% 
  pivot_longer(cols = axonDev:ossification,
               names_to = 'BP', 
               values_to = 'presence') %>% 
  ggplot(aes(x=BP, y=factor(EntrezGeneID, levels = rev(cma_res_top_hits$EntrezGeneID)), fill = factor(presence))) + 
  geom_tile(color = 'grey') + 
  ylab('') + 
  scale_fill_manual(values=c('white', 'black'))

ggsave('/003_module/BP_presence_001.pdf',
       height = 7,
       width = 3.3)
 
write.csv(cma_res_top_hits, '003_module/cma_res_top30_hits.csv')

############################################################################
# Create supplemental table for CMA (Add to the WLZ association table)
############################################################################
# Load in the WLZ associations
MAM_SAM_pM_compare <- read.csv("MAM_SAM_pM_compare_0.05fdr.csv") %>% 
  select(-X, -elisa_deltaMAM, -elisa_deltaMAM_fdr, -elisa_deltaSAM, -elisa_deltaSAM_fdr, -pSM_elisa_MAM, -pSM_elisa_SAM, -pSM_SOMA_MAM, -pSM_SOMA_SAM, -pM_SOMA,
         -soma_deltaSAM, -soma_deltaSAM_sde, -soma_deltaSAM_pval_anova, -soma_deltaSAM_anova_fdr)

CMA_output <- read.csv('001_CMA_res/CMA_meta.csv')

# Write the results for a supplemental table:
supp_table <- CMA_output %>% 
  left_join(soma_anno %>% 
              dplyr::rename(markname = 'SeqId') %>% 
              dplyr::select(TargetFullName, EntrezGeneSymbol, markname, SeqId.version) %>% 
              arrange(markname) %>% 
              filter(duplicated(markname) == F) %>% 
              unique())  %>% 
  left_join(pM %>% 
              dplyr::rename(pM_coeff = 'pearsonRho') %>% 
              select(-pval, -X, -EntrezGeneSymbol)) %>% 
  left_join(pSM_MAM %>% 
              dplyr::rename(pSM_MAM_coeff = 'soma_deltaMAM') %>% 
              select(markname, pSM_MAM_coeff)) %>% 
  dplyr::rename(pSM_MAM_pval = 'pSM_MAM') %>% 
  dplyr::rename(pM_pval = 'pM') %>% # mutate(meta_pval_BH = p.adjust(meta_p, method = 'BH')) %>% filter(sign(pSM_MAM_coeff) == sign(pM_coeff) & meta_pval_BH < 0.05) %>% dim() # 262!? 
  mutate(meta_pval_BH = p.adjust(meta_p, method = 'BH'),
         direction_MAMphase = ifelse(pSM_MAM_coeff >= 0, 'pos', 'neg'),
         direction_pM = ifelse(pM_coeff >= 0, 'pos', 'neg')) %>% 
  mutate(direction_MAMphase_pM = paste0(direction_MAMphase,"/",direction_pM)) %>% 
  select(markname, SeqId.version, EntrezGeneSymbol, TargetFullName, 
         pSM_MAM_coeff, pSM_MAM_pval, 
         pM_coeff, pM_pval,
         meta_p,direction_MAMphase_pM, meta_pval_BH) %>% 
  mutate(meta_pval_BH = ifelse(sign(pM_coeff) == sign(pSM_MAM_coeff), meta_pval_BH, '-'))# %>% 
  mutate(meta_pval_BH = ifelse(pM_coeff > 0 & pSM_MAM_coeff > 0, meta_pval_BH, '-'))# %>% 
  mutate(meta_pval_BH = ifelse(direction_MAMphase_pM == 'neg/pos', '-', meta_pval_BH))
  
sum(is.na(supp_table$pM_coeff))
sum(supp_table$meta_pval_BH < 0.05 & supp_table$meta_pval_BH != '-')

write.csv(supp_table, '003_module/cma_res_supplemental.csv')

res <- supp_table %>% 
  filter(meta_pval_BH!='-') %>% 
  mutate(meta_pval_BH = as.numeric(meta_pval_BH))

sum(res$meta_pval_BH < 0.05)

############################################################################
# Identify proteins associated with WLZ across MAM-phases (CMA) that overlap with the SAM phase linear models
# Can run as standalone
############################################################################
library(tidyverse)

wlz_cma <- read.csv('003_module/cma_res_supplemental.csv') %>% 
  drop_na()

sig_cma <- wlz_cma %>% 
  filter(meta_pval_BH != '-') %>% 
  mutate(meta_pval_BH = as.numeric(meta_pval_BH)) %>% # view()
  filter(meta_pval_BH < 0.05)
 
sum(sig_cma$direction_MAMphase_pM == 'pos/pos') # 222 <- this is the unique aptamers
sum(sig_cma$direction_MAMphase_pM == 'neg/neg') # 44

sig_cma %>% 
  select(EntrezGeneSymbol, direction_MAMphase_pM) %>% 
  unique() %>% 
  filter(direction_MAMphase_pM == 'pos/pos') %>% 
  dim() # 215 unique proteins are positive

wlz_sam <- read.csv('pSM_MAM_SAM_compare.csv')

sig_sam <- wlz_sam %>% 
  filter(soma_deltaSAM_anova_fdr < 0.05)

# Find positive overlap
pos_overlap_aptamers <- intersect(sig_sam[sig_sam$soma_deltaSAM > 0,]$SeqId, sig_cma[sig_cma$direction_MAMphase_pM == 'pos/pos',]$markname) # 28 aptamers

# Find negatifve overlap
neg_overlap_aptamers <- intersect(sig_sam[sig_sam$soma_deltaSAM > 0,]$SeqId, sig_cma[sig_cma$direction_MAMphase_pM == 'neg/neg',]$markname) # empty

# Create sigfig function
sigfigs <- function(vec, num_figs = 3) {
  formatC(signif(vec, digits=num_figs), digits=num_figs, format="fg", flag="#")
}

# Create table of the overlapping aptamers:
overlap_table <- wlz_sam %>% 
  filter(SeqId %in% pos_overlap_aptamers) %>% 
  select(SeqId, SeqIdVersion, TargetFullName, EntrezGeneID, EntrezGeneSymbol, soma_deltaSAM, soma_deltaSAM_anova_fdr) %>% 
  left_join(sig_cma %>% 
              rename(SeqId = 'markname') %>% 
              select(SeqId, meta_pval_BH, direction_MAMphase_pM)) %>% 
  mutate(soma_deltaSAM = sigfigs(soma_deltaSAM),
         soma_deltaSAM_anova_fdr= sigfigs(soma_deltaSAM_anova_fdr),
         meta_pval_BH = sigfigs(meta_pval_BH))

# Write table to file
write.csv(overlap_table, '003_module/CMA_overlap_SAMphase_WLZassociations.csv')

















