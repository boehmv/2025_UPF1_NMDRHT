#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Transcript_Level
# Objective: Perform transcript-level data preparation and initial analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)
library(tximeta)

##
# Load GENCODE-based Salmon QC data -----------------------------------------------------
##

# Use TidyMultiqc for importing multiple multiqc output files
GENCODE_UPF1_Salmon_QC_forComparison <- TidyMultiqc::load_multiqc(UPF1_NMDRHT_datasets %>% 
                                                                    filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                                                                                "HCT116_UPF1_AID_recovery_this_Study",
                                                                                                "HCT116_UPF1_AID_degradation_riboMinus_this_Study",
                                                                                                "HCT116_UPF1_FKBP_degradation_this_Study",
                                                                                                "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
                                                                    mutate(QC_path = paste0(path,"/QC/QC_salmon/multiqc_data/multiqc_data.json")) %>% 
                                                                    pull(QC_path)) %>% 
  left_join(UPF1_NMDRHT_datasets_experiments,
            by = c("metadata.sample_id" = "sample")) %>% 
  filter(experimentSet != "HEK293_NMD_Boehm_2018")

# Save as csv
GENCODE_UPF1_Salmon_QC_forComparison %>% write_csv("Resources/QC/GENCODE_UPF1_Salmon_QC_forComparison.csv")

##
# Load NMDRHT-based Salmon QC data -----------------------------------------------------
##

# Use TidyMultiqc for importing multiple multiqc output files
NMDRHT_UPF1_Salmon_QC_forComparison <- TidyMultiqc::load_multiqc(UPF1_NMDRHT_datasets %>% 
                                                                   filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                                                                               "HCT116_UPF1_AID_recovery_this_Study",
                                                                                               "HCT116_UPF1_AID_degradation_riboMinus_this_Study",
                                                                                               "HCT116_UPF1_FKBP_degradation_this_Study",
                                                                                               "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
                                                                   mutate(QC_path = paste0(path,"/QC/QC_salmon_NMD/multiqc_data/multiqc_data.json")) %>% 
                                                                   pull(QC_path)) %>% 
  left_join(UPF1_NMDRHT_datasets_experiments,
            by = c("metadata.sample_id" = "sample")) %>% 
  filter(experimentSet != "HEK293_NMD_Boehm_2018")

# Save as csv
NMDRHT_UPF1_Salmon_QC_forComparison %>% write_csv("Resources/QC/NMDRHT_UPF1_Salmon_QC_forComparison.csv")



###
# Cluster Transcript -----------------------------------------------------------
###

# Get List of significantly upregulated transcripts from edgeR DTE - NMDRHT analysis
# this means: retain any transcript_id which is at least once (in >= 1 condition) significantly upregulated
Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  group_by(transcript_id) %>% 
  filter(FDR < 0.0001) %>% 
  filter(logFC > 1) %>% 
  distinct(transcript_id) %>% 
  pull(transcript_id)

# Downregulated as well
Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  group_by(transcript_id) %>% 
  filter(FDR < 0.0001) %>% 
  filter(logFC < -1) %>% 
  distinct(transcript_id) %>% 
  pull(transcript_id)

### Clustering transcript-level ---------------------------------------------------

#
# Note: first check how the DESeq2-"based" vsd method for normalization is correlated with the edgeR-based "logcpm" method
#

# Start with DESeq2-approach - Note: no count scaling based on overdispersion was performed which might explain slight differences later on
Rev_1_F5_mydir1="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"
Rev_1_F5_mydir2="/home/volker/gencode.v42.datasets/2023_UPF1_AID_recovery_DW"

# Get samples and check if all files are present
Rev_1_F5_samples1 <- read.table(file.path(Rev_1_F5_mydir1, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) 

Rev_1_F5_samples2 <- read.table(file.path(Rev_1_F5_mydir2, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "Recovery_control",
                               TRUE ~ condition))

Rev_1_F5_samples <- as_tibble(bind_rows(Rev_1_F5_samples1,
                                        Rev_1_F5_samples2))

# Get unique conditions
Rev_1_F5_condition <- Rev_1_F5_samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate data frame used for tximeta
Rev_1_F5_coldata <- tibble(files = file.path(Rev_1_F5_mydir1, "Salmon_NMD", Rev_1_F5_samples1$sample, "quant.sf")) %>% 
  bind_rows(tibble(files = file.path(Rev_1_F5_mydir2, "Salmon_NMD", Rev_1_F5_samples2$sample, "quant.sf"))) 

# Supplement with sample IDs as "names"
Rev_1_F5_coldata$names <- Rev_1_F5_samples$sample

# Join with samples to get condition
Rev_1_F5_coldata <- Rev_1_F5_coldata %>% 
  left_join(Rev_1_F5_samples,
            by = c("names" = "sample")) %>% 
  mutate(condition = fct_relevel(as_factor(condition),
                                 "control"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(Rev_1_F5_coldata$files)))

# Use tximport to import salmon data
Rev_1_F5_all_y <- tximeta(Rev_1_F5_coldata,
                          type = "salmon",
                          txOut = TRUE,
                          useHub = FALSE) # reads in counts and inf reps

# Check which assays are available
assayNames(Rev_1_F5_all_y)

# Label spike-ins (ERCC and SIRV) as FALSE in keep column -> do not analyze them!
mcols(Rev_1_F5_all_y)$keep <- !str_detect(rownames(assays(Rev_1_F5_all_y)[["counts"]]), 'ERCC') & !str_detect(rownames(assays(Rev_1_F5_all_y)[["counts"]]), 'SIRV')

# Remove spike-in genes from analysis (not just labelled as "keep == FALSE", but really remove)
Rev_1_F5_all_y <- Rev_1_F5_all_y[mcols(Rev_1_F5_all_y)$keep,]

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSet(Rev_1_F5_all_y,
                       design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Extracting count data log2-transformed values 
# According to https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog
dds_vsd <- vst(ddsTxi, blind=FALSE)

# Which dimensions does the matrix have?
dim(dds_vsd)

### Normalize counts --------------------------------------------------------

# extract the normalized, but not z-scaled counts of NMD:
Rev_1_F5_transcript_cts_up <- assay(dds_vsd)[rownames(assay(dds_vsd)) %in% Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list,]
Rev_1_F5_transcript_cts_down <- assay(dds_vsd)[rownames(assay(dds_vsd)) %in% Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list,]

# Z-scale the log2-transformed count matrix
Rev_1_F5_transcript_test_up <- (t(scale(t(Rev_1_F5_transcript_cts_up))))
Rev_1_F5_transcript_test_down <- (t(scale(t(Rev_1_F5_transcript_cts_down))))

##
#### edgeR variant -----------------------------------------------------------
##

# Generate data frame used for tximeta
Rev_1_F5_coldata_edgeR <- tibble(files = file.path(Rev_1_F5_mydir1, "Salmon_NMD", Rev_1_F5_samples1$sample)) %>% 
  bind_rows(tibble(files = file.path(Rev_1_F5_mydir2, "Salmon_NMD", Rev_1_F5_samples2$sample))) 

# Supplement with sample IDs as "names"
Rev_1_F5_coldata_edgeR$names <- Rev_1_F5_samples$sample

# Join with samples to get condition
Rev_1_F5_coldata_edgeR <- Rev_1_F5_coldata_edgeR %>% 
  left_join(Rev_1_F5_samples,
            by = c("names" = "sample")) %>% 
  mutate(condition = fct_relevel(as_factor(condition),
                                 "control"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(Rev_1_F5_coldata_edgeR$files)))

# Import transcript-level counts 
catch <- edgeR::catchSalmon(paths = Rev_1_F5_coldata_edgeR$files)

# Account for the mapping ambiguity
scaled.counts <- catch$counts/catch$annotation$Overdispersion

# Create DGEList object
DGEList <- edgeR::DGEList(counts = scaled.counts,
                   samples = Rev_1_F5_samples,
                   group = Rev_1_F5_samples$condition,
                   genes = catch$annotation)

# Filter out spike-ins
keep_woSpike <- !str_detect(rownames(DGEList$genes), 'ERCC') & !str_detect(rownames(DGEList$genes), 'SIRV')
names(keep_woSpike) <- rownames(DGEList$genes)

print("Kept transcripts")
table(keep_woSpike)

DGEList_filt <- DGEList[keep_woSpike, , keep.lib.sizes=FALSE]

# Perform normalization
DGEList_filt_norm <- edgeR::normLibSizes(DGEList_filt)

# see e.g. https://support.bioconductor.org/p/130683/
logcpm <- edgeR::cpm(DGEList_filt_norm, log=TRUE)

# save data to allow reproducing plots
save(Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list,
     Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list,
     dds_vsd,
     Rev_1_F5_transcript_test_up,
     logcpm,
     file = paste0("Resources/Cluster/Rev_1_F5_Raw_forClustering.rds"))

# load if necessary
# load("Resources/Cluster/Rev_1_F5_Raw_forClustering.rds")

# scatter plot of two samples - from both edgeR and DESeq2 (vsd)
bind_rows(as_tibble(logcpm[, 1:2]) %>%
            mutate(transformation = "edgeR") %>% 
            magrittr::set_colnames(c("x", "y", "transformation")),
          as_tibble(assay(dds_vsd)[rownames(assay(dds_vsd)) %in% c(Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list, Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list),1:2]) %>%
            mutate(transformation = "DESeq2") %>% 
            magrittr::set_colnames(c("x", "y", "transformation"))) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_hex(bins = 80) +
  coord_fixed() + 
  facet_grid(~transformation)  

# Z-scale the log2-transformed count matrix
# edgeR
Rev_1_F5_transcript_test_edgeR <- (t(scale(t(logcpm))))

# rename columns to make consistent with DESeq2 approach
colnames(Rev_1_F5_transcript_test_edgeR) <- colnames(Rev_1_F5_transcript_test_up)

# Z-scale and select only significant up- and downregulated transcripts
Rev_1_F5_transcript_test_DESeq2 <- (t(scale(t(assay(dds_vsd)[rownames(assay(dds_vsd)) %in% c(Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list, Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list),]))))

# Final verdict: both DESeq2 and edgeR-based methods are highly correlated 
# due to edgeR dealing with read-to-transcript ambiguity -> we will select the edgeR-based method

### Normalize counts --------------------------------------------------------

# extract the normalized, but not z-scaled counts of NMD:
Rev_1_F5_transcript_cts_up <- logcpm[rownames(logcpm) %in% Rev_1_F5_edgeR_DTE_NMDRHT_sig_upregulated_list,]
Rev_1_F5_transcript_cts_down <- logcpm[rownames(logcpm) %in% Rev_1_F5_edgeR_DTE_NMDRHT_sig_downregulated_list,]

# Z-scale the log2-transformed count matrix
Rev_1_F5_transcript_test_up <- (t(scale(t(Rev_1_F5_transcript_cts_up))))
Rev_1_F5_transcript_test_down <- (t(scale(t(Rev_1_F5_transcript_cts_down))))

###  upregulated - Hierarchical clustering --------------------------------------------------------
set.seed(321)
Rev_1_F5_transcript_dist_up <- dist(Rev_1_F5_transcript_test_up)

Rev_1_F5_transcript_hclust_up <- hclust(Rev_1_F5_transcript_dist_up, method = "ward.D2")

# # Elbow method
# set.seed(321)
# Rev_1_F5_transcript_hclust_wss_plot <- fviz_nbclust(Rev_1_F5_transcript_test_up, hcut, method = "wss") +
#   labs(subtitle = "Hierarchical - Elbow method")
# 
# # # Silhouette method
# set.seed(321)
# Rev_1_F5_transcript_hclust_silhouette_plot <- fviz_nbclust(Rev_1_F5_transcript_test_up, hcut, method = "silhouette") +
#   labs(subtitle = "Hierarchical - Silhouette method")
# 
# # # Gap statistic
# set.seed(321)
# Rev_1_F5_transcript_hclust_gap_plot <- fviz_nbclust(Rev_1_F5_transcript_test_up, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
#   labs(subtitle = "Hierarchical - Gap statistic method")

Rev_1_F5_transcript_hclust_up_transcript_cluster <- cutree(Rev_1_F5_transcript_hclust_up, k = 4) %>% 
  enframe() %>% 
  dplyr::rename(transcript_id = name, cluster = value)

head(Rev_1_F5_transcript_hclust_up_transcript_cluster)

table(cutree(Rev_1_F5_transcript_hclust_up, k = 4))

Rev_1_F5_transcript_hclust_up_transcript_cluster <- Rev_1_F5_transcript_hclust_up_transcript_cluster %>% 
  mutate(label=case_when(cluster == 1 ~ "1:early",
                         cluster == 2 ~ "2:delayed",
                         cluster == 3 ~ "3:late",
                         cluster == 4 ~ "4:inverse")) %>% 
  pull(label, name=transcript_id)

###  downregulated - Hierarchical clustering --------------------------------------------------------
Rev_1_F5_transcript_dist_down <- dist(Rev_1_F5_transcript_test_down)

Rev_1_F5_transcript_hclust_down <- hclust(Rev_1_F5_transcript_dist_down, method = "ward.D2")

# # Elbow method
# set.seed(321)
# Rev_1_F5_transcript_down_hclust_wss_plot <- fviz_nbclust(Rev_1_F5_transcript_test_down, hcut, method = "wss") +
#   labs(subtitle = "Hierarchical - Elbow method")
# 
# # Silhouette method
# set.seed(321)
# Rev_1_F5_transcript_down_hclust_silhouette_plot <- fviz_nbclust(Rev_1_F5_transcript_test_down, hcut, method = "silhouette") +
#   labs(subtitle = "Hierarchical - Silhouette method")
# 
# # Gap statistic
# set.seed(321)
# Rev_1_F5_transcript_down_hclust_gap_plot <- fviz_nbclust(Rev_1_F5_transcript_test_down, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
#   labs(subtitle = "Hierarchical - Gap statistic method")

Rev_1_F5_transcript_hclust_down_transcript_cluster <- cutree(Rev_1_F5_transcript_hclust_down, k = 4) %>% 
  enframe() %>% 
  dplyr::rename(transcript_id = name, cluster = value)

head(Rev_1_F5_transcript_hclust_down_transcript_cluster)

table(cutree(Rev_1_F5_transcript_hclust_down, k = 4))

Rev_1_F5_transcript_hclust_down_transcript_cluster <- Rev_1_F5_transcript_hclust_down_transcript_cluster %>% 
  mutate(label=case_when(cluster == 1 ~ "1:early",
                         cluster == 2 ~ "3:late",
                         cluster == 3 ~ "4:inverse",
                         cluster == 4 ~ "2:delayed"
  )) %>% 
  pull(label, name=transcript_id)

#### Generate data frames -------------------------------------------------------------------------

# Get transcript_id to DGE_cluster combination for up- and downregulated
Rev_1_F5_transcript_klus_cluster_up_df <- as_tibble(Rev_1_F5_transcript_hclust_up_transcript_cluster, rownames = "transcript_id") %>% 
  dplyr::rename("DTE_cluster_up" = "value") %>% 
  mutate(DTE_cluster_up = paste0("up ",DTE_cluster_up))

Rev_1_F5_transcript_klus_cluster_down_df <- as_tibble(Rev_1_F5_transcript_hclust_down_transcript_cluster, rownames = "transcript_id") %>% 
  dplyr::rename("DTE_cluster_down" = "value") %>% 
  mutate(DTE_cluster_down = paste0("down ",DTE_cluster_down))

### Merge with NMDRHT ------------------------------------------------------

# List of expressed transcripts
Rev_1_F5_edgeR_DTE_NMDRHT_expressed_transcript_ids <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  dplyr::select(transcript_id) %>% 
  dplyr::distinct(transcript_id)

# Merge NMDRHT annotation with cluster from edgeR DTE analysis
NMDRHT.v1.2_tbl_DTE_cluster <- NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  left_join(Rev_1_F5_transcript_klus_cluster_up_df,
            by = c("transcript_id" = "transcript_id")) %>% 
  left_join(Rev_1_F5_transcript_klus_cluster_down_df,
            by = c("transcript_id" = "transcript_id")) %>% 
  mutate(DTE_cluster = case_when(is.na(DTE_cluster_up) & is.na(DTE_cluster_down) ~ "NA",
                                 !is.na(DTE_cluster_up) & is.na(DTE_cluster_down) ~ DTE_cluster_up,
                                 is.na(DTE_cluster_up) & !is.na(DTE_cluster_down) ~ DTE_cluster_down,
                                 !is.na(DTE_cluster_up) & !is.na(DTE_cluster_down) ~ paste0("complex"))) %>% 
  mutate(DTE_cluster = case_when(DTE_cluster == "NA" & transcript_id %in% Rev_1_F5_edgeR_DTE_NMDRHT_expressed_transcript_ids$transcript_id ~ "expressed",
                                 DTE_cluster == "NA" & !transcript_id %in% Rev_1_F5_edgeR_DTE_NMDRHT_expressed_transcript_ids$transcript_id ~ "not_expressed",
                                 TRUE ~ DTE_cluster)) %>% 
  mutate(DTE_cluster = fct_relevel(DTE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "up 4:inverse",
                                   "expressed",
                                   "not_expressed",
                                   "complex",
                                   "down 4:inverse",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early"))

### Fix "complex" transcripts (with up- and down-clusters) ------------------------------------------------------
# Slice by FDR (e.g. most significant transcript expression changes)

edgeR_DTE_NMDRHT_combined_transcript_fix <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  left_join(NMDRHT.v1.2_tbl_DTE_cluster) %>% 
  filter(DTE_cluster == "complex") %>%
  group_by(transcript_id) %>% 
  mutate(min_logFC = min(logFC),
         max_logFC = max(logFC),
         min_FDR = min(FDR),
         median(logFC)) %>% 
  slice_min(FDR, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(DTE_cluster_fix= case_when(logFC > 0 ~ DTE_cluster_up,
                                    logFC < 0 ~ DTE_cluster_down),
         .after=DTE_cluster)

# Fix transcripts in relevant dataframe
NMDRHT.v1.2_tbl_DTE_cluster <- NMDRHT.v1.2_tbl_DTE_cluster %>% 
  left_join(edgeR_DTE_NMDRHT_combined_transcript_fix %>% dplyr::select(transcript_id, transcript_name, DTE_cluster, DTE_cluster_fix)) %>% 
  mutate(DTE_cluster = case_when(DTE_cluster == "complex" ~ DTE_cluster_fix,
                                 TRUE ~ DTE_cluster)) %>% 
  dplyr::select(-DTE_cluster_fix) %>% 
  mutate(DTE_cluster = fct_relevel(DTE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "up 4:inverse",
                                   "expressed",
                                   "not_expressed",
                                   "down 4:inverse",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early"))

# Check DTE_cluster count
NMDRHT.v1.2_tbl_DTE_cluster %>% 
  dplyr::count(DTE_cluster)

# Fix complex genes in cluster data
## Up
Rev_1_F5_transcript_klus_cluster_up_df_fix <- Rev_1_F5_transcript_klus_cluster_up_df %>% 
  filter(transcript_id %in% (NMDRHT.v1.2_tbl_DTE_cluster %>% 
                               filter(DTE_cluster %in% c("up 1:early",
                                                         "up 2:delayed",
                                                         "up 3:late",
                                                         "up 4:inverse")) %>% 
                               pull(transcript_id)))

# Counts
Rev_1_F5_transcript_klus_cluster_up_df_fix %>% 
  dplyr::count()

Rev_1_F5_transcript_klus_cluster_up_df_fix %>% 
  dplyr::count(DTE_cluster_up)

### in cluster vector
Rev_1_F5_transcript_hclust_up_transcript_cluster_fix <- Rev_1_F5_transcript_klus_cluster_up_df_fix %>% 
  pull(DTE_cluster_up, name=transcript_id)

## Down
Rev_1_F5_transcript_klus_cluster_down_df_fix <- Rev_1_F5_transcript_klus_cluster_down_df %>% 
  filter(transcript_id %in% (NMDRHT.v1.2_tbl_DTE_cluster %>% 
                               filter(DTE_cluster %in% c("down 4:inverse",
                                                         "down 3:late",
                                                         "down 2:delayed",
                                                         "down 1:early")) %>% 
                               pull(transcript_id)))

# Counts
Rev_1_F5_transcript_klus_cluster_down_df_fix %>% 
  dplyr::count()

Rev_1_F5_transcript_klus_cluster_down_df_fix %>% 
  dplyr::count(DTE_cluster_down)

### in cluster vector
Rev_1_F5_transcript_hclust_down_transcript_cluster_fix <- Rev_1_F5_transcript_klus_cluster_down_df_fix %>% 
  pull(DTE_cluster_down, name=transcript_id)

### Prepare for Heatmap
Rev_1_F5_transcript_test_up_fix <- subset(Rev_1_F5_transcript_test_up, rownames(Rev_1_F5_transcript_test_up) %in% Rev_1_F5_transcript_klus_cluster_up_df_fix$transcript_id)

Rev_1_F5_transcript_test_down_fix <- subset(Rev_1_F5_transcript_test_down, rownames(Rev_1_F5_transcript_test_down) %in% Rev_1_F5_transcript_klus_cluster_down_df_fix$transcript_id)

### Fit the log2-norm counts ------------------------------------------------

#### Up ------------------------------------------------

# Z-scale the log2-transformed count matrix
Rev_1_F5_transcript_test_up_klus_cond <- as_tibble(Rev_1_F5_transcript_test_up_fix,
                                                   rownames = "transcript_id") %>% 
  magrittr::set_colnames(c("transcript_id", Rev_1_F5_samples$sample)) %>% 
  janitor::clean_names() %>% 
  left_join(Rev_1_F5_transcript_klus_cluster_up_df_fix) %>% 
  dplyr::rowwise() %>% 
  mutate(control_0h = mean(c(x191495, x191498, x191500)),
         control_48h = mean(c(x191502, x191504, x191506)),
         UPF1_Nter_0h = mean(c(x191508, x191510, x191512)),
         UPF1_Nter_2h = mean(c(x191514, x191516, x191518)),
         UPF1_Nter_4h = mean(c(x191520, x191522, x191524)),
         UPF1_Nter_8h = mean(c(x191526, x191528, x191530)),
         UPF1_Nter_12h = mean(c(x191532, x191534, x191536)),
         UPF1_Nter_24h = mean(c(x191538, x191540, x191542)),
         UPF1_Nter_48h = mean(c(x191544, x191546, x191548)),
         UPF1_Nter_48h = mean(c(x191544, x191546, x191548)),
         control_R0h = mean(c(x200688, x200691, x200693)),
         UPF1_Nter_24h_R0h = mean(c(x200695, x200697, x200699)),
         UPF1_Nter_24h_R2h = mean(c(x200701, x200703, x200705)),
         UPF1_Nter_24h_R4h = mean(c(x200707, x200709, x200711)),
         UPF1_Nter_24h_R8h = mean(c(x200713, x200715, x200717)),
         UPF1_Nter_24h_R12h = mean(c(x200719, x200721, x200723)),
         UPF1_Nter_24h_R24h = mean(c(x200725, x200727, x200729)),
         UPF1_Nter_24h_R48h = mean(c(x200731, x200733, x200735))) %>% 
  dplyr::select(transcript_id,
                control_0h,
                control_48h,
                UPF1_Nter_0h,
                UPF1_Nter_2h,
                UPF1_Nter_4h,
                UPF1_Nter_8h,
                UPF1_Nter_12h,
                UPF1_Nter_24h,
                UPF1_Nter_48h,
                control_R0h,
                UPF1_Nter_24h_R0h,
                UPF1_Nter_24h_R2h,
                UPF1_Nter_24h_R4h,
                UPF1_Nter_24h_R8h,
                UPF1_Nter_24h_R12h,
                UPF1_Nter_24h_R24h,
                UPF1_Nter_24h_R48h,
                DTE_cluster_up) %>% 
  ungroup() %>% 
  mutate(DTE_cluster_up = (fct_relevel(DTE_cluster_up,
                                       "up 1:early",
                                       "up 2:delayed",
                                       "up 3:late",
                                       "up 4:inverse")))

#### Combined Up/Down ------------------------------------------------

# Z-scale the log2-transformed count matrix
Rev_1_F5_transcript_test_down_klus_cond <- as_tibble(Rev_1_F5_transcript_test_down_fix,
                                                     rownames = "transcript_id") %>% 
  magrittr::set_colnames(c("transcript_id", Rev_1_F5_samples$sample)) %>% 
  janitor::clean_names() %>% 
  left_join(Rev_1_F5_transcript_klus_cluster_down_df_fix) %>% 
  dplyr::rowwise() %>% 
  mutate(control_0h = mean(c(x191495, x191498, x191500)),
         control_48h = mean(c(x191502, x191504, x191506)),
         UPF1_Nter_0h = mean(c(x191508, x191510, x191512)),
         UPF1_Nter_2h = mean(c(x191514, x191516, x191518)),
         UPF1_Nter_4h = mean(c(x191520, x191522, x191524)),
         UPF1_Nter_8h = mean(c(x191526, x191528, x191530)),
         UPF1_Nter_12h = mean(c(x191532, x191534, x191536)),
         UPF1_Nter_24h = mean(c(x191538, x191540, x191542)),
         UPF1_Nter_48h = mean(c(x191544, x191546, x191548)),
         UPF1_Nter_48h = mean(c(x191544, x191546, x191548)),
         control_R0h = mean(c(x200688, x200691, x200693)),
         UPF1_Nter_24h_R0h = mean(c(x200695, x200697, x200699)),
         UPF1_Nter_24h_R2h = mean(c(x200701, x200703, x200705)),
         UPF1_Nter_24h_R4h = mean(c(x200707, x200709, x200711)),
         UPF1_Nter_24h_R8h = mean(c(x200713, x200715, x200717)),
         UPF1_Nter_24h_R12h = mean(c(x200719, x200721, x200723)),
         UPF1_Nter_24h_R24h = mean(c(x200725, x200727, x200729)),
         UPF1_Nter_24h_R48h = mean(c(x200731, x200733, x200735))) %>% 
  dplyr::select(transcript_id,
                control_0h,
                control_48h,
                UPF1_Nter_0h,
                UPF1_Nter_2h,
                UPF1_Nter_4h,
                UPF1_Nter_8h,
                UPF1_Nter_12h,
                UPF1_Nter_24h,
                UPF1_Nter_48h,
                control_R0h,
                UPF1_Nter_24h_R0h,
                UPF1_Nter_24h_R2h,
                UPF1_Nter_24h_R4h,
                UPF1_Nter_24h_R8h,
                UPF1_Nter_24h_R12h,
                UPF1_Nter_24h_R24h,
                UPF1_Nter_24h_R48h,
                DTE_cluster_down) %>% 
  ungroup() %>% 
  mutate(DTE_cluster_down = (fct_relevel(DTE_cluster_down,
                                         "down 1:early",
                                         "down 2:delayed",
                                         "down 3:late",
                                         "down 4:inverse")))

### Save cluster information -------------------------------------------
save(catch,
     logcpm,
     Rev_1_F5_coldata,
     Rev_1_F5_transcript_cts_up, 
     Rev_1_F5_transcript_cts_down,
     Rev_1_F5_transcript_test_up_fix,
     Rev_1_F5_transcript_test_down_fix,
     Rev_1_F5_transcript_hclust_up,
     Rev_1_F5_transcript_hclust_down,
     Rev_1_F5_transcript_klus_cluster_up_df_fix,
     Rev_1_F5_transcript_klus_cluster_down_df_fix,
     Rev_1_F5_transcript_test_up_klus_cond,
     Rev_1_F5_transcript_test_down_klus_cond,
     Rev_1_F5_transcript_hclust_up_transcript_cluster_fix,
     Rev_1_F5_transcript_hclust_down_transcript_cluster_fix,
     NMDRHT.v1.2_tbl_DTE_cluster,
     file = paste0("Resources/Cluster/Rev_1_F5_transcript_cluster_information.rds"))

# load("Resources/Cluster/Rev_1_F5_transcript_cluster_information.rds")

##
# ImpulseDE2 ----------------------------------------------------------
##

### Samples---------------------------------------------------

Rev_1_F5_mydir1="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"
Rev_1_F5_mydir2="/home/volker/gencode.v42.datasets/2023_UPF1_AID_recovery_DW"

# Get samples and check if all files are present
Rev_1_F5_samples1 <- read.table(file.path(Rev_1_F5_mydir1, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) 

Rev_1_F5_samples2 <- read.table(file.path(Rev_1_F5_mydir2, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "Recovery_control",
                               TRUE ~ condition))

Rev_1_F5_samples <- as_tibble(bind_rows(Rev_1_F5_samples1,
                                        Rev_1_F5_samples2))

# Define Samples
Rev_1_F5_samples_for_ImpulseDE2 <- Rev_1_F5_samples %>% 
  mutate(Time = case_when(condition == "control" ~ 0,
                          condition == "control_48h" ~ 48,
                          condition == "UPF1_Nter_0h" ~ 0,
                          condition == "UPF1_Nter_2h" ~ 2,
                          condition == "UPF1_Nter_4h" ~ 4,
                          condition == "UPF1_Nter_8h" ~ 8,
                          condition == "UPF1_Nter_12h" ~ 12,
                          condition == "UPF1_Nter_24h" ~ 24,
                          condition == "UPF1_Nter_48h" ~ 48,
                          condition == "Recovery_control" ~ 0,
                          condition == "UPF1_Nter_24h_R0h" ~ 24,
                          condition == "UPF1_Nter_24h_R2h" ~ 26,
                          condition == "UPF1_Nter_24h_R4h" ~ 28,
                          condition == "UPF1_Nter_24h_R8h" ~ 32,
                          condition == "UPF1_Nter_24h_R12h" ~ 36,
                          condition == "UPF1_Nter_24h_R24h" ~ 48,
                          condition == "UPF1_Nter_24h_R48h" ~ 72)) %>% 
  mutate(Condition = case_when(condition == "control" ~ "control",
                               condition == "control_48h" ~ "control",
                               condition == "UPF1_Nter_0h" ~ "case",
                               condition == "UPF1_Nter_2h" ~ "case",
                               condition == "UPF1_Nter_4h" ~ "case",
                               condition == "UPF1_Nter_8h" ~ "case",
                               condition == "UPF1_Nter_12h" ~ "case",
                               condition == "UPF1_Nter_24h" ~ "case",
                               condition == "UPF1_Nter_48h" ~ "case",
                               condition == "Recovery_control" ~ "control",
                               condition == "UPF1_Nter_24h_R0h" ~ "case",
                               condition == "UPF1_Nter_24h_R2h" ~ "case",
                               condition == "UPF1_Nter_24h_R4h" ~ "case",
                               condition == "UPF1_Nter_24h_R8h" ~ "case",
                               condition == "UPF1_Nter_24h_R12h" ~ "case",
                               condition == "UPF1_Nter_24h_R24h" ~ "case",
                               condition == "UPF1_Nter_24h_R48h" ~ "case"),
         .after = sample) %>% 
  mutate(Batch = case_when(condition %in% c("control",
                                            "control_48h",
                                            "UPF1_Nter_0h",
                                            "UPF1_Nter_2h",
                                            "UPF1_Nter_4h",
                                            "UPF1_Nter_8h",
                                            "UPF1_Nter_12h",
                                            "UPF1_Nter_24h",
                                            "UPF1_Nter_48h") ~ "A",
                           condition %in% c("Recovery_control",
                                            "UPF1_Nter_24h_R0h",
                                            "UPF1_Nter_24h_R2h",
                                            "UPF1_Nter_24h_R4h",
                                            "UPF1_Nter_24h_R8h",
                                            "UPF1_Nter_24h_R12h",
                                            "UPF1_Nter_24h_R24h",
                                            "UPF1_Nter_24h_R48h") ~ "B")) %>% 
  dplyr::rename("Sample" = "sample") %>% 
  dplyr::select(-c(condition)) %>% 
  mutate(RowNames = Sample) %>% 
  column_to_rownames(var="RowNames")

# Filter Samples for just case samples, exclude 48h degradation time point
Rev_1_F5_samples_for_ImpulseDE2_Combined <- Rev_1_F5_samples_for_ImpulseDE2 %>% 
  filter(Condition == "case") %>% 
  filter(!Sample %in% c("191544",
                        "191546",
                        "191548"))

### Count data---------------------------------------------------

# Import pre-computed counts
load("Resources/Cluster/Rev_1_F5_transcript_cluster_information.rds")

# Import transcript-level counts 
catch_ImpulseDE2 <- catch

# Account for the mapping ambiguity
scaled.counts_ImpulseDE2 <- catch_ImpulseDE2$counts/catch_ImpulseDE2$annotation$Overdispersion

# Fix column names
colnames(scaled.counts_ImpulseDE2) <- Rev_1_F5_samples$sample

# Select only columns/conditions for ImpulseDE
scaled.counts_ImpulseDE2_filt <- scaled.counts_ImpulseDE2[, Rev_1_F5_samples_for_ImpulseDE2_Combined$Sample]

# Convert to integer
mode(scaled.counts_ImpulseDE2_filt) <- "integer"

# Run ImpulseDE2 - just case - Batch effects sensitive! - check for transient models
Rev_1_F5_ImpulseDE2_Combined <- ImpulseDE2::runImpulseDE2(
  matCountData    = scaled.counts_ImpulseDE2_filt, 
  dfAnnotation    = Rev_1_F5_samples_for_ImpulseDE2_Combined,
  boolCaseCtrl    = FALSE,
  vecConfounders  = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc        = 15 )

# Get results dataframe
Rev_1_F5_ImpulseDE2_Combined_df <- Rev_1_F5_ImpulseDE2_Combined$dfImpulseDE2Results

# Combine with annotation - define significance and best fit  
Rev_1_F5_ImpulseDE2_Combined_df_Cluster <- NMDRHT.v1.2_tbl_DTE_cluster %>% 
  left_join(Rev_1_F5_ImpulseDE2_Combined_df,
            by=c("transcript_id" = "Gene")) %>% 
  dplyr::select(-c("converge_impulse",
                   "converge_const",
                   "converge_sigmoid")) %>% 
  mutate(sigImpulseDE2 = case_when(padj < 0.0001 ~ "sig.",
                                   TRUE ~ "n.s."),
         .after = padj) %>% 
  mutate(bestFit=case_when(isMonotonous == TRUE ~ "Sigmoid",
                           isTransient == TRUE ~ "Impulse",
                           TRUE ~ "noFit"),
         .after = sigImpulseDE2)

# Obtain 1st Impulse parameters
Rev_1_F5_ImpulseDE2_Combined_Impulse_Parameter <- do.call(rbind, lapply(
  Rev_1_F5_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Impulse") %>% pull(transcript_id), function(x) {
    vecImpulseParam = get_lsModelFits(obj=Rev_1_F5_ImpulseDE2_Combined)[["case"]][[x]]$lsImpulseFit$vecImpulseParam
  }))

rownames(Rev_1_F5_ImpulseDE2_Combined_Impulse_Parameter) <- Rev_1_F5_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Impulse") %>% pull(transcript_id)

# Obtain 2nd Sigmoid parameters
Rev_1_F5_ImpulseDE2_Combined_Sigmoid_Parameter <- do.call(rbind, lapply(
  Rev_1_F5_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Sigmoid") %>% pull(transcript_id), function(x) {
    vecImpulseParam = get_lsModelFits(obj=Rev_1_F5_ImpulseDE2_Combined)[["case"]][[x]]$lsSigmoidFit$vecSigmoidParam
  }))

rownames(Rev_1_F5_ImpulseDE2_Combined_Sigmoid_Parameter) <- Rev_1_F5_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Sigmoid") %>% pull(transcript_id)

# Combine big dataframe with both information - add t2 and h2 (both NA) for sigmoid just for completeness
Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param <- Rev_1_F5_ImpulseDE2_Combined_df_Cluster %>% 
  left_join(bind_rows(as_tibble(Rev_1_F5_ImpulseDE2_Combined_Impulse_Parameter,
                                rownames="transcript_id"),
                      as_tibble(Rev_1_F5_ImpulseDE2_Combined_Sigmoid_Parameter,
                                rownames="transcript_id") %>% 
                        dplyr::rename("t1" = "t") %>% 
                        mutate(h2=as.numeric(NA),
                               .after=h1) %>% 
                        mutate(t2=as.numeric(NA),
                               .after=t1)
  ))

# How many significant transcripts
Rev_1_F5_ImpulseDE2_Combined_df %>% 
  filter(padj < 0.0001) %>% 
  dplyr::count()

#### Save ImpulseDE2 information -------------------------------------------
save(Rev_1_F5_ImpulseDE2_Combined_df,
     Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param,
     Rev_1_F5_ImpulseDE2_Combined,
     file = paste0("Resources/ImpulseDE2/Rev_1_F5_ImpulseDE2.rds"))

##
# NMD relevance - Cluster NMDRHT ----------------------------------------------------------
##

# Determine for each significantly up-/downregulated transcript in how many conditions it is significant
# Consider only concordant events (e.g. Up-cluster with up-regulated events | Down-cluster with down-regulated events)
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c(
    # "HEK293_UPF1_KD_Kishor_2019",
    "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
    # "HeLa_UPF1_KD_Longman_2020",
    # "HEK293_UPF1_KD_Fritz_2022",
    # "K562_UPF1_KD_Hug_2022",
    # "HUH7_UPF1_KD_Lee_2022",
    # "HepG2_UPF1_KD_He_2023",
    "HEK293_SMG567_KO_KD_Boehm_2021",
    "HEK293_SMG67_KD_Britto_Borges_2024",
    "HeLa_SMG67_KD_Britto_Borges_2024",
    "MCF7_SMG67_KD_Britto_Borges_2024",
    "U2OS_SMG67_KD_Britto_Borges_2024",
    "HEK293_UPF3_dKO_Wallmeroth_2022",
    # "HEK293_CASC3_KO_Gerbracht_2020",
    # "HCT116_UPF1_AID_degradation_this_Study",
    # "HCT116_UPF1_AID_recovery_this_Study",
    "HCT116_UPF1_FKBP_degradation_this_Study",
    "HEK293_UPF1_FKBP_degradation_this_Study",
    "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
    "HFF_SMG1i_this_Study",
    "HUVEC_SMG1i_this_Study")) %>% 
  filter(!condition_2 %in% c("SMG5_KD",
                             "SMG7_KO_2",
                             "SMG7_KO_34",
                             "UPF3BKO",
                             "control_01uM",
                             "SMG8_KO_0uM",
                             "SMG9_KO_0uM",
                             "SMG8_delKID_0uM",
                             "SMG8_KO_01uM",
                             "SMG9_KO_01uM",
                             "SMG8_delKID_01uM",
                             "SMG8_KO_1uM",
                             "SMG9_KO_1uM",
                             "SMG8_delKID_1uM",
                             "UPF1_FKBP_HEK_0h",
                             "UPF1_FKBP_HCT_0h")) %>% 
  filter(FDR < 0.0001) %>%
  filter(abs(logFC) > 1) %>% 
  mutate(UpDown = case_when(logFC > 1 & FDR < 0.0001 ~ "up",
                            logFC < -1 & FDR < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  left_join(NMDRHT.v1.2_tbl_DTE_cluster) %>% 
  filter((DTE_cluster %in% c("up 1:early",
                             "up 2:delayed",
                             "up 3:late",
                             "up 4:inverse") & UpDown == "up") | (DTE_cluster %in% c("down 1:early",
                                                                                     "down 2:delayed",
                                                                                     "down 3:late",
                                                                                     "down 4:inverse") & UpDown == "down" ) ) %>%  
  group_by(transcript_id, transcript_name, DTE_cluster, UpDown) %>% 
  summarize(NMD_n_sig_tx=n(),
            L2FC_median_NMD_tx = median(logFC),
            padj_median_NMD_tx = median(FDR)) %>% 
  ungroup() 

# Fill dataframe up with those transcripts that were never significant in any of the other NMD-compromised conditions
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_fill <- edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance %>% 
  bind_rows(NMDRHT.v1.2_tbl_DTE_cluster %>% 
              filter(DTE_cluster %in% c("up 1:early",
                                        "down 1:early",
                                        "up 2:delayed",
                                        "down 2:delayed",
                                        "up 3:late",
                                        "down 3:late",
                                        "up 4:inverse",
                                        "down 4:inverse"
              )) %>% 
              filter(!transcript_id %in% edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance$transcript_id) %>% 
              dplyr::select(transcript_id, transcript_name, DTE_cluster)) %>% 
  replace_na(list(NMD_n_sig_tx = 0, L2FC_median_NMD_tx = 0, padj_median_NMD_tx = 1)) %>% 
  mutate(UpDown=case_when(is.na(UpDown) & DTE_cluster %in% c("up 1:early",
                                                             "up 2:delayed",
                                                             "up 3:late",
                                                             "up 4:inverse") ~ "up",
                          is.na(UpDown) & DTE_cluster %in% c("down 4:inverse",
                                                             "down 3:late",
                                                             "down 2:delayed",
                                                             "down 1:early") ~ "down",
                          TRUE ~ UpDown)
  )

# Determine for each transcript the significant conditions in wide format
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_wide <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c(
    # "HEK293_UPF1_KD_Kishor_2019",
    "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
    # "HeLa_UPF1_KD_Longman_2020",
    # "HEK293_UPF1_KD_Fritz_2022",
    # "K562_UPF1_KD_Hug_2022",
    # "HUH7_UPF1_KD_Lee_2022",
    # "HepG2_UPF1_KD_He_2023",
    "HEK293_SMG567_KO_KD_Boehm_2021",
    "HEK293_SMG67_KD_Britto_Borges_2024",
    "HeLa_SMG67_KD_Britto_Borges_2024",
    "MCF7_SMG67_KD_Britto_Borges_2024",
    "U2OS_SMG67_KD_Britto_Borges_2024",
    "HEK293_UPF3_dKO_Wallmeroth_2022",
    # "HEK293_CASC3_KO_Gerbracht_2020",
    # "HCT116_UPF1_AID_degradation_this_Study",
    # "HCT116_UPF1_AID_recovery_this_Study",
    "HCT116_UPF1_FKBP_degradation_this_Study",
    "HEK293_UPF1_FKBP_degradation_this_Study",
    "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
    "HFF_SMG1i_this_Study",
    "HUVEC_SMG1i_this_Study")) %>% 
  filter(!condition_2 %in% c("SMG5_KD",
                             "SMG7_KO_2",
                             "SMG7_KO_34",
                             "UPF3BKO",
                             "control_01uM",
                             "SMG8_KO_0uM",
                             "SMG9_KO_0uM",
                             "SMG8_delKID_0uM",
                             "SMG8_KO_01uM",
                             "SMG9_KO_01uM",
                             "SMG8_delKID_01uM",
                             "SMG8_KO_1uM",
                             "SMG9_KO_1uM",
                             "SMG8_delKID_1uM",
                             "UPF1_FKBP_HEK_0h",
                             "UPF1_FKBP_HCT_0h")) %>% 
  mutate(condition_2 = case_when(condition_2 == "UPF1" ~ experimentSet,
                                 condition_2 == "SMG1i" ~ experimentSet,
                                 TRUE ~ condition_2)) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  filter(FDR < 0.0001) %>%
  filter(abs(logFC) > 1) %>% 
  mutate(UpDown = case_when(logFC > 1 & FDR < 0.0001 ~ "up",
                            logFC < -1 & FDR < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  dplyr::select(transcript_id, transcript_name, transcript_biotype, type, condition_2, UpDown) %>% 
  pivot_wider(names_from = condition_2, values_from = UpDown) 

# Join both dataframes
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete <- edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_fill %>% 
  left_join(edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_wide) %>% 
  relocate(transcript_biotype, type, .after=transcript_name) %>% 
  mutate(NMD_n_sig_tx_perc = 100*NMD_n_sig_tx/length(colnames(edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_wide %>% dplyr::select(-c(transcript_id, transcript_name, transcript_biotype, type)))),
         .after=NMD_n_sig_tx)

# Which transcripts are counted in up- and down?
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(transcript_id) %>% 
  dplyr::count() %>% 
  filter(n==2)

# How many transcripts per UpDown per NMD relevance
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  dplyr::count(UpDown, NMD_n_sig_tx_perc) %>% 
  arrange(desc(UpDown), desc(NMD_n_sig_tx_perc))

# Median NMD relevance per DGE_cluster
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  group_by(DTE_cluster) %>% 
  summarize(mean_NMD = mean(NMD_n_sig_tx_perc))

### NMD relevance bins ------------------------------------------------------

edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  dplyr::count(UpDown, DTE_cluster) %>% 
  arrange(desc(n))

# Per up-/downregulated -> cut into 4 equal width bins according to NMD relevance
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin <- edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late") & UpDown == "up" | DTE_cluster %in% c("down 1:early",
                                                                               "down 2:delayed",
                                                                               "down 3:late") & UpDown == "down"  ) %>%  
  group_by(UpDown) %>% 
  mutate(NMD_bin_tx = cut_width(NMD_n_sig_tx_perc, 25, boundary = 0), .after=NMD_n_sig_tx_perc) %>% 
  ungroup()

##
# NMD relevance - all transcripts ----------------------------------------------------------
##

# Obtain NMD relevance for *all* NMDRHT transcripts
NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c(
    # "HEK293_UPF1_KD_Kishor_2019",
    "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
    # "HeLa_UPF1_KD_Longman_2020",
    # "HEK293_UPF1_KD_Fritz_2022",
    # "K562_UPF1_KD_Hug_2022",
    # "HUH7_UPF1_KD_Lee_2022",
    # "HepG2_UPF1_KD_He_2023",
    "HEK293_SMG567_KO_KD_Boehm_2021",
    "HEK293_SMG67_KD_Britto_Borges_2024",
    "HeLa_SMG67_KD_Britto_Borges_2024",
    "MCF7_SMG67_KD_Britto_Borges_2024",
    "U2OS_SMG67_KD_Britto_Borges_2024",
    "HEK293_UPF3_dKO_Wallmeroth_2022",
    # "HEK293_CASC3_KO_Gerbracht_2020",
    # "HCT116_UPF1_AID_degradation_this_Study",
    # "HCT116_UPF1_AID_recovery_this_Study",
    "HCT116_UPF1_FKBP_degradation_this_Study",
    "HEK293_UPF1_FKBP_degradation_this_Study",
    "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
    "HFF_SMG1i_this_Study",
    "HUVEC_SMG1i_this_Study")) %>% 
  filter(!condition_2 %in% c("SMG5_KD",
                             "SMG7_KO_2",
                             "SMG7_KO_34",
                             "UPF3BKO",
                             "control_01uM",
                             "SMG8_KO_0uM",
                             "SMG9_KO_0uM",
                             "SMG8_delKID_0uM",
                             "SMG8_KO_01uM",
                             "SMG9_KO_01uM",
                             "SMG8_delKID_01uM",
                             "SMG8_KO_1uM",
                             "SMG9_KO_1uM",
                             "SMG8_delKID_1uM",
                             "UPF1_FKBP_HEK_0h",
                             "UPF1_FKBP_HCT_0h")) %>% 
  filter(FDR < 0.0001) %>%
  filter(abs(logFC) > 1) %>% 
  mutate(UpDown = case_when(logFC > 1 & FDR < 0.0001 ~ "up",
                            logFC < -1 & FDR < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  left_join(NMDRHT.v1.2_tbl_DTE_cluster) %>% 
  filter((DTE_cluster %in% c("up 1:early",
                             "up 2:delayed",
                             "up 3:late",
                             "up 4:inverse") & UpDown == "up") | 
           (DTE_cluster %in% c("down 1:early",
                               "down 2:delayed",
                               "down 3:late",
                               "down 4:inverse") & UpDown == "down" ) | 
           (DTE_cluster %in% c("expressed",
                               "not_expressed"))) %>% 
  group_by(transcript_id, transcript_name, DTE_cluster, UpDown) %>% 
  summarize(NMD_n_sig_tx=n(),
            L2FC_median_NMD_tx = median(logFC),
            padj_median_NMD_tx = median(FDR)) %>% 
  ungroup() 

# Fill dataframe up with those transcripts that were never significant in any of the other NMD-compromised conditions
NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill <- NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance %>% 
  bind_rows(NMDRHT.v1.2_tbl_DTE_cluster %>% 
              filter(!transcript_id %in% NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance$transcript_id) %>% 
              dplyr::select(transcript_id, transcript_name, DTE_cluster)) %>% 
  replace_na(list(NMD_n_sig_tx = 0, median_log2FC = 0, median_FDR = 1)) %>% 
  group_by(transcript_id, DTE_cluster) %>% 
  slice_max(NMD_n_sig_tx, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(UpDown=case_when(is.na(UpDown) & DTE_cluster %in% c("up 1:early",
                                                             "up 2:delayed",
                                                             "up 3:late",
                                                             "up 4:inverse") ~ "up",
                          is.na(UpDown) & DTE_cluster %in% c("down 4:inverse",
                                                             "down 3:late",
                                                             "down 2:delayed",
                                                             "down 1:early") ~ "down",
                          is.na(UpDown) ~ "n.s.",
                          TRUE ~ UpDown)
  ) %>%  
  mutate(NMD_n_sig_tx_perc = 100*NMD_n_sig_tx/max(NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance$NMD_n_sig_tx),
         .after=NMD_n_sig_tx)

# Determine NMD bin
NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill_bin <- NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill %>% 
  filter(UpDown != "n.s.") %>% 
  group_by(UpDown) %>% 
  mutate(NMD_bin_tx = cut_width(NMD_n_sig_tx_perc, 25, boundary = 0), .after=NMD_n_sig_tx_perc) %>% 
  ungroup() %>% 
  bind_rows(NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill %>% 
              filter(UpDown == "n.s.") %>% 
              mutate(NMD_bin_tx = "n.s."))

### Save NMD relevance information -------------------------------------------
save(edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete,
     edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin,
     NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill_bin,
     file = paste0("Resources/NMD_relevance/","edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance.rds"))

# Save as csv
NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill_bin %>% 
  write_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill_bin.csv")

##
# Transcript-MainTable ---------------------------------------------------
##

# Load required data
NMDRHT.v1.2_tbl_NMD <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl_NMD.csv")
load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_2_datasources.rds")
load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_3_datasources.rds")

# MainTable
NMDRHT.v1.2_MainTable <- NMDRHT.v1.2_tbl_DTE_cluster %>% 
  # remove unnecessary information (mostly from NMDRHT annotation generation) - can be recovered by joining the respective tbl
  dplyr::select(-c(GENCODE_gene_id,
                   GENCODE_gene_name,
                   GENCODE_transcript_name,
                   GENCODE_gene_type,
                   GENCODE_transcript_type,
                   GTF_transcript_id,
                   study,
                   cell_line,
                   method,
                   sequencing,
                   annotation,
                   gene_id_conflict,
                   ORFquant_Distance_to_lastExEx,
                   ORFquant_NMD_tx_status_combined,
                   ORFquant_NMD_reason,
                   ORFanage_status,
                   ORFanage_duplicity,
                   ORFanage_template_source,
                   ORFanage_NMD_status,
                   ORFanage_stop_to_lastEJ,
                   ORFanage_num_of_downEJs,
                   ORFanage_UTR3_length,
                   utr3_length_str
                   )) %>% 
  left_join(NMDRHT.v1.2_tbl_DTE_cluster_NMD_relevance_fill_bin) %>% 
  left_join(Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
              dplyr::select(transcript_id,padj,sigImpulseDE2, bestFit, beta, h0, h1, h2, t1, t2) %>% 
              dplyr::rename("bestFit_ImpulseDE2_tx" = "bestFit",
                            "padj_ImpulseDE2_tx" = "padj")) %>% 
  left_join(NMDRHT.v1.2_tbl_NMD %>% dplyr::select(transcript_id, NMD_50nt_rule, stop_to_lastEJ, num_of_downEJs, UTR3_length)) %>% 
  left_join(NMDRHT_tracking_filtered %>% 
              dplyr::select(NMDRHT_transcript_id,UIC_filter_level) %>% 
              dplyr::rename("transcript_id"="NMDRHT_transcript_id")) %>% 
  left_join(NMDRHT_final_selection %>% 
              dplyr::select(NMDRHT_transcript_id, dist_to_CAGE_peak, dist_to_polyA_site) %>% 
              dplyr::rename("transcript_id"="NMDRHT_transcript_id")) %>% 
  relocate(dist_to_CAGE_peak, dist_to_polyA_site, .after=within_polyA_site) %>% 
  relocate(UIC_filter_level, .after=UIC_SR_support) %>% 
  relocate("NMD_50nt_rule",
           "stop_to_lastEJ",
           "num_of_downEJs",
           "UTR3_length", .after=NMD_tx_reason) %>% 
  relocate(DTE_cluster_up, DTE_cluster_down, .after = "DTE_cluster") 

# Export NMDRHT main table
NMDRHT.v1.2_MainTable %>% 
  write_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

save(NMDRHT.v1.2_MainTable,
     file="Resources/NMDRHT/NMDRHT.v1.2_MainTable.rds")

# Read in - if necessary
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

load("Resources/NMDRHT/NMDRHT.v1.2_MainTable.rds")

##
# PARE analysis (Schmidt et al. 2015) --------------------------------------------
##

# Import GTF file
gtf <- rtracklayer::import(file.path("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), format="gtf")

# Generate TxDb from annotation
NMDRHT_txdb1 <- makeTxDbFromGRanges(gtf, drop.stop.codons=FALSE)

# Load ORFik package
library(ORFik)

# Define transcript regions
mrna <- loadRegion(NMDRHT_txdb1, "mrna")

# Extract regions around stop codons (99nt upstream, 100nt downstream)
NMDRHT_txdb1_200nt_stop <- stopRegion(loadRegion(NMDRHT_txdb1, "cds"), mrna, upstream = 99, downstream = 100, is.sorted = TRUE)

# Check width distribution
as_tibble(sum(width(NMDRHT_txdb1_200nt_stop)), rownames = "transcript_id") %>% 
  ggplot(aes(x=value)) +
  geom_histogram()

# Export as unlisted bed file
export(NMDRHT_txdb1_200nt_stop@unlistData, "Resources/NMDRHT/NMDRHT.v1.2_stop_200nt_region.bed")

##
## Count PARE reads in 200nt window around stop ----------------------------
##

# Convert Bam to ofst - done only once
ofst_out_dir <- file.path("/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/ofst")
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574720.bam", ofst_out_dir)
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574721.bam", ofst_out_dir)
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574722.bam", ofst_out_dir)
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574726.bam", ofst_out_dir)
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574727.bam", ofst_out_dir)
# convert_bam_to_ofst(NULL, "/home/volker/gencode.v42.datasets/2014_Endocleavage_Schmidt/Fragments/BAM/SRR1574728.bam", ofst_out_dir)

# Find the file again
ofst_files <- list.files(ofst_out_dir, full.names = TRUE)

##
# Load it
##

# XRN1_rep1
SRR1574720_ofst <- fimport(ofst_files[1])
SRR1574720_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574720_ofst)

# XRN1_rep2
SRR1574726_ofst <- fimport(ofst_files[4])
SRR1574726_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574726_ofst)

# XRN1_SMG6_rep1
SRR1574721_ofst <- fimport(ofst_files[2])
SRR1574721_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574721_ofst)

# XRN1_SMG6_rep2
SRR1574727_ofst <- fimport(ofst_files[5])
SRR1574727_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574727_ofst)

# XRN1_UPF1_rep1
SRR1574722_ofst <- fimport(ofst_files[3])
SRR1574722_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574722_ofst)

# XRN1_UPF1_rep2
SRR1574728_ofst <- fimport(ofst_files[6])
SRR1574728_200nt_stop <- countOverlapsW(NMDRHT_txdb1_200nt_stop, SRR1574728_ofst)

##
## Combine outputs ---------------------------------------------------------
##

PARE_analysis_combined <- as_tibble(SRR1574720_200nt_stop, rownames="transcript_id") %>% 
  dplyr::rename("XRN1_rep1" = "value") %>% 
  left_join(as_tibble(SRR1574726_200nt_stop, rownames="transcript_id") %>% 
              dplyr::rename("XRN1_rep2" = "value")) %>% 
  left_join(as_tibble(SRR1574721_200nt_stop, rownames="transcript_id") %>% 
              dplyr::rename("XRN1_SMG6_rep1" = "value"))  %>% 
  left_join(as_tibble(SRR1574727_200nt_stop, rownames="transcript_id") %>% 
              dplyr::rename("XRN1_SMG6_rep2" = "value")) %>% 
  left_join(as_tibble(SRR1574722_200nt_stop, rownames="transcript_id") %>% 
              dplyr::rename("XRN1_UPF1_rep1" = "value"))  %>% 
  left_join(as_tibble(SRR1574728_200nt_stop, rownames="transcript_id") %>% 
              dplyr::rename("XRN1_UPF1_rep2" = "value"))

PARE_analysis_combined_scaled <- PARE_analysis_combined %>% 
  mutate(XRN1_rep1_CPM = XRN1_rep1/length(SRR1574720_ofst)*10^6,
         XRN1_rep2_CPM = XRN1_rep2/length(SRR1574726_ofst)*10^6,
         XRN1_SMG6_rep1_CPM = XRN1_SMG6_rep1/length(SRR1574721_ofst)*10^6,
         XRN1_SMG6_rep2_CPM = XRN1_SMG6_rep2/length(SRR1574727_ofst)*10^6,
         XRN1_UPF1_rep1_CPM = XRN1_UPF1_rep1/length(SRR1574722_ofst)*10^6,
         XRN1_UPF1_rep2_CPM = XRN1_UPF1_rep2/length(SRR1574728_ofst)*10^6) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name)) %>% 
  relocate(gene_id, gene_name, transcript_id, transcript_name)

### Save PARE analysis information -------------------------------------------
save(SRR1574720_ofst,
     SRR1574726_ofst,
     SRR1574721_ofst,
     SRR1574727_ofst,
     SRR1574722_ofst,
     SRR1574728_ofst,
     PARE_analysis_combined,
     PARE_analysis_combined_scaled,
     file = paste0("Resources/PARE/","PARE_analysis.rds"))

##
## edgeR PARE analysis -----------------------------------------------------
##

# Define DGEList
PARE_y <- edgeR::DGEList(PARE_analysis_combined,
                  group=c("XRN1","XRN1","XRN1_SMG6","XRN1_SMG6","XRN1_UPF1","XRN1_UPF1"))

# Supply total mapped reads as lib.size - important since coverage was only computed for +/- 100 nt around stop!
PARE_y$samples$lib.size <- c(length(SRR1574720_ofst),
                             length(SRR1574726_ofst),
                             length(SRR1574721_ofst),
                             length(SRR1574727_ofst),
                             length(SRR1574722_ofst),
                             length(SRR1574728_ofst))

# Filter out low "expressed" stop codon posiitons
PARE_keep <- edgeR::filterByExpr(PARE_y)

PARE_y_filt <- PARE_y[PARE_keep, , keep.lib.sizes=TRUE]

# Normalize data
PARE_y_filt_norm <- edgeR::normLibSizes(PARE_y_filt)

# Desing matrix
design <- model.matrix(~ 0 + group,data = PARE_y_filt_norm$samples)
colnames(design) <- gsub("group", "", colnames(design))

# Dispersion estimation
PARE_y_filt_norm_NB <- edgeR::estimateDisp(PARE_y_filt_norm, design, robust=TRUE)
print("Common dispersion")
PARE_y_filt_norm_NB$common.dispersion

# Plot BCV
edgeR::plotBCV(PARE_y_filt_norm_NB)

# Quasi-likelihood (QL) pipeline
PARE_fit <- edgeR::glmQLFit(PARE_y_filt_norm_NB, design, robust=TRUE)

# Plot QLDisp
edgeR::plotQLDisp(PARE_fit)

##
### XRN1+SMG6 ---------------------------------------------------------------
##

# Make contrasts
contrast_XS6 <- makeContrasts(paste("XRN1", "-", "XRN1_SMG6"),levels=design)

# Test using QL F-test
qlf_XS6 <- edgeR::glmQLFTest(PARE_fit, contrast=contrast_XS6)

# Obtain summary
is.de_XS6 <- decideTests(qlf_XS6, p.value=0.1)
summary(is.de_XS6)

# Get data frame and join with annotation
qlf_XS6_df <- as.data.frame(topTags(qlf_XS6, n = Inf, sort.by = "none")) %>% 
  # rownames_to_column(var="transcript_id") %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name)) %>% 
  relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  mutate(condition_2 = "XRN1_SMG6")

##
### XRN1+UPF1 ---------------------------------------------------------------
##

# Make contrasts
contrast_XU1 <- makeContrasts(paste("XRN1", "-", "XRN1_UPF1"),levels=design)

# Test using QL F-test
qlf_XU1 <- edgeR::glmQLFTest(PARE_fit, contrast=contrast_XU1)

# Obtain summary
is.de_XU1 <- decideTests(qlf_XU1, p.value=0.1)
summary(is.de_XU1)

# Get data frame and join with annotation
qlf_XU1_df <- as.data.frame(topTags(qlf_XU1, n = Inf, sort.by = "none")) %>% 
  # rownames_to_column(var="transcript_id") %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name)) %>% 
  relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  mutate(condition_2 = "XRN1_UPF1")

###
### Combined ----------------------------------------------------------------
###

PARE_edgeR_combined <- qlf_XS6_df %>% 
  dplyr::select(gene_id, gene_name, transcript_id, transcript_name, logFC, logCPM) %>% 
  dplyr::rename("L2FC_XS6" = "logFC",
                "logCPM_XS6" = "logCPM") %>% 
  left_join(qlf_XU1_df %>% 
              dplyr::select(gene_id, gene_name, transcript_id, transcript_name, logFC, logCPM) %>% 
              dplyr::rename("L2FC_XU1" = "logFC",
                            "logCPM_XU1" = "logCPM")) %>% 
  left_join(as_tibble(PARE_fit$fitted.values) %>% mutate(transcript_id = PARE_fit$genes$transcript_id)) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(transcript_id, NMD_tx_status, NMD_tx_reason, NMD_50nt_rule, NMD_bin_tx, DTE_cluster, UpDown))

# Save to disk 
PARE_edgeR_combined %>% write_csv("Resources/PARE/PARE_edgeR_combined.csv")

##
# Imamachi 2017 -----------------------------------------------------------
##

## Table 3 - Parameters ----------------------------------------------------

# Initialize bioMart | GENCODE v.42 is equal to ensembl version 108
ensembl_108 = biomaRt::useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 108)

# Data from Supplemental Table S3 | expand refseq column
Imamachi_Supplemental_Table_S3 <- read_excel("Resources/External/Imamachi2017/Supplemental_Table_S3.xlsx", 
                                             skip = 2) %>% 
  clean_names() %>% 
  dplyr::select(gene_symbol, 
                ref_seq_id_human, 
                predicted_upf1_motifs_top4_ccugggg_ccuggga_ccuggaa_ccugaga,
                predicted_upf1_motifs_top1_ccugggg,
                fc_rna_half_life_log2_si_upf1_si_ctrl_bric_seq_tani_et_al,
                fold_enrichment_log2_ha_upf1_ip_input_rip_seq,
                fold_enrichment_log2_p_upf1_ip_ctrl_rna_footprint_kurosaki_et_al,
                fc_rna_half_life_log2_si_stau1_si_ctrl_bric_seq_maekawa_et_al) %>% 
  separate_rows(ref_seq_id_human, sep =",", convert = TRUE)

# How many distinct gene_names = 8426
Imamachi_Supplemental_Table_S3 %>% 
  distinct(gene_symbol) %>% 
  dplyr::count()


# Perform bioMart request
Imamachi_Supplemental_Table_S3_ENSG <- getBM(c("ensembl_gene_id_version","hgnc_symbol", "refseq_mrna"), "refseq_mrna", Imamachi_Supplemental_Table_S3$ref_seq_id_human, ensembl_108)

# Merge and manually select the matching GENCODE gene id
Imamachi_Supplemental_Table_S3_ENSG_distinct <- Imamachi_Supplemental_Table_S3_ENSG %>% 
  left_join(Imamachi_Supplemental_Table_S3,
            by=c("refseq_mrna" = "ref_seq_id_human")) %>% 
  # distinct(ensembl_gene_id_version, .keep_all = TRUE) %>%
  dplyr::rename("gene_id" = "ensembl_gene_id_version",
                "gene_name" = "hgnc_symbol") %>%
  distinct(gene_id, gene_name, .keep_all = TRUE) %>% 
  group_by(gene_symbol) %>% 
  mutate(n_hits = n()) %>% 
  ungroup()

# check for available GENCODE information and filter for those with information | check again for multi-hits
Imamachi_Supplemental_Table_S3_ENSG_distinct_filt <- Imamachi_Supplemental_Table_S3_ENSG_distinct %>% 
  left_join(gtf_gencode_df_short_DGE_cluster %>% dplyr::select(gene_id, gene_name, gene_type)) %>% 
  filter(!is.na(gene_type)) %>% 
  group_by(gene_name) %>% 
  mutate(n_hits = n()) %>% 
  ungroup()

##### Join with global GENCODE gene-level table  ----------------------------------------------------
Imamachi_Supplemental_Table_S3_GENCODE <- Imamachi_Supplemental_Table_S3_ENSG_distinct_filt %>% 
  left_join(GENCODE_v42_MainTable) %>% 
  dplyr::rename("upf1_motifs_top4" = "predicted_upf1_motifs_top4_ccugggg_ccuggga_ccuggaa_ccugaga",
                "upf1_motifs_top1" = "predicted_upf1_motifs_top1_ccugggg",
                "L2FC_hl_siUPF1_Tani" = "fc_rna_half_life_log2_si_upf1_si_ctrl_bric_seq_tani_et_al",
                "L2FC_UPF1_RIP" = "fold_enrichment_log2_ha_upf1_ip_input_rip_seq",
                "L2FC_pUPF1_Kurosaki" = "fold_enrichment_log2_p_upf1_ip_ctrl_rna_footprint_kurosaki_et_al",
                "L2FC_hl_siSTAU1_Maekawa" = "fc_rna_half_life_log2_si_stau1_si_ctrl_bric_seq_maekawa_et_al") %>% 
  dplyr::select(-n_hits)

# Save as csv
Imamachi_Supplemental_Table_S3_GENCODE %>% write_csv("Resources/External/Imamachi2017/Imamachi_Supplemental_Table_S3_GENCODE.csv")

#
##
###
# Boruta RF ---------------------------------------------------------------
###
##
#

# Objective: determine important parameters that explain upregulation of seemingly non-NMD transcripts
# Approach: random forest classification-based feature selection using the Boruta R package 

# Load Boruta package
library(Boruta)

# Load required data
NMDRHT_utr3_GC_200nt <- read_csv("Resources/NMDRHT/NMDRHT_utr3_GC_200nt.csv")
Imamachi_Supplemental_Table_S3_GENCODE <- read_csv("Resources/External/Imamachi2017_PMID_27940950/Imamachi_Supplemental_Table_S3_GENCODE.csv")
Kurosaki2014_tbl <- read_csv("Resources/External/Kurosaki2014_PMID_25184677/Kurosaki2014_tbl.csv")
GSE84722_ST3 <- read_csv("Resources/External/Mukherjee2017_PMID_27870833/GSE84722_ST3.csv")

# Define and select parameters to test
NMDRHT.v1.2_MainTable_forBoruta <- NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>% 
  left_join(NMDRHT_utr3_GC_200nt %>% dplyr::rename("transcript_id" = "tx_id")) %>% 
  left_join(Imamachi_Supplemental_Table_S3_GENCODE %>% 
              dplyr::select(gene_id, gene_name, upf1_motifs_top4, upf1_motifs_top1, L2FC_hl_siUPF1_Tani, L2FC_UPF1_RIP, L2FC_pUPF1_Kurosaki, L2FC_hl_siSTAU1_Maekawa)) %>% 
  left_join(Kurosaki2014_tbl %>% filter(condition_2 == "p_UPF1") %>% dplyr::select(gene_id, log2FoldChange) %>% dplyr::rename("L2FC_pUPF1_Kurosaki_new" = "log2FoldChange")) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                      L2FC_kdeg,
                                                      padj_kdeg,
                                                      L2FC_ksyn,
                                                      padj_ksyn,
                                                      Mech_score,
                                                      kdeg_conclusion,
                                                      RNA_conclusion,
                                                      Mech_conclusion)) %>% 
  left_join(GSE84722_ST3 %>% dplyr::select(gene_id, Syn, Proc, Deg, CytNuc, PolyCyt, TrP, Copies)) %>% 
  mutate(utr3_mfe_nt = -utr3_mfe_nt) %>% 
  dplyr::select(gene_name,
                gene_id,
                transcript_name,
                DTE_cluster,
                logFC,
                NMD_50nt_rule,
                stop_to_lastEJ,
                num_of_downEJs,
                ORFquant_log2FC_tx_ORF_pct_P_sites_pN,
                nexon,
                tx_len,
                cds_len,
                utr5_len,
                utr3_len,
                GC_tx,
                GC_utr5,
                GC_cds,
                GC_utr3,
                GC_utr3_200nt,
                utr3_mfe_nt,
                # NMD_bin,
                NMD_n_sig_tx_perc,
                upf1_motifs_top4,
                upf1_motifs_top1,
                L2FC_hl_siUPF1_Tani,
                L2FC_kdeg,
                Mech_score,
                Mech_conclusion,
                L2FC_UPF1_RIP,
                L2FC_pUPF1_Kurosaki_new,
                Syn,
                Proc,
                Deg,
                CytNuc,
                PolyCyt,
                TrP,
                Copies) %>% 
  mutate(
    # NMD_bin = as_factor(NMD_bin),
    NMD_50nt_rule = as_factor(NMD_50nt_rule),
    DTE_cluster = as_factor(DTE_cluster),
    Mech_conclusion = as_factor(Mech_conclusion))

# Get overview of individual parameters (important: number of NAs)
summary(NMDRHT.v1.2_MainTable_forBoruta)

# Comment: we define three parameter sets:
# Full = all parameters
# Reduced = removing parameters with many missing values
# Minimal = most essential parameters only

###
## early up - all - FULL PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_all <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule %in% c("FALSE", "TRUE")) %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_all)

# Perform Boruta!
set.seed(123)
boruta_fullParam_earlyUP_all_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_all), doTrace=3)

# Output
print(boruta_fullParam_earlyUP_all_output)

# extract Boruta's Importance history
boruta_fullParam_earlyUP_all_output_ImpHis <- as_tibble(boruta_fullParam_earlyUP_all_output$ImpHistory)

# extract Boruta's Decision
boruta_fullParam_earlyUP_all_output_Decision <- as_tibble(boruta_fullParam_earlyUP_all_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_fullParam_earlyUP_all_output_ImpHis_wide <- boruta_fullParam_earlyUP_all_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_fullParam_earlyUP_all_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_fullParam_earlyUP_all_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE+FALSE\nfull parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(34*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_all_fullParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - all - Reduced PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_all <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule %in% c("FALSE", "TRUE")) %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Copies, L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_all)

# Perform Boruta!
set.seed(123)
boruta_redParam_earlyUP_all_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_all), doTrace=3)

# Output
print(boruta_redParam_earlyUP_all_output)

# extract Boruta's Importance history
boruta_redParam_earlyUP_all_output_ImpHis <- as_tibble(boruta_redParam_earlyUP_all_output$ImpHistory)

# extract Boruta's Decision
boruta_redParam_earlyUP_all_output_Decision <- as_tibble(boruta_redParam_earlyUP_all_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_redParam_earlyUP_all_output_ImpHis_wide <- boruta_redParam_earlyUP_all_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_redParam_earlyUP_all_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_redParam_earlyUP_all_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE+FALSE\nreduced parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(28*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_all_reducedParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - all - Minimal PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_all <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule %in% c("FALSE", "TRUE")) %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Syn,
                   Proc,
                   Deg,
                   CytNuc,
                   PolyCyt,
                   TrP,
                   Copies,
                   L2FC_kdeg,
                   Mech_score,
                   Mech_conclusion,
                   L2FC_pUPF1_Kurosaki_new,L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_all)

# Perform Boruta!
set.seed(123)
boruta_minParam_earlyUP_all_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_all), doTrace=3)

# Output
print(boruta_minParam_earlyUP_all_output)

# extract Boruta's Importance history
boruta_minParam_earlyUP_all_output_ImpHis <- as_tibble(boruta_minParam_earlyUP_all_output$ImpHistory)

# extract Boruta's Decision
boruta_minParam_earlyUP_all_output_Decision <- as_tibble(boruta_minParam_earlyUP_all_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_minParam_earlyUP_all_output_ImpHis_wide <- boruta_minParam_earlyUP_all_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_minParam_earlyUP_all_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_minParam_earlyUP_all_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE+FALSE\nminimal parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(18*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_all_minimalParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - Non-NMD - FULL PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_no50nt <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "FALSE") %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_no50nt)

# Perform Boruta!
set.seed(123)
boruta_fullParam_earlyUP_no50nt_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_no50nt), doTrace=3)

# Output
print(boruta_fullParam_earlyUP_no50nt_output)

# extract Boruta's Importance history
boruta_fullParam_earlyUP_no50nt_output_ImpHis <- as_tibble(boruta_fullParam_earlyUP_no50nt_output$ImpHistory)

# extract Boruta's Decision
boruta_fullParam_earlyUP_no50nt_output_Decision <- as_tibble(boruta_fullParam_earlyUP_no50nt_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_fullParam_earlyUP_no50nt_output_ImpHis_wide <- boruta_fullParam_earlyUP_no50nt_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_fullParam_earlyUP_no50nt_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_fullParam_earlyUP_no50nt_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=none\nfull parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(34*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_noNMD_fullParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - Non-NMD - Reduced PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_no50nt <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "FALSE") %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Copies, L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_no50nt)

# Perform Boruta!
set.seed(123)
boruta_redParam_earlyUP_no50nt_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_no50nt), doTrace=3)

# Output
print(boruta_redParam_earlyUP_no50nt_output)

# extract Boruta's Importance history
boruta_redParam_earlyUP_no50nt_output_ImpHis <- as_tibble(boruta_redParam_earlyUP_no50nt_output$ImpHistory)

# extract Boruta's Decision
boruta_redParam_earlyUP_no50nt_output_Decision <- as_tibble(boruta_redParam_earlyUP_no50nt_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_redParam_earlyUP_no50nt_output_ImpHis_wide <- boruta_redParam_earlyUP_no50nt_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_redParam_earlyUP_no50nt_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_redParam_earlyUP_no50nt_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=none\nreduced parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(28*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_noNMD_reducedParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - Non-NMD - Minimal PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_no50nt <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "FALSE" & Mech_conclusion == "Degradation") %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Syn,
                   Proc,
                   Deg,
                   CytNuc,
                   PolyCyt,
                   TrP,
                   Copies,
                   L2FC_kdeg,
                   Mech_score,
                   Mech_conclusion,
                   L2FC_pUPF1_Kurosaki_new,L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_no50nt)

# Perform Boruta!
set.seed(123)
boruta_minParam_earlyUP_no50nt_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_no50nt), doTrace=3)

# Output
print(boruta_minParam_earlyUP_no50nt_output)

# extract Boruta's Importance history
boruta_minParam_earlyUP_no50nt_output_ImpHis <- as_tibble(boruta_minParam_earlyUP_no50nt_output$ImpHistory)

# extract Boruta's Decision
boruta_minParam_earlyUP_no50nt_output_Decision <- as_tibble(boruta_minParam_earlyUP_no50nt_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_minParam_earlyUP_no50nt_output_ImpHis_wide <- boruta_minParam_earlyUP_no50nt_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_minParam_earlyUP_no50nt_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_minParam_earlyUP_no50nt_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=none\nminimal parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(18*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_noNMD_minimalParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - NMD - FULL PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_NMD <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "TRUE") %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_NMD)

# Perform Boruta!
set.seed(123)
boruta_fullParam_earlyUP_NMD_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_fullParam_earlyUp_NMD), doTrace=3)

# Output
print(boruta_fullParam_earlyUP_NMD_output)

# extract Boruta's Importance history
boruta_fullParam_earlyUP_NMD_output_ImpHis <- as_tibble(boruta_fullParam_earlyUP_NMD_output$ImpHistory)

# extract Boruta's Decision
boruta_fullParam_earlyUP_NMD_output_Decision <- as_tibble(boruta_fullParam_earlyUP_NMD_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_fullParam_earlyUP_NMD_output_ImpHis_wide <- boruta_fullParam_earlyUP_NMD_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_fullParam_earlyUP_NMD_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_fullParam_earlyUP_NMD_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE\nfull parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(34*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_onlyNMD_fullParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - NMD - Reduced PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_NMD <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "TRUE") %>% 
  
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Copies, L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_NMD)

# Perform Boruta!
set.seed(123)
boruta_redParam_earlyUP_NMD_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_redParam_earlyUp_NMD), doTrace=3)

# Output
print(boruta_redParam_earlyUP_NMD_output)

# extract Boruta's Importance history
boruta_redParam_earlyUP_NMD_output_ImpHis <- as_tibble(boruta_redParam_earlyUP_NMD_output$ImpHistory)

# extract Boruta's Decision
boruta_redParam_earlyUP_NMD_output_Decision <- as_tibble(boruta_redParam_earlyUP_NMD_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_redParam_earlyUP_NMD_output_ImpHis_wide <- boruta_redParam_earlyUP_NMD_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_redParam_earlyUP_NMD_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_redParam_earlyUP_NMD_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE\nreduced parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(28*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_onlyNMD_reducedParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## early up - NMD - Minimal PARAMETER SET  ------------------------------------
###

# Select transcripts
NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_NMD <- NMDRHT.v1.2_MainTable_forBoruta %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == "TRUE" & Mech_conclusion == "Degradation") %>% 
  dplyr::select(-c(gene_id, DTE_cluster, gene_name, transcript_name, Syn,
                   Proc,
                   Deg,
                   CytNuc,
                   PolyCyt,
                   TrP,
                   Copies,
                   L2FC_kdeg,
                   Mech_score,
                   Mech_conclusion,
                   L2FC_pUPF1_Kurosaki_new,L2FC_UPF1_RIP, L2FC_hl_siUPF1_Tani, upf1_motifs_top4, upf1_motifs_top1,ORFquant_log2FC_tx_ORF_pct_P_sites_pN))

# Get summary of parameters
summary(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_NMD)

# Perform Boruta!
set.seed(123)
boruta_minParam_earlyUP_NMD_output <- Boruta(logFC ~ ., data=na.omit(NMDRHT.v1.2_MainTable_forBoruta_minParam_earlyUp_NMD), doTrace=3)

# Output
print(boruta_minParam_earlyUP_NMD_output)

# extract Boruta's Importance history
boruta_minParam_earlyUP_NMD_output_ImpHis <- as_tibble(boruta_minParam_earlyUP_NMD_output$ImpHistory)

# extract Boruta's Decision
boruta_minParam_earlyUP_NMD_output_Decision <- as_tibble(boruta_minParam_earlyUP_NMD_output$finalDecision, rownames = "parameter")

# Join both - pivot longer and tidy up
boruta_minParam_earlyUP_NMD_output_ImpHis_wide <- boruta_minParam_earlyUP_NMD_output_ImpHis %>% 
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "importance") %>% 
  left_join(boruta_minParam_earlyUP_NMD_output_Decision) %>% 
  dplyr::rename("decision" = "value") %>% 
  mutate(decision = case_when(!is.na(decision) ~ decision,
                              TRUE ~ "Shadow")) %>% 
  group_by(parameter) %>% 
  mutate(median_imp = median(importance)) %>% 
  ungroup() %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "shadowMax",
                                 "shadowMean",
                                 "shadowMin",
                                 after = Inf)) %>% 
  mutate(decision = fct_relevel(decision,
                                "Confirmed",
                                "Tentative",
                                "Rejected",
                                "Shadow"))

# Plot
boruta_minParam_earlyUP_NMD_output_ImpHis_wide %>%  
  ggplot(aes(x=parameter,
             y=importance,
             fill=decision)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("Confirmed" = "#5E8C61",
                             "Tentative" = "#80A1C1",
                             "Rejected" = "#EEE3AB",
                             "Shadow" = "#C94277")) +
  labs(title="cluster 1:up | NMD reason=TRUE\nminimal parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(18*2+0.4, "mm"))

ggsave(file.path("Resources/Boruta", "Rev_1_Boruta_earlyUp_onlyNMD_minimalParameter_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#
## Combined analyses -------------------------------------------------------
#

combined_boruta_analysis <- bind_rows(boruta_fullParam_earlyUP_all_output_ImpHis_wide %>% 
                                        mutate(analysis = "all_FP")) %>% 
  bind_rows(boruta_redParam_earlyUP_all_output_ImpHis_wide %>% 
              mutate(analysis = "all_RP")) %>% 
  bind_rows(boruta_minParam_earlyUP_all_output_ImpHis_wide %>% 
              mutate(analysis = "all_MP")) %>% 
  bind_rows(boruta_fullParam_earlyUP_no50nt_output_ImpHis_wide %>% 
              mutate(analysis = "noNMD_FP")) %>% 
  bind_rows(boruta_redParam_earlyUP_no50nt_output_ImpHis_wide %>% 
              mutate(analysis = "noNMD_RP")) %>% 
  bind_rows(boruta_minParam_earlyUP_no50nt_output_ImpHis_wide %>% 
              mutate(analysis = "noNMD_MP")) %>% 
  bind_rows(boruta_fullParam_earlyUP_NMD_output_ImpHis_wide %>% 
              mutate(analysis = "NMD_FP")) %>% 
  bind_rows(boruta_redParam_earlyUP_NMD_output_ImpHis_wide %>% 
              mutate(analysis = "NMD_RP")) %>% 
  bind_rows(boruta_minParam_earlyUP_NMD_output_ImpHis_wide %>% 
              mutate(analysis = "NMD_MP")) %>% 
  mutate(parameter = case_when(parameter == "NMD_50nt_rule" ~ "50-nucleotide rule",
                               parameter == "ORFquant_log2FC_tx_ORF_pct_P_sites_pN" ~ "log2FC(P-sites)",
                               parameter == "num_of_downEJs" ~ "num. downstream EJs",
                               parameter == "stop_to_lastEJ" ~ "distance stop-lastEJ",
                               parameter == "NMD_n_sig_tx_perc" ~ "NMD relevance (%)",
                               parameter == "Mech_score" ~ "(d) mech. Z-score",
                               parameter == "utr3_mfe_nt" ~ paste0("-","\u0394","G/nt"),
                               parameter == "GC_utr3" ~ "GC% 3'UTR",
                               parameter == "cds_len" ~ "length CDS",
                               parameter == "nexon" ~ "num. exons",
                               parameter == "GC_utr3_200nt" ~ "GC% 200nt after stop",
                               parameter == "GC_tx" ~ "GC% transcript",
                               parameter == "GC_cds" ~ "GC% CDS",
                               parameter == "L2FC_UPF1_RIP" ~ "(b) log2FC(UPF1-RIP)",
                               parameter == "utr3_len" ~ "length 3'UTR",
                               parameter == "tx_len" ~ "length transcript",
                               parameter == "TrP" ~ "(a) Ribo-Seq",
                               parameter == "PolyCyt" ~ "(a) polyribo./cyto.",
                               parameter == "L2FC_kdeg" ~ "(d) log2FC(kdeg)",
                               parameter == "Copies" ~ "(a) average tx copies/cell",
                               parameter == "Deg" ~ "(a) degradation rate",
                               parameter == "L2FC_hl_siUPF1_Tani" ~ "(b) log2FC(t1/2)-UPF1-KD",
                               parameter == "Proc" ~ "(a) processing rate",
                               parameter == "Syn" ~ "(a) synthesis rate",
                               parameter == "CytNuc" ~ "(a) cyto./nuclear",
                               parameter == "upf1_motifs_top4" ~ "(b) UPF1-top4-bindMotif",
                               parameter == "utr5_len" ~ "length 5'UTR",
                               parameter == "GC_utr5" ~ "GC% 5'UTR",
                               parameter == "upf1_motifs_top1" ~ "(b) UPF1-top1-bindMotif",
                               parameter == "Mech_conclusion" ~ "(d) mech. conclusion",
                               parameter == "L2FC_pUPF1_Kurosaki_new" ~ "(c) log2FC(pUPF1)",
                               parameter == "shadowMax" ~ "(f) Max.",
                               parameter == "shadowMean" ~ "(f) Mean",
                               parameter == "shadowMin" ~ "(f) Min."))

### Save Boruta Data --------------------------------------------------------

save(boruta_fullParam_earlyUP_all_output,
     boruta_redParam_earlyUP_all_output,
     boruta_minParam_earlyUP_all_output,
     boruta_fullParam_earlyUP_no50nt_output,
     boruta_redParam_earlyUP_no50nt_output,
     boruta_minParam_earlyUP_no50nt_output,
     boruta_fullParam_earlyUP_NMD_output,
     boruta_redParam_earlyUP_NMD_output,
     boruta_minParam_earlyUP_NMD_output,
     combined_boruta_analysis,
     file = paste0("Resources/Boruta/",
                   "Boruta_data.rds")
)

