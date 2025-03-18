#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Gene_Level
# Objective: Perform gene-level data preparation and initial analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)
library(DESeq2)              # To work with imported counts
library(tximeta)             # To import salmon counts for Heatmap-Clustering

##
# rRNA quantification -----------------------------------------------------
##

rRNA_QC <- TidyMultiqc::load_multiqc(UPF1_NMDRHT_datasets %>% 
                                       mutate(QC_path = paste0(path,"/QC/QC_rRNA/multiqc_data/multiqc_data.json")) %>% 
                                       pull(QC_path)) %>% 
  separate(metadata.sample_id, c("metadata.sample_id", NA), sep = "_bowtie2_") %>% 
  left_join(UPF1_NMDRHT_datasets_experiments,
            by = c("metadata.sample_id" = "sample")) %>% 
  filter(experimentSet != "HEK293_NMD_Boehm_2018")

# Save pre-parsed datasources - as csv
rRNA_QC  %>%  write_csv(file.path("Resources/QC", "rRNA_QC.csv"))

##
# Load Salmon QC data ---------------------------------------------------------------
##

# Use TidyMultiqc for importing multiple multiqc output files
UPF1_Salmon_QC <- TidyMultiqc::load_multiqc(UPF1_NMDRHT_datasets %>% 
                                              mutate(QC_path = paste0(path,"/QC/QC_salmon/multiqc_data/multiqc_data.json")) %>% 
                                              pull(QC_path)) %>% 
  left_join(UPF1_NMDRHT_datasets_experiments,
            by = c("metadata.sample_id" = "sample")) %>% 
  filter(experimentSet != "HEK293_NMD_Boehm_2018")

# Save pre-parsed datasources - as csv
UPF1_Salmon_QC  %>%  write_csv(file.path("Resources/QC", "UPF1_Salmon_QC.csv"))

##
# PCA gene-level ----------------------------------------------------------
##

### PCA Degradation -----------------------------------------------------------

# Define directory
mydir = "/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")

# Supplement with sample IDs as "names"
names(files) <- samples$sample

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v42.SIRV.ERCC.tsv"))
txi <- tximport::tximport(files, 
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Pre-filter more stringently: require at least half of the samples to have equal or more than 10 counts per gene
# Also remove ERCC and SIRV spike-in genes
keep_string <- rowSums(counts(ddsTxi) >= 10 ) >= nrow(samples)/2 & !str_detect(rownames(counts(ddsTxi)), 'ERCC') & !str_detect(rownames(counts(ddsTxi)), 'SIRV')

ddsTxi_string <- ddsTxi[keep_string,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
dds_string <- DESeq(ddsTxi_string)

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Principal component plot of the samples
plt_deg <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE) %>% 
  mutate(condition = case_when(condition == "control" ~ "control\n0h",
                               condition == "control_48h" ~ "control\n48h",
                               condition == "UPF1_Nter_0h" ~ "N-AID-\nUPF1\n0h",
                               condition == "UPF1_Nter_2h" ~ "N-AID-\nUPF1\n2h",
                               condition == "UPF1_Nter_4h" ~ "N-AID-\nUPF1\n4h",
                               condition == "UPF1_Nter_8h" ~ "N-AID-\nUPF1\n8h",
                               condition == "UPF1_Nter_12h" ~ "N-AID-\nUPF1\n12h",
                               condition == "UPF1_Nter_24h" ~ "N-AID-\nUPF1\n24h",
                               condition == "UPF1_Nter_48h" ~ "N-AID-\nUPF1\n48h",))

percentVar_deg <- round(100 * attr(plt_deg, "percentVar"))

# save data
save(plt_deg,
     percentVar_deg,
     file = paste0("Resources/QC/PCA_UPF1_degradation.rds"))

### PCA Recovery -----------------------------------------------------------

# Define directory
mydir = "/home/volker/gencode.v42.datasets/2023_UPF1_AID_recovery_DW"

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")

# Supplement with sample IDs as "names"
names(files) <- samples$sample

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v42.SIRV.ERCC.tsv"))
txi <- tximport::tximport(files, 
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Pre-filter more stringently: require at least half of the samples to have equal or more than 10 counts per gene
# Also remove ERCC and SIRV spike-in genes
keep_string <- rowSums(counts(ddsTxi) >= 10 ) >= nrow(samples)/2 & !str_detect(rownames(counts(ddsTxi)), 'ERCC') & !str_detect(rownames(counts(ddsTxi)), 'SIRV')

ddsTxi_string <- ddsTxi[keep_string,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
dds_string <- DESeq(ddsTxi_string)

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Principal component plot of the samples
plt_rec <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE) %>% 
  mutate(condition = case_when(condition == "control" ~ "control\n0h",
                               condition == "UPF1_Nter_24h_R0h" ~ "N-AID-\nUPF1\n24h-R0h",
                               condition == "UPF1_Nter_24h_R2h" ~ "N-AID-\nUPF1\n24h-R2h",
                               condition == "UPF1_Nter_24h_R4h" ~ "N-AID-\nUPF1\n24h-R4h",
                               condition == "UPF1_Nter_24h_R8h" ~ "N-AID-\nUPF1\n24h-R8h",
                               condition == "UPF1_Nter_24h_R12h" ~ "N-AID-\nUPF1\n24h-R12h",
                               condition == "UPF1_Nter_24h_R24h" ~ "N-AID-\nUPF1\n24h-R24h",
                               condition == "UPF1_Nter_24h_R48h" ~ "N-AID-\nUPF1\n24h-R48h"))

percentVar_rec <- round(100 * attr(plt_rec, "percentVar"))

# save data
save(plt_rec,
     percentVar_rec,
     file = paste0("Resources/QC/PCA_UPF1_recovery.rds"))

### PCA FKBP-HCT -----------------------------------------------------------

# Define directory
mydir = "/home/volker/gencode.v42.datasets/2023_UPF1_FKBP_HCT_DW"

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")

# Supplement with sample IDs as "names"
names(files) <- samples$sample

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v42.SIRV.ERCC.tsv"))
txi <- tximport::tximport(files, 
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Pre-filter more stringently: require at least half of the samples to have equal or more than 10 counts per gene
# Also remove ERCC and SIRV spike-in genes
keep_string <- rowSums(counts(ddsTxi) >= 10 ) >= nrow(samples)/2 & !str_detect(rownames(counts(ddsTxi)), 'ERCC') & !str_detect(rownames(counts(ddsTxi)), 'SIRV')

ddsTxi_string <- ddsTxi[keep_string,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
dds_string <- DESeq(ddsTxi_string)

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Principal component plot of the samples
plt_FKBP_HCT <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE) %>% 
  mutate(condition = case_when(condition == "control" ~ "control\n0h",
                               condition == "UPF1_FKBP_HCT_0h" ~ "HCT116\nFKBP-UPF1\n0h",
                               condition == "UPF1_FKBP_HCT_12h" ~ "HCT116\nFKBP-UPF1\n12h"))

percentVar_FKBP_HCT <- round(100 * attr(plt_FKBP_HCT, "percentVar"))

# save data
save(plt_FKBP_HCT,
     percentVar_FKBP_HCT,
     file = paste0("Resources/QC/PCA_UPF1_FKBP_HCT.rds"))

### PCA FKBP-HEK -----------------------------------------------------------

# Define directory
mydir = "/home/volker/gencode.v42.datasets/2023_UPF1_FKBP_HEK_DW"

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")

# Supplement with sample IDs as "names"
names(files) <- samples$sample

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v42.SIRV.ERCC.tsv"))
txi <- tximport::tximport(files, 
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Pre-filter more stringently: require at least half of the samples to have equal or more than 10 counts per gene
# Also remove ERCC and SIRV spike-in genes
keep_string <- rowSums(counts(ddsTxi) >= 10 ) >= nrow(samples)/2 & !str_detect(rownames(counts(ddsTxi)), 'ERCC') & !str_detect(rownames(counts(ddsTxi)), 'SIRV')

ddsTxi_string <- ddsTxi[keep_string,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
dds_string <- DESeq(ddsTxi_string)

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Principal component plot of the samples
plt_FKBP_HEK <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE) %>% 
  mutate(condition = case_when(condition == "control" ~ "control\n0h",
                               condition == "UPF1_FKBP_HEK_0h" ~ "HEK293\nFKBP-UPF1\n0h",
                               condition == "UPF1_FKBP_HEK_12h" ~ "HEK293\nFKBP-UPF1\n12h"))

percentVar_FKBP_HEK <- round(100 * attr(plt_FKBP_HEK, "percentVar"))

# save data
save(plt_FKBP_HEK,
     percentVar_FKBP_HEK,
     file = paste0("Resources/QC/PCA_UPF1_FKBP_HEK.rds"))

##
# Gene Cluster ----------------------------------------------------------
##

# Get List of significantly upregulated genes from DESeq2 DGE analysis - means: retain any gene_id which is at least once (in >= 1 condition) significantly upregulated
Rev_1_F2_B_DESeq2_DGE_sig_upregulated_list <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  group_by(gene_id) %>% 
  filter(padj < 0.0001) %>% 
  filter(log2FoldChange > 1) %>% 
  distinct(gene_id) %>% 
  pull(gene_id)

# Downregulated as well
Rev_1_F2_B_DESeq2_DGE_sig_downregulated_list <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  group_by(gene_id) %>% 
  filter(padj < 0.0001) %>% 
  filter(log2FoldChange < -1) %>% 
  distinct(gene_id) %>% 
  pull(gene_id)

### Clustering gene-level ---------------------------------------------------

Rev_1_F2_B_mydir1="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"
Rev_1_F2_B_mydir2="/home/volker/gencode.v42.datasets/2023_UPF1_AID_recovery_DW"

# Get samples and check if all files are present
Rev_1_F2_B_samples1 <- read.table(file.path(Rev_1_F2_B_mydir1, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) 

Rev_1_F2_B_samples2 <- read.table(file.path(Rev_1_F2_B_mydir2, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "Recovery_control",
                               TRUE ~ condition))

Rev_1_F2_B_samples <- as_tibble(bind_rows(Rev_1_F2_B_samples1,
                                          Rev_1_F2_B_samples2))

# Get unique conditions
Rev_1_F2_B_condition <- Rev_1_F2_B_samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate data frame used for tximeta
Rev_1_F2_B_coldata <- tibble(files = file.path(Rev_1_F2_B_mydir1, "Salmon", Rev_1_F2_B_samples1$sample, "quant.sf")) %>% 
  bind_rows(tibble(files = file.path(Rev_1_F2_B_mydir2, "Salmon", Rev_1_F2_B_samples2$sample, "quant.sf"))) 

# Supplement with sample IDs as "names"
Rev_1_F2_B_coldata$names <- Rev_1_F2_B_samples$sample

# Join with samples to get condition
Rev_1_F2_B_coldata <- Rev_1_F2_B_coldata %>% 
  left_join(Rev_1_F2_B_samples,
            by = c("names" = "sample")) %>% 
  mutate(condition = fct_relevel(as_factor(condition),
                                 "control"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(Rev_1_F2_B_coldata$files)))

# Use tximport to import salmon data
Rev_1_F2_B_all_y <- tximeta(Rev_1_F2_B_coldata,
                            type = "salmon",
                            txOut = TRUE,
                            useHub = FALSE) # reads in counts and inf reps

Rev_1_F2_B_all_y_gene <- summarizeToGene(Rev_1_F2_B_all_y)

# Check which assays are available
assayNames(Rev_1_F2_B_all_y_gene)

# Label spike-ins (ERCC and SIRV) as FALSE in keep column -> do not analyze them!
mcols(Rev_1_F2_B_all_y_gene)$keep <- !str_detect(rownames(assays(Rev_1_F2_B_all_y_gene)[["counts"]]), 'ERCC') & !str_detect(rownames(assays(Rev_1_F2_B_all_y_gene)[["counts"]]), 'SIRV')

# Remove spike-in genes from analysis (not just labelled as "keep == FALSE", but really remove)
Rev_1_F2_B_all_y_gene <- Rev_1_F2_B_all_y_gene[mcols(Rev_1_F2_B_all_y_gene)$keep,]

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSet(Rev_1_F2_B_all_y_gene,
                       design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Extracting count data log2-transformed values 
# According to https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog
dds_vsd <- vst(ddsTxi, blind=FALSE)

# Which dimensions does the matrix have?
dim(dds_vsd)

### Normalize counts --------------------------------------------------------

# extract the normalized, but not z-scaled counts:
Rev_1_F2_B_gene_cts_up <- assay(dds_vsd)[rownames(assay(dds_vsd)) %in% Rev_1_F2_B_DESeq2_DGE_sig_upregulated_list,]
Rev_1_F2_B_gene_cts_down <- assay(dds_vsd)[rownames(assay(dds_vsd)) %in% Rev_1_F2_B_DESeq2_DGE_sig_downregulated_list,]

# Z-scale the log2-transformed count matrix
Rev_1_F2_B_gene_test_up <- (t(scale(t(Rev_1_F2_B_gene_cts_up))))
Rev_1_F2_B_gene_test_down <- (t(scale(t(Rev_1_F2_B_gene_cts_down))))

###  UP - Hierarchical clustering --------------------------------------------------------
set.seed(321)
Rev_1_F2_B_gene_dist_up <- dist(Rev_1_F2_B_gene_test_up)

Rev_1_F2_B_gene_hclust_up <- hclust(Rev_1_F2_B_gene_dist_up, method = "ward.D2")

# Elbow method
# set.seed(321)
# Rev_1_F2_B_gene_hclust_wss_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_up, hcut, method = "wss") +
#   labs(subtitle = "Hierarchical - Elbow method")
# 
# # Silhouette method
# set.seed(321)
# Rev_1_F2_B_gene_hclust_silhouette_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_up, hcut, method = "silhouette") +
#   labs(subtitle = "Hierarchical - Silhouette method")
# 
# # Gap statistic
# set.seed(321)
# Rev_1_F2_B_gene_hclust_gap_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_up, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
#   labs(subtitle = "Hierarchical - Gap statistic method")

Rev_1_F2_B_gene_hclust_up_gene_cluster <- cutree(Rev_1_F2_B_gene_hclust_up, k = 4) %>% 
  enframe() %>% 
  dplyr::rename(gene_id = name, cluster = value)

head(Rev_1_F2_B_gene_hclust_up_gene_cluster)

table(cutree(Rev_1_F2_B_gene_hclust_up, k = 4))

Rev_1_F2_B_gene_hclust_up_gene_cluster <- Rev_1_F2_B_gene_hclust_up_gene_cluster %>% 
  mutate(label=case_when(cluster == 1 ~ "1:early",
                         cluster == 2 ~ "3:late",
                         cluster == 3 ~ "2:delayed",
                         cluster == 4 ~ "4:inverse")) %>% 
  pull(label, name=gene_id)

###  DOWN - Hierarchical clustering --------------------------------------------------------
set.seed(321)
Rev_1_F2_B_gene_dist_down <- dist(Rev_1_F2_B_gene_test_down)

Rev_1_F2_B_gene_hclust_down <- hclust(Rev_1_F2_B_gene_dist_down, method = "ward.D2")

# Elbow method
# set.seed(321)
# Rev_1_F2_B_gene_down_hclust_wss_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_down, hcut, method = "wss") +
#   labs(subtitle = "Hierarchical - Elbow method")
# 
# # Silhouette method
# set.seed(321)
# Rev_1_F2_B_gene_down_hclust_silhouette_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_down, hcut, method = "silhouette") +
#   labs(subtitle = "Hierarchical - Silhouette method")
# 
# # Gap statistic
# set.seed(321)
# Rev_1_F2_B_gene_down_hclust_gap_plot <- fviz_nbclust(Rev_1_F2_B_gene_test_down, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
#   labs(subtitle = "Hierarchical - Gap statistic method")

Rev_1_F2_B_gene_hclust_down_gene_cluster <- cutree(Rev_1_F2_B_gene_hclust_down, k = 4) %>% 
  enframe() %>% 
  dplyr::rename(gene_id = name, cluster = value)

head(Rev_1_F2_B_gene_hclust_down_gene_cluster)

table(cutree(Rev_1_F2_B_gene_hclust_down, k = 4))

Rev_1_F2_B_gene_hclust_down_gene_cluster <- Rev_1_F2_B_gene_hclust_down_gene_cluster %>% 
  mutate(label=case_when(cluster == 1 ~ "4:inverse",
                         cluster == 2 ~ "1:early",
                         cluster == 3 ~ "3:late",
                         cluster == 4 ~ "2:delayed")) %>% 
  pull(label, name=gene_id)

#### Generate data frames -------------------------------------------------------------------------

# Get gene_id to DGE_cluster combination for up- and downregulated
Rev_1_F2_B_gene_klus_cluster_up_df <- as_tibble(Rev_1_F2_B_gene_hclust_up_gene_cluster, rownames = "gene_id") %>% 
  dplyr::rename("DGE_cluster_up" = "value") %>% 
  mutate(DGE_cluster_up = paste0("up ",DGE_cluster_up))

Rev_1_F2_B_gene_klus_cluster_down_df <- as_tibble(Rev_1_F2_B_gene_hclust_down_gene_cluster, rownames = "gene_id") %>% 
  dplyr::rename("DGE_cluster_down" = "value") %>% 
  mutate(DGE_cluster_down = paste0("down ",DGE_cluster_down))

### Merge with GENCODE ------------------------------------------------------

# Load simplified gencode version 42 annotation separately
gtf_gencode_df_short <- read_csv("Resources/GENCODE/gencode.v42.gtf_df_short.csv")

# List of expressed genes
Rev_1_F2_B_DESeq2_DGE_expressed_gene_ids <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  dplyr::select(gene_id) %>% 
  dplyr::distinct(gene_id)

# Merge annotation with cluster from Gencode Swish DTE analysis
gtf_gencode_df_short_DGE_cluster <- gtf_gencode_df_short %>% 
  dplyr::filter(type=="gene") %>% 
  dplyr::select(gene_id, gene_name, gene_type) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(Rev_1_F2_B_gene_klus_cluster_up_df,
            by = c("gene_id" = "gene_id")) %>% 
  left_join(Rev_1_F2_B_gene_klus_cluster_down_df,
            by = c("gene_id" = "gene_id")) %>% 
  mutate(DGE_cluster = case_when(is.na(DGE_cluster_up) & is.na(DGE_cluster_down) ~ "NA",
                                 !is.na(DGE_cluster_up) & is.na(DGE_cluster_down) ~ DGE_cluster_up,
                                 is.na(DGE_cluster_up) & !is.na(DGE_cluster_down) ~ DGE_cluster_down,
                                 !is.na(DGE_cluster_up) & !is.na(DGE_cluster_down) ~ paste0("complex"))) %>% 
  mutate(DGE_cluster = case_when(DGE_cluster == "NA" & gene_id %in% Rev_1_F2_B_DESeq2_DGE_expressed_gene_ids$gene_id ~ "expressed",
                                 DGE_cluster == "NA" & !gene_id %in% Rev_1_F2_B_DESeq2_DGE_expressed_gene_ids$gene_id ~ "not_expressed",
                                 TRUE ~ DGE_cluster)) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
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

# Which fraction of genes per UpDown?
gtf_gencode_df_short_DGE_cluster  %>% 
  filter(!DGE_cluster %in% c("expressed", "not_expressed")) %>% 
  mutate(UpDown = case_when(str_detect(DGE_cluster, "up") ~ "up",
                            str_detect(DGE_cluster, "down") ~ "down")) %>% 
  dplyr::count(DGE_cluster, UpDown) %>% 
  group_by(UpDown) %>% 
  mutate(n_per = n/sum(n))

### Fix "complex" genes (with up- and down-clusters) ------------------------------------------------------
# Slice by padj (e.g. most significant gene expression changes)

DGE_complex_genes_fix <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  left_join(gtf_gencode_df_short_DGE_cluster) %>% 
  filter(DGE_cluster == "complex") %>% 
  group_by(gene_id) %>% 
  mutate(min_log2FC = min(log2FoldChange),
         max_log2FC = max(log2FoldChange),
         min_padj = min(padj),
         median(log2FoldChange))  %>% 
  slice_min(padj, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(DGE_cluster_fix= case_when(log2FoldChange > 0 ~ DGE_cluster_up,
                                    log2FoldChange < 0 ~ DGE_cluster_down),
         .after=DGE_cluster)

# Fix genes in relevant dataframe
gtf_gencode_df_short_DGE_cluster <- gtf_gencode_df_short_DGE_cluster %>% 
  left_join(DGE_complex_genes_fix %>% dplyr::select(gene_id, gene_name, DGE_cluster, DGE_cluster_fix)) %>% 
  mutate(DGE_cluster = case_when(DGE_cluster == "complex" ~ DGE_cluster_fix,
                                 TRUE ~ DGE_cluster)) %>% 
  dplyr::select(-DGE_cluster_fix) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
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

# Check DGE_cluster count
gtf_gencode_df_short_DGE_cluster %>% 
  dplyr::count(DGE_cluster)

# Save as csv
gtf_gencode_df_short_DGE_cluster %>% 
  write_csv("Resources/GENCODE/gtf_gencode_df_short_DGE_cluster.csv")

# Fix complex genes in cluster data
## Up
Rev_1_F2_B_gene_klus_cluster_up_df_fix <- Rev_1_F2_B_gene_klus_cluster_up_df %>% 
  filter(gene_id %in% (gtf_gencode_df_short_DGE_cluster %>% 
                         filter(DGE_cluster %in% c("up 1:early",
                                                   "up 2:delayed",
                                                   "up 3:late",
                                                   "up 4:inverse")) %>% 
                         pull(gene_id)))

# Counts
Rev_1_F2_B_gene_klus_cluster_up_df_fix %>% 
  dplyr::count()

Rev_1_F2_B_gene_klus_cluster_up_df_fix %>% 
  dplyr::count(DGE_cluster_up)

### in cluster vector
Rev_1_F2_B_gene_hclust_up_gene_cluster_fix <- Rev_1_F2_B_gene_klus_cluster_up_df_fix %>% 
  pull(DGE_cluster_up, name=gene_id)

## Down
Rev_1_F2_B_gene_klus_cluster_down_df_fix <- Rev_1_F2_B_gene_klus_cluster_down_df %>% 
  filter(gene_id %in% (gtf_gencode_df_short_DGE_cluster %>% 
                         filter(DGE_cluster %in% c("down 4:inverse",
                                                   "down 3:late",
                                                   "down 2:delayed",
                                                   "down 1:early")) %>% 
                         pull(gene_id)))

# Counts
Rev_1_F2_B_gene_klus_cluster_down_df_fix %>% 
  dplyr::count()

Rev_1_F2_B_gene_klus_cluster_down_df_fix %>% 
  dplyr::count(DGE_cluster_down)

### in cluster vector
Rev_1_F2_B_gene_hclust_down_gene_cluster_fix <- Rev_1_F2_B_gene_klus_cluster_down_df_fix %>% 
  pull(DGE_cluster_down, name=gene_id)

### Prepare for Heatmap
Rev_1_F2_B_gene_test_up_fix <- subset(Rev_1_F2_B_gene_test_up, rownames(Rev_1_F2_B_gene_test_up) %in% Rev_1_F2_B_gene_klus_cluster_up_df_fix$gene_id)

Rev_1_F2_B_gene_test_down_fix <- subset(Rev_1_F2_B_gene_test_down, rownames(Rev_1_F2_B_gene_test_down) %in% Rev_1_F2_B_gene_klus_cluster_down_df_fix$gene_id)

### Fit the log2-norm counts ------------------------------------------------

#### Up ------------------------------------------------

# Z-scale the log2-transformed count matrix
Rev_1_F2_B_gene_test_up_klus_cond <- as_tibble(Rev_1_F2_B_gene_test_up_fix,
                                               rownames = "gene_id") %>% 
  janitor::clean_names() %>% 
  left_join(Rev_1_F2_B_gene_klus_cluster_up_df_fix) %>% 
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
  dplyr::select(gene_id,
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
                DGE_cluster_up) %>% 
  ungroup() %>% 
  mutate(DGE_cluster_up = (fct_relevel(DGE_cluster_up,
                                       "up 1:early",
                                       "up 2:delayed",
                                       "up 3:late",
                                       "up 4:inverse")))

#### Down ------------------------------------------------

# Z-scale the log2-transformed count matrix
Rev_1_F2_B_gene_test_down_klus_cond <- as_tibble(Rev_1_F2_B_gene_test_down_fix,
                                                 rownames = "gene_id") %>% 
  janitor::clean_names() %>% 
  left_join(Rev_1_F2_B_gene_klus_cluster_down_df_fix) %>% 
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
  dplyr::select(gene_id,
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
                DGE_cluster_down) %>% 
  ungroup() %>% 
  mutate(DGE_cluster_down = (fct_relevel(DGE_cluster_down,
                                         "down 1:early",
                                         "down 2:delayed",
                                         "down 3:late",
                                         "down 4:inverse")))

### Save cluster information -------------------------------------------
save(Rev_1_F2_B_all_y_gene,
     ddsTxi,
     Rev_1_F2_B_coldata,
     Rev_1_F2_B_gene_cts_up, 
     Rev_1_F2_B_gene_cts_down,
     Rev_1_F2_B_gene_test_up_fix,
     Rev_1_F2_B_gene_test_down_fix,
     Rev_1_F2_B_gene_hclust_up,
     Rev_1_F2_B_gene_hclust_down,
     Rev_1_F2_B_gene_klus_cluster_up_df_fix,
     Rev_1_F2_B_gene_klus_cluster_down_df_fix,
     Rev_1_F2_B_gene_test_up_klus_cond,
     Rev_1_F2_B_gene_test_down_klus_cond,
     Rev_1_F2_B_gene_hclust_up_gene_cluster_fix,
     Rev_1_F2_B_gene_hclust_down_gene_cluster_fix,
     gtf_gencode_df_short_DGE_cluster,
     file = paste0("Resources/Cluster/","Rev_1_F2_B_gene_cluster_information.rds"))

# Save gtf_gencode_df_short_DGE_cluster as csv
gtf_gencode_df_short_DGE_cluster %>% write_csv("Resources/GENCODE/gtf_gencode_df_short_DGE_cluster.csv")

#
##
###
# Additional Analyses ---------------------------------------------------------------
###
##
#

#
## UPF1 RNA-Seq count distribution -----------------------------------------------------------------
#

# Objective: check how the RNA-Seq count distribution looks 
# Potential problem: UPF1 depletion might lead to severely distorted distributions 
# Potential consequence: This might be problematic for normalization and general DGE, etc. analyses
# Approach: check raw, DESeq2-normalized and qsmooth-normalized (PMID: 29036413) counts

library(qsmooth)
library(quantro)

# Use tximport to import salmon data
UPF1_y_qsmooth <- tximeta(Rev_1_F2_B_coldata,
                          type = "salmon",
                          txOut = TRUE,
                          useHub = FALSE) # reads in counts and inf reps

UPF1_y_qsmooth_gene <- summarizeToGene(UPF1_y_qsmooth)

# Check which assays are available
assayNames(UPF1_y_qsmooth_gene)

# Generate DESeqDataSet using samples and condition as parameters
UPF1_ddsTxi_qsmooth_gene <- DESeq2::DESeqDataSet(UPF1_y_qsmooth_gene,
                                         design = ~ condition)

# Filter step: keep only genes that have at least 10 counts in at least 3 samples!
keep <- rowSums(counts(UPF1_ddsTxi_qsmooth_gene) >= 10) >= 3
UPF1_ddsTxi_qsmooth_gene_filt <- UPF1_ddsTxi_qsmooth_gene[keep,]

UPF1_coldata <- Rev_1_F2_B_coldata
# 
# Estimate Size Factors for DESeq2 normalization
UPF1_ddsTxi_qsmooth_gene_filt <- estimateSizeFactors(UPF1_ddsTxi_qsmooth_gene_filt)

# Get raw counts
UPF1_counts_df <- counts(UPF1_ddsTxi_qsmooth_gene_filt)

# qsmooth
qs_norm_e1 <- qsmooth(object = UPF1_counts_df, group_factor = Rev_1_F2_B_coldata$condition)
counts_df_qsmooth <- as_tibble(qsmoothData(qs_norm_e1),
                               rownames = "gene")
##
### Pivot longer  -----------------------------------------------------------------
##

# Raw
UPF1_counts_df_long <- as_tibble(UPF1_counts_df,
                                 rownames = "gene") %>% 
  pivot_longer(!gene,
               names_to = "names",
               values_to = "counts") %>% 
  left_join(UPF1_coldata %>% dplyr::select(names, condition)) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "gene") %>% dplyr::select(gene_id, gene_name),
            by = c("gene" = "gene_id"))
# qsmooth
counts_df_qsmooth_long <- counts_df_qsmooth %>% 
  pivot_longer(!gene,
               names_to = "names",
               values_to = "counts") %>% 
  left_join(UPF1_coldata %>% dplyr::select(names, condition)) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "gene") %>% dplyr::select(gene_id, gene_name),
            by = c("gene" = "gene_id"))

# DESeq2-norm
counts_df_DESeq2_long <- as_tibble(counts(UPF1_ddsTxi_qsmooth_gene_filt, normalized=T),
                                   rownames = "gene") %>% 
  pivot_longer(!gene,
               names_to = "names",
               values_to = "counts") %>% 
  left_join(UPF1_coldata %>% dplyr::select(names, condition)) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "gene") %>% dplyr::select(gene_id, gene_name),
            by = c("gene" = "gene_id"))

##
### Plot basic --------------------------------------------------------
##

# ggplot - without ERCC
density_raw_woERCC <- UPF1_counts_df_long %>% 
  filter(!str_detect(gene, c("ERCC"))) %>% 
  ggplot(aes(x=log2(counts+1), group = names, color=condition)) +
  theme_bw() +
  geom_density() +
  scale_color_viridis_d() +
  labs(title = "raw density")

density_qsmooth_woERCC <- counts_df_qsmooth_long %>% 
  filter(!str_detect(gene, "ERCC")) %>% 
  ggplot(aes(x=log2(counts+1), group = names, color=condition)) +
  theme_bw() +
  geom_density() +
  scale_color_viridis_d() +
  labs(title = "qsmooth density")

boxplot_raw_woERCC <- UPF1_counts_df_long %>% 
  filter(!str_detect(gene, "ERCC")) %>% 
  ggplot(aes(y=log2(counts+1), group = names, x = condition, fill=condition)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "raw boxplot")

boxplot_qsmooth_woERCC <- counts_df_qsmooth_long %>% 
  filter(!str_detect(gene, "ERCC")) %>% 
  ggplot(aes(y=log2(counts+1), group = names, x = condition, fill=condition)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "qsmooth boxplot")


final_plot_woERCC <- cowplot::plot_grid(boxplot_raw_woERCC + theme(legend.position="none"), 
                                        density_raw_woERCC + theme(legend.position="none"),
                                        boxplot_qsmooth_woERCC + theme(legend.position="none"), 
                                        density_qsmooth_woERCC + theme(legend.position="none"),
                                        nrow = 2,
                                        labels = c('A', 'B', 'C', 'D'))

ggsave(file.path("Resources/QC/Expression_qsmooth_woERCC.pdf"),
       final_plot_woERCC,
       width = 25,
       height = 25,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# ggplot - ONLY ERCC
density_raw_ERCC <- UPF1_counts_df_long %>% 
  filter(str_detect(gene, "ERCC")) %>% 
  ggplot(aes(x=log2(counts+1), group = names, color=condition)) +
  theme_bw() +
  geom_density() +
  scale_color_viridis_d() +
  labs(title = "raw density")

density_qsmooth_ERCC <- counts_df_qsmooth_long %>% 
  filter(str_detect(gene, "ERCC")) %>% 
  ggplot(aes(x=log2(counts+1), group = names, color=condition)) +
  theme_bw() +
  geom_density() +
  scale_color_viridis_d() +
  labs(title = "qsmooth density")

boxplot_raw_ERCC <- UPF1_counts_df_long %>% 
  filter(str_detect(gene, "ERCC")) %>% 
  ggplot(aes(y=log2(counts+1), group = names, x = condition, fill=condition)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "raw boxplot")

boxplot_qsmooth_ERCC <- counts_df_qsmooth_long %>% 
  filter(str_detect(gene, "ERCC")) %>% 
  ggplot(aes(y=log2(counts+1), group = names, x = condition, fill=condition)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "qsmooth boxplot")


final_plot_ERCC <- cowplot::plot_grid(boxplot_raw_ERCC + theme(legend.position="none"), 
                                      density_raw_ERCC + theme(legend.position="none"),
                                      boxplot_qsmooth_ERCC + theme(legend.position="none"), 
                                      density_qsmooth_ERCC + theme(legend.position="none"),
                                      nrow = 2,
                                      labels = c('A', 'B', 'C', 'D'))

ggsave(file.path("Resources/QC/Expression_qsmooth_onlyERCC.pdf"),
       final_plot_ERCC,
       width = 25,
       height = 25,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##
### Plot raw vs qsmooth vs DESeq2 -------------------------------------------------------------------------
##

density_raw_norm_woERCC <- counts_df_DESeq2_long %>% 
  filter(!str_detect(gene, c("ERCC"))) %>% 
  ggplot(aes(x=log2(counts+1), group = names, color=condition)) +
  theme_bw() +
  geom_density() +
  scale_color_viridis_d() +
  labs(title = "DESeq2 normalized density")

boxplot_raw_norm_woERCC <- counts_df_DESeq2_long %>% 
  filter(!str_detect(gene, c("ERCC"))) %>% 
  ggplot(aes(y=log2(counts+1), group = names, x = condition, fill=condition)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "DESeq2 normalized boxplot")

final_plot_DESeq2_norm <- cowplot::plot_grid(boxplot_raw_woERCC + theme(legend.position="none"), 
                                             density_raw_woERCC + theme(legend.position="none"),
                                             boxplot_qsmooth_woERCC + theme(legend.position="none"), 
                                             density_qsmooth_woERCC + theme(legend.position="none"),
                                             boxplot_raw_norm_woERCC + theme(legend.position="none"), 
                                             density_raw_norm_woERCC + theme(legend.position="none"),
                                             nrow = 3,
                                             labels = c('A', 'B', 'C', 'D', 'E', 'F'))

ggsave(file.path("Resources/QC/Expression_qsmooth_DESeq2.pdf"),
       final_plot_DESeq2_norm,
       width = 25,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Verdict: The distributions also of raw counts are not that different and normal DESeq2-normalization
# seems to get the job done very well -> no need for qsmooth or ERCC-only normalization

##
# ImpulseDE2 ----------------------------------------------------------
##


### Samples---------------------------------------------------

Rev_1_F2_C_mydir1="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"
Rev_1_F2_C_mydir2="/home/volker/gencode.v42.datasets/2023_UPF1_AID_recovery_DW"

# Get samples and check if all files are present
Rev_1_F2_C_samples1 <- read.table(file.path(Rev_1_F2_C_mydir1, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) 

Rev_1_F2_C_samples2 <- read.table(file.path(Rev_1_F2_C_mydir2, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "Recovery_control",
                               TRUE ~ condition))

Rev_1_F2_C_samples <- as_tibble(bind_rows(Rev_1_F2_C_samples1,
                                          Rev_1_F2_C_samples2))

# Define Samples
Rev_1_F2_C_samples_for_ImpulseDE2 <- Rev_1_F2_C_samples %>% 
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
Rev_1_F2_C_samples_for_ImpulseDE2_Combined <- Rev_1_F2_C_samples_for_ImpulseDE2 %>% 
  filter(Condition == "case") %>% 
  filter(!Sample %in% c("191544",
                        "191546",
                        "191548"))

### Count data---------------------------------------------------

# Modify counts - from dbl to int
Rev_1_F2_C_all_y_gene_mean_ImpulseDE2 <- assay(Rev_1_F2_B_all_y_gene)
mode(Rev_1_F2_C_all_y_gene_mean_ImpulseDE2) <- "integer"

# Run ImpulseDE2 - just case - Batch effects sensitive! - check for transient models
Rev_1_F2_C_ImpulseDE2_Combined <- ImpulseDE2::runImpulseDE2(
  matCountData    = Rev_1_F2_C_all_y_gene_mean_ImpulseDE2, 
  dfAnnotation    = Rev_1_F2_C_samples_for_ImpulseDE2_Combined,
  boolCaseCtrl    = FALSE,
  vecConfounders  = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc        = 15 )

# Get results dataframe
Rev_1_F2_C_ImpulseDE2_Combined_df <- Rev_1_F2_C_ImpulseDE2_Combined$dfImpulseDE2Results

# Combine with annotation - define significance and best fit  
Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster <- gtf_gencode_df_short_DGE_cluster %>% 
  left_join(Rev_1_F2_C_ImpulseDE2_Combined_df,
            by=c("gene_id" = "Gene")) %>% 
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
Rev_1_F2_C_ImpulseDE2_Combined_Impulse_Parameter <- do.call(rbind, lapply(
  Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Impulse") %>% pull(gene_id), function(x) {
    vecImpulseParam = ImpulseDE2::get_lsModelFits(obj=Rev_1_F2_C_ImpulseDE2_Combined)[["case"]][[x]]$lsImpulseFit$vecImpulseParam
  }))

rownames(Rev_1_F2_C_ImpulseDE2_Combined_Impulse_Parameter) <- Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Impulse") %>% pull(gene_id)

# Obtain 2nd Sigmoid parameters
Rev_1_F2_C_ImpulseDE2_Combined_Sigmoid_Parameter <- do.call(rbind, lapply(
  Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Sigmoid") %>% pull(gene_id), function(x) {
    vecImpulseParam = ImpulseDE2::get_lsModelFits(obj=Rev_1_F2_C_ImpulseDE2_Combined)[["case"]][[x]]$lsSigmoidFit$vecSigmoidParam
  }))

rownames(Rev_1_F2_C_ImpulseDE2_Combined_Sigmoid_Parameter) <- Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster %>% filter(bestFit == "Sigmoid") %>% pull(gene_id)

# Combine big dataframe with both information - add t2 and h2 (both NA) for sigmoid just for completeness
Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param <- Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster %>% 
  left_join(bind_rows(as_tibble(Rev_1_F2_C_ImpulseDE2_Combined_Impulse_Parameter,
                                rownames="gene_id"),
                      as_tibble(Rev_1_F2_C_ImpulseDE2_Combined_Sigmoid_Parameter,
                                rownames="gene_id") %>% 
                        dplyr::rename("t1" = "t") %>% 
                        mutate(h2=as.numeric(NA),
                               .after=h1) %>% 
                        mutate(t2=as.numeric(NA),
                               .after=t1)
  ))

# How many significant genes
Rev_1_F2_C_ImpulseDE2_Combined_df %>% 
  filter(padj < 0.0001) %>% 
  dplyr::count()


#### Save impulseDE2 information -------------------------------------------
save(Rev_1_F2_C_ImpulseDE2_Combined_df,
     Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param,
     Rev_1_F2_C_ImpulseDE2_Combined,
     file = paste0("Resources/ImpulseDE2/","Rev_1_F2_C_ImpulseDE2.rds"))

##
# NMD relevance - Cluster ----------------------------------------------------------
##

# Determine for each significantly up-/downregulated gene in how many conditions it is significant
# Consider only concordant events (e.g. Up-cluster with up-regulated events | Down-cluster with down-regulated events)
DESeq2_DGE_combined_DGE_cluster_NMD_relevance <- DESeq2_DGE_combined %>% 
  # filter for relevant datasets
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
  # remove non-relevant or redundant conditions
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
  filter(padj < 0.0001) %>%
  filter(abs(log2FoldChange) > 1) %>% 
  mutate(UpDown = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up",
                            log2FoldChange < -1 & padj < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  left_join(gtf_gencode_df_short_DGE_cluster) %>% 
  filter((DGE_cluster %in% c("up 1:early",
                             "up 2:delayed",
                             "up 3:late") & UpDown == "up") | (DGE_cluster %in% c("down 1:early",
                                                                                  "down 2:delayed",
                                                                                  "down 3:late") & UpDown == "down" ) ) %>%  
  group_by(gene_id, gene_name, DGE_cluster, UpDown) %>% 
  summarize(NMD_n_sig=n(),
            median_log2FC = median(log2FoldChange),
            median_padj = median(padj)) %>% 
  ungroup() 

# Fill dataframe up with those genes that were never significant in any of the other NMD-compromised conditions
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_fill <- DESeq2_DGE_combined_DGE_cluster_NMD_relevance %>% 
  bind_rows(gtf_gencode_df_short_DGE_cluster %>% 
              filter(DGE_cluster %in% c("up 1:early",
                                        "down 1:early",
                                        "up 2:delayed",
                                        "down 2:delayed",
                                        "up 3:late",
                                        "down 3:late"
              )) %>% 
              filter(!gene_id %in% DESeq2_DGE_combined_DGE_cluster_NMD_relevance$gene_id) %>% 
              dplyr::select(gene_id, gene_name, DGE_cluster)) %>% 
  replace_na(list(NMD_n_sig = 0, median_log2FC = 0, median_padj = 1)) %>% 
  mutate(UpDown=case_when(is.na(UpDown) & DGE_cluster %in% c("up 1:early",
                                                             "up 2:delayed",
                                                             "up 3:late",
                                                             "up 4:inverse") ~ "up",
                          is.na(UpDown) & DGE_cluster %in% c("down 4:inverse",
                                                             "down 3:late",
                                                             "down 2:delayed",
                                                             "down 1:early") ~ "down",
                          TRUE ~ UpDown)
  )

# Determine for each gene the significant conditions in wide format
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_wide <- DESeq2_DGE_combined %>% 
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
  filter(padj < 0.0001) %>%
  filter(abs(log2FoldChange) > 1) %>% 
  mutate(UpDown = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up",
                            log2FoldChange < -1 & padj < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  dplyr::select(gene_id, gene_name, gene_type, type, condition_2, UpDown) %>% 
  pivot_wider(names_from = condition_2, values_from = UpDown) 

# Join both dataframes
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete <- DESeq2_DGE_combined_DGE_cluster_NMD_relevance_fill %>% 
  left_join(DESeq2_DGE_combined_DGE_cluster_NMD_relevance_wide) %>% 
  relocate(gene_type, type, .after=gene_name) %>% 
  mutate(NMD_n_sig_perc = 100*NMD_n_sig/length(colnames(DESeq2_DGE_combined_DGE_cluster_NMD_relevance_wide %>% dplyr::select(-c(gene_id, gene_name, gene_type, type)))),
         .after=NMD_n_sig)

# Sanity check: Which genes are counted in up- and down? Should be empty
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(gene_id) %>% 
  dplyr::count() %>% 
  filter(n==2)

# How many genes per UpDown per NMD relevance
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  dplyr::count(UpDown, NMD_n_sig_perc) %>% 
  arrange(desc(UpDown), desc(NMD_n_sig_perc))

# Median NMD relevance per DGE_cluster
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  group_by(DGE_cluster) %>% 
  summarize(mean_NMD = mean(NMD_n_sig_perc))

### NMD relevance bins ------------------------------------------------------

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  dplyr::count(UpDown, DGE_cluster) %>% 
  arrange(desc(n))

# Per up-/downregulated -> cut into 4 equal width bins according to NMD relevance
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin <- DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late") & UpDown == "up" | DGE_cluster %in% c("down 1:early",
                                                                               "down 2:delayed",
                                                                               "down 3:late") & UpDown == "down"  ) %>%  
  group_by(UpDown) %>% 
  mutate(NMD_bin = cut_width(NMD_n_sig_perc, 25, boundary = 0), .after=NMD_n_sig_perc) %>% 
  ungroup()

##
# NMD relevance - all genes ----------------------------------------------------------
##

# Obtain NMD relevance for *all* GENCODE genes
gtf_gencode_df_short_DGE_cluster_NMD_relevance <- DESeq2_DGE_combined %>% 
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
  filter(padj < 0.0001) %>%
  filter(abs(log2FoldChange) > 1) %>% 
  mutate(UpDown = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up",
                            log2FoldChange < -1 & padj < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  left_join(gtf_gencode_df_short_DGE_cluster) %>% 
  filter((DGE_cluster %in% c("up 1:early",
                             "up 2:delayed",
                             "up 3:late",
                             "up 4:inverse") & UpDown == "up") | 
           (DGE_cluster %in% c("down 1:early",
                               "down 2:delayed",
                               "down 3:late",
                               "down 4:inverse") & UpDown == "down" ) | 
           (DGE_cluster %in% c("expressed",
                               "not_expressed"))) %>% 
  group_by(gene_id, gene_name, DGE_cluster, UpDown) %>% 
  summarize(NMD_n_sig=n(),
            median_log2FC = median(log2FoldChange),
            median_padj = median(padj)) %>% 
  ungroup() 

# Fill dataframe up with those genes that were never significant in any of the other NMD-compromised conditions
gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill <- gtf_gencode_df_short_DGE_cluster_NMD_relevance %>% 
  bind_rows(gtf_gencode_df_short_DGE_cluster %>% 
              filter(!gene_id %in% gtf_gencode_df_short_DGE_cluster_NMD_relevance$gene_id) %>% 
              dplyr::select(gene_id, gene_name, DGE_cluster)) %>% 
  replace_na(list(NMD_n_sig = 0, median_log2FC = 0, median_padj = 1)) %>% 
  group_by(gene_id, DGE_cluster) %>% 
  slice_max(NMD_n_sig, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  mutate(UpDown=case_when(is.na(UpDown) & DGE_cluster %in% c("up 1:early",
                                                             "up 2:delayed",
                                                             "up 3:late",
                                                             "up 4:inverse") ~ "up",
                          is.na(UpDown) & DGE_cluster %in% c("down 4:inverse",
                                                             "down 3:late",
                                                             "down 2:delayed",
                                                             "down 1:early") ~ "down",
                          is.na(UpDown) ~ "n.s.",
                          TRUE ~ UpDown)
  ) %>%  
  mutate(NMD_n_sig_perc = 100*NMD_n_sig/max(gtf_gencode_df_short_DGE_cluster_NMD_relevance$NMD_n_sig),
         .after=NMD_n_sig)

# Determine NMD bin
gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin <- gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill %>% 
  filter(UpDown != "n.s.") %>% 
  group_by(UpDown) %>% 
  mutate(NMD_bin = cut_width(NMD_n_sig_perc, 25, boundary = 0), .after=NMD_n_sig_perc) %>% 
  ungroup() %>% 
  bind_rows(gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill %>% 
              filter(UpDown == "n.s.") %>% 
              mutate(NMD_bin = "n.s."))

### Save NMD relevance information -------------------------------------------
save(DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete,
     DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin,
     gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin,
     file = paste0("Resources/NMD_relevance/","DESeq2_DGE_combined_DGE_cluster_NMD_relevance.rds"))

# Save as csv
gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin %>% 
  write_csv("Resources/GENCODE/gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin.csv")

##
# Phospho-UPF1 Kurosaki 2014 ----------------------------------------------------------
##

# mydir is the first argument = working directory
mydir_Kurosaki2014="/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki"	

# Get samples and check if all files are present
samples_Kurosaki2014 <- read.table(file.path(mydir_Kurosaki2014, "Samples.txt"), header = TRUE) 

##
### featureCounts -----------------------------------------------------------
##

# reason: salmon performs very bad with those very short reads
# reads were a) adapter-trimmed (cutadapt); b) filtered for >= 15 nt and c) aligned against rRNA, tRNA, miRNA and snoRNA -> unaligned used for further analyses

# Load Reference
Reference <- "Resources/GENCODE/gencode.v42.SIRVomeERCCome.annotation.gtf"

# Define BAM files
# Exclude IgG
BAMs_Kurosaki2014 <- c("/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535652_unaligned.bam",
                       "/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535653_unaligned.bam",
                       "/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535654_unaligned.bam",
                       "/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535655_unaligned.bam",
                       "/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535656_unaligned.bam",
                       "/home/volker/gencode.v42.datasets/2014_p_UPF1_Kurosaki/BAM/SRR1535657_unaligned.bam")

# Perform unstranded featureCounts
featureCounts_Out_Kurosaki2014 <- Rsubread::featureCounts(files = BAMs_Kurosaki2014,
                                                          annot.ext = Reference,
                                                          isGTFAnnotationFile=TRUE,
                                                          GTF.featureType="exon",
                                                          GTF.attrType="gene_id",
                                                          nthreads=10,
                                                          strandSpecific=0)

# Extract counts and fix names
counts_Kurosaki2014 <- as.data.frame(featureCounts_Out_Kurosaki2014$counts) %>% 
  rename_with(~str_remove(., '_unaligned.bam'))

# Generate GENCODE-merged table
counts_Kurosaki2014_tbl <- as_tibble(counts_Kurosaki2014 %>% 
                                       rownames_to_column(var="gene_id")) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "gene") %>% dplyr::select("gene_name", "gene_id"))

# Deisng DESeq2DataSet
ddsMat_Kurosaki2014 <- DESeqDataSetFromMatrix(countData = counts_Kurosaki2014,
                                              colData = samples_Kurosaki2014,
                                              design = ~ condition)

# Set control condition as reference
ddsMat_Kurosaki2014$condition <- relevel(ddsMat_Kurosaki2014$condition, ref = "control")

# Pre-filter 
# Also remove ERCC and SIRV spike-in genes
keep_string_Kurosaki2014 <- !str_detect(rownames(counts(ddsMat_Kurosaki2014)), 'ERCC') & !str_detect(rownames(counts(ddsMat_Kurosaki2014)), 'SIRV')

ddsMat_Kurosaki2014 <- ddsMat_Kurosaki2014[keep_string_Kurosaki2014,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
ddsMat_Kurosaki2014_analysis <- DESeq(ddsMat_Kurosaki2014)

##
### p_UPF1 versus control  ------------------------------------------------------------------
##
res_Kurosaki2014 <- results(ddsMat_Kurosaki2014_analysis,
                            # alpha=0.0001,
                            contrast=c("condition","p_UPF1","control"))	

# Convert results file to tibble
res_Kurosaki2014_tbl <- as_tibble(res_Kurosaki2014,
                                  rownames = "gene_id") %>% 
  mutate(condition_2 = "p_UPF1")

# Get additional IDs from Gencode annotation
res_Kurosaki2014_tbl_final <- left_join(res_Kurosaki2014_tbl,
                                        gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin)

##
### IgG versus control  ------------------------------------------------------------------
##

res_Kurosaki2014_IgG <- results(ddsMat_Kurosaki2014_analysis,
                                # alpha=0.0001,
                                contrast=c("condition","IgG","control"))	

# Convert results file to tibble
res_Kurosaki2014_IgG_tbl <- as_tibble(res_Kurosaki2014_IgG,
                                      rownames = "gene_id") %>% 
  mutate(condition_2 = "IgG")

# Get additional IDs from Gencode annotation
res_Kurosaki2014_IgG_tbl_final <- left_join(res_Kurosaki2014_IgG_tbl,
                                            gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin)

# Merge 
Kurosaki2014_tbl <- res_Kurosaki2014_tbl_final %>% 
  bind_rows(res_Kurosaki2014_IgG_tbl_final)

### Save Kurosaki2014_phospho-UPF1 information -------------------------------------------
save(featureCounts_Out_Kurosaki2014,
     ddsMat_Kurosaki2014_analysis,
     Kurosaki2014_tbl,
     file = paste0("Resources/External/Kurosaki2014_PMID_25184677/","Kurosaki2014_phospho_UPF1.rds"))

Kurosaki2014_tbl %>% write_csv("Resources/External/Kurosaki2014_PMID_25184677/Kurosaki2014_tbl.csv")

##
# Gene-MainTable -------------------------------------------------------------
##

GENCODE_v42_MainTable <- gtf_gencode_df_short %>% 
  filter(type=="gene") %>% 
  dplyr::select(-c(type)) %>% 
  dplyr::select(-c(transcript_id, transcript_type, transcript_name, transcript_support_level, ensembl_canonical, ensembl_basic, appris_principal)) %>% 
  left_join(gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin %>% 
              dplyr::rename("L2FC_median_NMD" = "median_log2FC",
                            "padj_median_NMD" = "median_padj")) %>% 
  left_join(Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
              dplyr::select(gene_id,padj,sigImpulseDE2, bestFit, beta, h0, h1, h2, t1, t2) %>% 
              dplyr::rename("bestFit_ImpulseDE2" = "bestFit",
                            "padj_ImpulseDE2" = "padj")) %>% 
  left_join(combined_TE_results %>% 
              dplyr::rename("L2FC_RNA_RiboSeq" = "L2FC_RNA",
                            "padj_RNA_RiboSeq" = "padj_RNA",
                            "class_RiboSeq" = "class") %>% 
              relocate(L2FC_TE, 
                       padj_TE, .after="padj_Ribo")) %>% 
  left_join(Mechs_UPF1_combined %>% 
              filter(condition == "Nter_12h") %>% 
              dplyr::select(-c(condition, type, transcript_id, transcript_type, transcript_name, transcript_support_level, ensembl_canonical, ensembl_basic, appris_principal,
                               bakR_se, bakR_pval, bakR_score, DE_score, DE_se, DE_pval, mech_pval, meta_pval, ksyn_score, ksyn_pval, f_deg)) %>% 
              dplyr::rename("gene_id" = "XF",
                            "padj_kdeg" = "bakR_padj",
                            "padj_RNA" = "DE_padj",
                            "stat_mech" = "mech_stat",
                            "padj_mech" = "mech_padj",
                            "padj_meta" = "meta_padj",
                            "padj_ksyn" = "ksyn_padj")) %>% 
  left_join(Kurosaki2014_tbl %>% 
              filter(condition_2 == "p_UPF1") %>% 
              dplyr::select(gene_id, log2FoldChange, padj) %>% 
              dplyr::rename("L2FC_p_UPF1" = "log2FoldChange",
                            "padj_p_UPF1" = "padj")) %>% 
  # dplyr::select(-c(NMD_bin)) %>% 
  relocate(type, .after = "gene_type") %>% 
  relocate(DGE_cluster_up, DGE_cluster_down, .after = "DGE_cluster") 

# Export NMDRHT global table
GENCODE_v42_MainTable %>% 
  write_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

# Read data - if necessary
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")
