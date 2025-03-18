#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Translation
# Objective: Quantify Ribo-Seq data (0h & 12h IAA on HCT116 N-AID-UPF1) and determine translation efficiency (TE)
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Setup -------------------------------------------------------------------
library(Rsubread)
library(ggpointdensity)
library(DESeq2)

# RNA counts via Salmon->Swish --------------------------------------------------------------

### Import counts from Salmon -----------------------------------------------

mydir="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
condition <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate data frame used for tximeta
coldata <- tibble(files = file.path(mydir, "Salmon", samples$sample, "quant.sf"))

# Supplement with sample IDs as "names"
coldata$names <- samples$sample

# Join with samples to get condition
coldata <- coldata %>% 
  left_join(samples,
            by = c("names" = "sample")) %>% 
  mutate(condition = fct_relevel(as_factor(condition),
                                 "control"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(coldata$files)))

# Use tximeta to import salmon data
all_y <- tximeta::tximeta(coldata,
                 type = "salmon",
                 txOut = TRUE,
                 useHub = FALSE) # reads in counts and inf reps

all_y_gene <- tximeta::summarizeToGene(all_y)

# Scale inferential replicate counts
all_y_gene <- fishpond::scaleInfReps(all_y_gene) # scales counts

# Apply filtering step (default parameters for one condition)
all_y_gene <- fishpond::labelKeep(all_y_gene, x = "condition") # labels features to keep

# Label spike-ins (ERCC and SIRV) as FALSE in keep column -> do not analyze them!
mcols(all_y_gene)$keep <- !str_detect(rownames(assays(all_y_gene)[["counts"]]), 'ERCC') & !str_detect(rownames(assays(all_y_gene)[["counts"]]), 'SIRV')

# Remove spike-in genes from analysis (not just labelled as "keep == FALSE", but really remove)
all_y_gene <- all_y_gene[mcols(all_y_gene)$keep,]

raw_counts <- assay(all_y_gene, "counts")

RNA_raw_counts <- as.data.frame(raw_counts) %>% 
  dplyr::select("191508", "191510", "191512", "191532", "191534", "191536")

coldata_RNA <- read_delim("/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW/experiment.txt", 
                           delim = "\t",
                           escape_double = FALSE, 
                           col_names = c("samples", "condition"),
                           trim_ws = TRUE) %>% 
  column_to_rownames(var = "samples") %>% 
  filter(condition %in% c("UPF1_Nter_0h", "UPF1_Nter_12h")) %>% 
  mutate(condition = as.factor(condition)) %>% 
  mutate(seqtype = "RNA")

# Save data for re-use
save(RNA_raw_counts,
     coldata_RNA,
  file = paste0("Resources/Translation/RNA_raw_counts.rds"))
  
# Load data for re-use - if necessary
load(file = paste0("Resources/Translation/RNA_raw_counts.rds"))

# Ribo-Seq ----------------------------------------------------------------

Reference <- "Resources/GENCODE/gencode.v42.SIRVomeERCCome.annotation.gtf"

BAMs <- c("/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_0h_IAA_1.bam",
          "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_0h_IAA_2.bam",
          "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_0h_IAA_3.bam",
          "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_12h_IAA_1.bam",
          "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_12h_IAA_2.bam",
          "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM/HCT_N_AID_UPF1_12h_IAA_3.bam")

featureCounts_Out <- Rsubread::featureCounts(files = BAMs,
                                   annot.ext = Reference,
                                   isGTFAnnotationFile=TRUE,
                                   GTF.featureType="exon",
                                   GTF.attrType="gene_id",
                                   nthreads=10,
                                   strandSpecific=1)

counts_RiboSeq <- as.data.frame(featureCounts_Out$counts) %>% 
  rename_with(~str_remove(., '.bam')) 

coldata_Ribo <- read_delim("/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/experiment_forCRSA.txt", 
                      delim = "\t",
                      escape_double = FALSE, 
                      col_names = c("samples", "condition"),
                      trim_ws = TRUE) %>% 
  column_to_rownames(var = "samples") %>% 
  mutate(condition = case_when(condition == "control" ~ "UPF1_Nter_0h",
                               condition == "UPF1" ~ "UPF1_Nter_12h")) %>% 
  mutate(condition = as.factor(condition)) %>% 
  mutate(seqtype = "Ribo")

# Save data for re-use
save(counts_RiboSeq,
     coldata_Ribo,
     file = paste0("Resources/Translation/RiboSeq_counts.rds"))

# Load data for re-use - if necessary
load(file = paste0("Resources/Translation/RiboSeq_counts.rds"))

# Combined ----------------------------------------------------------------

# Merge Sample Info
sample_info <- bind_rows(coldata_RNA,coldata_Ribo)

# Merge count data
countData <- RNA_raw_counts %>% 
  mutate_if(is.numeric, as.integer) %>% 
  rownames_to_column(var="gene_id") %>% 
  left_join(counts_RiboSeq %>% 
              rownames_to_column(var="gene_id")) %>% 
  column_to_rownames(var = "gene_id")

countData_tbl <- countData %>% 
  rownames_to_column(var="gene_id")

# Create combined DESeq2 object
ddsMat <- DESeqDataSetFromMatrix(countData = countData,
                                 colData = sample_info,
                                 design = ~condition + seqtype + condition:seqtype)

ddsMat_keep_string <- rowSums(counts(ddsMat) >= 10 ) >= nrow(sample_info)/2
ddsMat <- ddsMat[ddsMat_keep_string,]

# Relevel reference settings
ddsMat$condition <- relevel(ddsMat$condition, ref = "UPF1_Nter_0h")
ddsMat$seqtype <- relevel(ddsMat$seqtype, ref = "RNA")

# Run DESeq2 for Translational Efficiency
ddsMat=DESeq(ddsMat)

resultsNames(ddsMat)

res=results(ddsMat, name="conditionUPF1_Nter_12h.seqtypeRibo")

res_TE_df <- as_tibble(res,
          rownames="gene_id")

# RNA ----------------------------------------------------------------

ddsMat_rna <- DESeqDataSetFromMatrix(countData= RNA_raw_counts %>% 
                                       mutate_if(is.numeric, as.integer),
                                     colData=coldata_RNA,
                                     design=~condition)

ddsMat_rna_keep_string <- rowSums(counts(ddsMat_rna) >= 10 ) >= nrow(coldata_RNA)/2
ddsMat_rna <- ddsMat_rna[ddsMat_rna_keep_string,]

ddsMat_rna=DESeq(ddsMat_rna)
ddsMat_rna$condition <- relevel(ddsMat_rna$condition, ref = "UPF1_Nter_0h")

resultsNames(ddsMat_rna)

res_rna=results(ddsMat_rna, name="condition_UPF1_Nter_12h_vs_UPF1_Nter_0h")

res_rna=lfcShrink(ddsMat_rna,coef="condition_UPF1_Nter_12h_vs_UPF1_Nter_0h",res=res_rna)

# Ribo ----------------------------------------------------------------

ddsMat_ribo <- DESeqDataSetFromMatrix(countData = counts_RiboSeq,
                                      colData = coldata_Ribo, 
                                      design = ~condition)

ddsMat_ribo_keep_string <- rowSums(counts(ddsMat_ribo) >= 10 ) >= nrow(coldata_Ribo)/2
ddsMat_ribo <- ddsMat_ribo[ddsMat_ribo_keep_string,]

ddsMat_ribo <- DESeq(ddsMat_ribo)
ddsMat_ribo$condition <- relevel(ddsMat_ribo$condition, ref = "UPF1_Nter_0h")

res_ribo <- results(ddsMat_ribo, name="condition_UPF1_Nter_12h_vs_UPF1_Nter_0h")
res_ribo <- lfcShrink(ddsMat_ribo, coef="condition_UPF1_Nter_12h_vs_UPF1_Nter_0h",res=res_ribo,type="apeglm")


# Combined TE results -----------------------------------------------------

# Load essential data for these analyses
load("Resources/NMD_relevance/DESeq2_DGE_combined_DGE_cluster_NMD_relevance.rds")
gtf_gencode_df_short_DGE_cluster <- read_csv("Resources/GENCODE/gtf_gencode_df_short_DGE_cluster.csv")

combined_TE_results <- as_tibble(res,
                                 rownames="gene_id") %>% 
  dplyr::select(gene_id,
                log2FoldChange,
                padj) %>% 
  dplyr::rename("L2FC_TE" = "log2FoldChange",
                "padj_TE" = "padj") %>% 
  left_join(as_tibble(res_rna,
                      rownames="gene_id") %>% 
              dplyr::select(gene_id,
                     log2FoldChange,
                     padj)) %>% 
  dplyr::rename("L2FC_RNA" = "log2FoldChange",
                "padj_RNA" = "padj") %>% 
  left_join(as_tibble(res_ribo,
                      rownames="gene_id") %>% 
              dplyr::select(gene_id,
                            log2FoldChange,
                            padj)) %>% 
  dplyr::rename("L2FC_Ribo" = "log2FoldChange",
                "padj_Ribo" = "padj") %>% 
  left_join(gtf_gencode_df_short_DGE_cluster) %>% 
  mutate(class = case_when(padj_Ribo < 0.01 & padj_RNA < 0.0001 & L2FC_RNA*L2FC_Ribo > 0 & L2FC_RNA > 1 ~ "Concordant_up",
                           padj_Ribo < 0.01 & padj_RNA < 0.0001 & L2FC_RNA*L2FC_Ribo > 0 & L2FC_RNA < -1 ~ "Concordant_down",
                           padj_TE < 0.01 & L2FC_TE > 1 ~ "Discordant_up",
                           padj_TE < 0.01 & L2FC_TE < -1 ~ "Discordant_down",
                           TRUE ~ "Not Sig.")) %>% 
  mutate(class = fct_relevel(class,
                             "Concordant_up",
                             "Concordant_down",
                             "Discordant_up",
                             "Discordant_down",
                             "Not Sig.")) %>% 
  left_join(DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% dplyr::select(gene_id, gene_name, gene_type, type, DGE_cluster, NMD_n_sig, NMD_n_sig_perc, UpDown, NMD_bin))

# Save pre-parsed datasources - as csv ---------------------------------------------

combined_TE_results  %>%  write_csv(file.path("Resources/Translation", "Rev_1_combined_TE_results.csv"))