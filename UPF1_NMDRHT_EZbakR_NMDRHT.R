#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_EZbakR_NMDRHT
# Objective: Differential kinetic analysis on the transcript isoform level via EZbakR on NMDRHT annotation for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(EZbakR)
library(tidyverse)

##
# Load data sources --------------------------------------------------------
##

# Load cB file from bam2bakR
UPF1_NMDRHT_EZbakR_cB_csv <- read_csv("Resources/EZbakR_NMDRHT/cB.csv.gz")

# How does it look?
head(UPF1_NMDRHT_EZbakR_cB_csv)

# Which samples do we have?
print(unique(UPF1_NMDRHT_EZbakR_cB_csv$sample))

# Load Metadata file
UPF1_NMDRHT_metadf <- read_delim("Resources/EZbakR_NMDRHT/metadf", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = TRUE, trim_ws = TRUE)

# Assemble EZbakR dataset -------------------------------------------------------
UPF1_NMDRHT_EZbakR <- EZbakRData(UPF1_NMDRHT_EZbakR_cB_csv, UPF1_NMDRHT_metadf)

# QC
UPF1_NMDRHT_EZbakR_QC <- EZQC(UPF1_NMDRHT_EZbakR)

##
# EZbakR analysis ---------------------------------------------------------
##

# Comment: analysis along tutorial: https://isaacvock.github.io/Isoform_Tutorial_Docs/EZbakR.html

## Estimate TEC fraction new  ---------------------------------------------------------
UPF1_NMDRHT_EZbakR_TEC <- EstimateFractions(UPF1_NMDRHT_EZbakR,
                                         pold_from_nolabel = TRUE,
                                         features = c("XF", "TEC"),
                                         filter_condition = `|`)

## Estimate transcript isoform fraction new ---------------------------------------------------------

### Load isoform quantification information
UPF1_NMDRHT_EZbakR_file_names <- list.files(path = "Resources/EZbakR_NMDRHT/rsem",
                         pattern = "isoform",
                         full.names = TRUE)

names(UPF1_NMDRHT_EZbakR_file_names) <- UPF1_NMDRHT_metadf$sample

### Deconvolve isoform fraction news

UPF1_NMDRHT_EZbakR_TEC <- ImportIsoformQuant(UPF1_NMDRHT_EZbakR_TEC,
                            files = UPF1_NMDRHT_EZbakR_file_names,
                            quant_tool = "rsem")

# Check how many isoforms are quantified
UPF1_NMDRHT_EZbakR_TEC$readcounts$isoform_quant_rsem %>% 
  left_join(UPF1_NMDRHT_metadf %>% dplyr::select(sample,group)) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  ggplot(aes(x=log2(TPM),
             fill=NMD_50nt_rule)) +
  geom_density(alpha=0.5) +
  geom_vline(xintercept = log2(0.2)) +
  scale_fill_manual(values=c("FALSE" = "#214D65",
                             "TRUE"="#E5BF86")) +
  facet_wrap(~group)

# EstimateIsoformFractions
UPF1_NMDRHT_EZbakR_TEC <- EstimateIsoformFractions(UPF1_NMDRHT_EZbakR_TEC,
                                                   TPM_min = 0.2)

## Estimate, average, and compare degradation rate constants ---------------------------------------------------------

UPF1_NMDRHT_EZbakR_TEC <- EstimateKinetics(UPF1_NMDRHT_EZbakR_TEC, features = "transcript_id",
                          exactMatch = FALSE)

### QC plots  ---------------------------------------------------------
# Check how many isoforms per 50-nt rule are quantified - per condition
UPF1_NMDRHT_EZbakR_TEC$kinetics$XF_transcriptid %>% 
  filter(!is.na(log_kdeg)) %>% 
  left_join(UPF1_NMDRHT_metadf %>% dplyr::select(sample,group)) %>% 
  distinct(transcript_id, group) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  dplyr::count(group, NMD_50nt_rule) %>% 
  ggplot(aes(y=group,
             x=n,
             fill=NMD_50nt_rule)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(colour = 'darkgray', linewidth = 0.1),
        panel.grid.major.y = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5), 
        legend.key.size = unit(0.15, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_col(position=position_dodge(width=0.8),
           width=0.75,
           color="black",
           linewidth=0.1) +
  geom_text(aes(label=case_when(NMD_50nt_rule == FALSE ~ paste(n))),
            hjust = 1.25,
            size = 4*0.36,
            color="white",
            position =position_dodge(width=0.8)
  ) +
  geom_text(aes(label=case_when(NMD_50nt_rule == TRUE ~ paste(n))),
            hjust = 1.25,
            size = 4*0.36,
            color="black",
            position =position_dodge(width=0.8)
  ) +
  scale_fill_manual(values=c("FALSE" = "#214D65",
                             "TRUE"="#E5BF86")) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

# Average rate constants
UPF1_NMDRHT_EZbakR_TEC <- AverageAndRegularize(UPF1_NMDRHT_EZbakR_TEC, features = "transcript_id",
                                            formula_mean = ~ group,
                              exactMatch = FALSE)

# Save analyzed EZbakR to file
save(UPF1_NMDRHT_EZbakR_TEC,
     file="Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_TPM02.rds")

# QC
# UPF1_NMDRHT_EZbakR_TEC_QC <- EZQC(UPF1_NMDRHT_EZbakR_TEC)

# load("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC.rds")

### Compare Ctrl-12h vs Ctrl-0h ---------------------------------------

UPF1_NMDRHT_EZbakR_TEC_Ctrl_12_vs_Ctrl_0 <- CompareParameters(UPF1_NMDRHT_EZbakR_TEC,
                                                              features = "transcript_id",
                                                              design_factor = "group",
                                                              reference = "Ctrl_0h_IAA",
                                                              experimental = "Ctrl_12h_IAA",
                                                              exactMatch = FALSE)

### Compare Ctrl-24h vs Ctrl-0h ---------------------------------------

UPF1_NMDRHT_EZbakR_TEC_Ctrl_24_vs_Ctrl_0 <- CompareParameters(UPF1_NMDRHT_EZbakR_TEC,
                                                              features = "transcript_id",
                                                              design_factor = "group",
                                                              reference = "Ctrl_0h_IAA",
                                                              experimental = "Ctrl_24h_IAA",
                                                              exactMatch = FALSE)

### Compare N-AID-UPF1-0h vs Ctrl-0h ---------------------------------------

UPF1_NMDRHT_EZbakR_TEC_Nter_0_vs_Ctrl_0 <- CompareParameters(UPF1_NMDRHT_EZbakR_TEC,
                                                              features = "transcript_id",
                                                              design_factor = "group",
                                                              reference = "Ctrl_0h_IAA",
                                                              experimental = "Nter_0h_IAA",
                                                              exactMatch = FALSE)

### Compare N-AID-UPF1-12h vs Ctrl-0h ---------------------------------------

UPF1_NMDRHT_EZbakR_TEC_Nter_12_vs_Ctrl_0 <- CompareParameters(UPF1_NMDRHT_EZbakR_TEC,
                           features = "transcript_id",
                           design_factor = "group",
                           reference = "Ctrl_0h_IAA",
                           experimental = "Nter_12h_IAA",
                           exactMatch = FALSE)

### Compare N-AID-UPF1-24h vs Ctrl-0h ---------------------------------------

UPF1_NMDRHT_EZbakR_TEC_Nter_24_vs_Ctrl_0 <- CompareParameters(UPF1_NMDRHT_EZbakR_TEC,
                                                              features = "transcript_id",
                                                              design_factor = "group",
                                                              reference = "Ctrl_0h_IAA",
                                                              experimental = "Nter_24h_IAA",
                                                              exactMatch = FALSE)

## Assess results ---------------------------------------

# Load NMDRHT main table - transcript-level information
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

# Load GENCODE main table - gene-level information
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

# Load edgeR DTE *NMDRHT*-annotation-based data
edgeR_DTE_NMDRHT_combined <- read_csv("Resources/edgeR_DTE_NMDRHT_combined.csv")

# Load Mechs_UPF1_combined bakR-gene-level results
Mechs_UPF1_combined <- read_csv(file=file.path("Resources/bakR/Mechs_UPF1_combined.csv"))

### Assess Ctrl-12h vs Ctrl-0h ---------------------------------------
comparison_TEC_Ctrl_12_vs_Ctrl_0 <- EZget(UPF1_NMDRHT_EZbakR_TEC_Ctrl_12_vs_Ctrl_0,
                    type = "comparisons") %>% 
  dplyr::rename("gene_id" = "XF",
                "L2FC_kdeg_tx" = "difference",
                "EZbakR_tx_uncertainty" = "uncertainty",
                "EZbakR_tx_stat" = "stat",
                "EZbakR_tx_pval" = "pval",
                "padj_kdeg_tx" = "padj",
                "EZbakR_tx_avg_coverage" = "avg_coverage") %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(Mechs_UPF1_combined %>%
              dplyr::select(XF, L2FC_kdeg, condition) %>% 
              dplyr::rename("gene_id" = "XF") %>% 
              filter(condition == "Ctrl_12h")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet == "HCT116_UPF1_AID_4SU_this_Study" & condition_2 == "control_12h") %>% 
              dplyr::select(transcript_id, logFC)) %>% 
  mutate(condition="control_12h")

# How many significant transcripts
comparison_TEC_Ctrl_12_vs_Ctrl_0 %>% 
  filter(padj_kdeg_tx < 0.01) %>% 
  dplyr::count()

### Assess Ctrl-24h vs Ctrl-0h ---------------------------------------
comparison_TEC_Ctrl_24_vs_Ctrl_0 <- EZget(UPF1_NMDRHT_EZbakR_TEC_Ctrl_24_vs_Ctrl_0,
                                          type = "comparisons") %>% 
  dplyr::rename("gene_id" = "XF",
                "L2FC_kdeg_tx" = "difference",
                "EZbakR_tx_uncertainty" = "uncertainty",
                "EZbakR_tx_stat" = "stat",
                "EZbakR_tx_pval" = "pval",
                "padj_kdeg_tx" = "padj",
                "EZbakR_tx_avg_coverage" = "avg_coverage") %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(Mechs_UPF1_combined %>%
              dplyr::select(XF, L2FC_kdeg, condition) %>% 
              dplyr::rename("gene_id" = "XF") %>% 
              filter(condition == "Ctrl_24h")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet == "HCT116_UPF1_AID_4SU_this_Study" & condition_2 == "control_24h") %>% 
              dplyr::select(transcript_id, logFC)) %>% 
  mutate(condition="control_24h")

# How many significant transcripts
comparison_TEC_Ctrl_24_vs_Ctrl_0 %>% 
  filter(padj_kdeg_tx < 0.01) %>% 
  dplyr::count()

### Assess N-AID-UPF1-0h vs Ctrl-0h ---------------------------------------
comparison_TEC_Nter_0_vs_Ctrl_0 <- EZget(UPF1_NMDRHT_EZbakR_TEC_Nter_0_vs_Ctrl_0,
                                          type = "comparisons") %>% 
  dplyr::rename("gene_id" = "XF",
                "L2FC_kdeg_tx" = "difference",
                "EZbakR_tx_uncertainty" = "uncertainty",
                "EZbakR_tx_stat" = "stat",
                "EZbakR_tx_pval" = "pval",
                "padj_kdeg_tx" = "padj",
                "EZbakR_tx_avg_coverage" = "avg_coverage") %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(Mechs_UPF1_combined %>%
              dplyr::select(XF, L2FC_kdeg, condition) %>% 
              dplyr::rename("gene_id" = "XF") %>% 
              filter(condition == "Nter_0h")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet == "HCT116_UPF1_AID_4SU_this_Study" & condition_2 == "UPF1_0h") %>% 
              dplyr::select(transcript_id, logFC)) %>% 
  mutate(condition="UPF1_0h")

# How many significant transcripts
comparison_TEC_Nter_0_vs_Ctrl_0 %>% 
  filter(padj_kdeg_tx < 0.01) %>% 
  dplyr::count()

### Assess N-AID-UPF1-12h vs Ctrl-0h ---------------------------------------
comparison_TEC_Nter_12_vs_Ctrl_0 <- EZget(UPF1_NMDRHT_EZbakR_TEC_Nter_12_vs_Ctrl_0,
                                         type = "comparisons") %>% 
  dplyr::rename("gene_id" = "XF",
                "L2FC_kdeg_tx" = "difference",
                "EZbakR_tx_uncertainty" = "uncertainty",
                "EZbakR_tx_stat" = "stat",
                "EZbakR_tx_pval" = "pval",
                "padj_kdeg_tx" = "padj",
                "EZbakR_tx_avg_coverage" = "avg_coverage") %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(Mechs_UPF1_combined %>%
              dplyr::select(XF, L2FC_kdeg, condition) %>% 
              dplyr::rename("gene_id" = "XF") %>% 
              filter(condition == "Nter_12h")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet == "HCT116_UPF1_AID_4SU_this_Study" & condition_2 == "UPF1_12h") %>% 
              dplyr::select(transcript_id, logFC)) %>% 
  mutate(condition="UPF1_12h")

# How many significant transcripts
comparison_TEC_Nter_12_vs_Ctrl_0 %>% 
  filter(padj_kdeg_tx < 0.01) %>% 
  dplyr::count()

### Assess N-AID-UPF1-24h vs Ctrl-0h ---------------------------------------
comparison_TEC_Nter_24_vs_Ctrl_0 <- EZget(UPF1_NMDRHT_EZbakR_TEC_Nter_24_vs_Ctrl_0,
                                          type = "comparisons") %>% 
  dplyr::rename("gene_id" = "XF",
                "L2FC_kdeg_tx" = "difference",
                "EZbakR_tx_uncertainty" = "uncertainty",
                "EZbakR_tx_stat" = "stat",
                "EZbakR_tx_pval" = "pval",
                "padj_kdeg_tx" = "padj",
                "EZbakR_tx_avg_coverage" = "avg_coverage") %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(Mechs_UPF1_combined %>%
              dplyr::select(XF, L2FC_kdeg, condition) %>% 
              dplyr::rename("gene_id" = "XF") %>% 
              filter(condition == "Nter_24h")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet == "HCT116_UPF1_AID_4SU_this_Study" & condition_2 == "UPF1_24h") %>% 
              dplyr::select(transcript_id, logFC)) %>% 
  mutate(condition="UPF1_24h")

# How many significant transcripts
comparison_TEC_Nter_24_vs_Ctrl_0 %>% 
  filter(padj_kdeg_tx < 0.01) %>% 
  dplyr::count()

### Combine -----------------------------------------------------------------
UPF1_NMDRHT_EZbakR_TEC_combined <- rbind(comparison_TEC_Ctrl_12_vs_Ctrl_0,
                                 comparison_TEC_Ctrl_24_vs_Ctrl_0,
                                 comparison_TEC_Nter_0_vs_Ctrl_0,
                                 comparison_TEC_Nter_12_vs_Ctrl_0,
                                 comparison_TEC_Nter_24_vs_Ctrl_0) %>% 
  mutate(kdeg_tx_conclusion = ifelse(padj_kdeg_tx < 0.01, case_when(L2FC_kdeg_tx < -1 ~ "Stabilized",
                                                                    L2FC_kdeg_tx > 1 ~ "Destabilized",
                                                            TRUE ~ "Not Sig."), "Not Sig.")) %>% 
  relocate(kdeg_tx_conclusion, condition, logFC, L2FC_kdeg, .after="EZbakR_tx_avg_coverage")

# Stats kdeg_tx_conclusion
UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  dplyr::count(condition, kdeg_tx_conclusion)

# save data as csv
UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  write_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))

# Read file again if needed
# UPF1_NMDRHT_EZbakR_TEC_combined <- read_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))