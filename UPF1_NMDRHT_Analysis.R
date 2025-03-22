#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Analysis
# Objective: Top-level analysis script for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: it is not recommended to run the whole script from top to bottom, but rather perform the analyses in chunks
# Comment: some data sources are loaded multiple times in different parts of the sub-analyses, hopefully allowing re-running also parts of the analysis

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)           # Generally required
library(patchwork)           # Allows easy plot stitching
library(extrafont)           # For using Arial in some plots
library(ggh4x)               # Allows fixing row/column height
library(viridis)             # Viridis color palette
library(ComplexHeatmap)      # Produce Heatmaps
library(circlize)            # Color for Heatmaps

##
# Setup ----------------------------------------------------------
##

# Set working directory
# Has to be changed when re-analyzed on a different machine!
setwd(file.path("/home/volker/2025_UPF1_NMDRHT"))

# Set column width for plotting
cw1 = 5.5
cw2 = 12
cw3 = 18

# Define Arial as standard for ComplexHeatmap
grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = "Arial")))

# Import custom functions
source("UPF1_NMDRHT_Functions.R")

#
##
# Import data sources and metadata --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run NMDRHT annotation generation

# source("UPF1_NMDRHT_Data_sources.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

## Import processed Metadata ----------------------------------------------------------------

# Load datasets for UPF1 NMDRHT study
UPF1_NMDRHT_datasets <- read_csv("Resources/Metadata/UPF1_NMDRHT_datasets.csv")

# Load individual experiments for UPF1 NMDRHT study
UPF1_NMDRHT_datasets_experiments <- read_csv("Resources/Metadata/UPF1_NMDRHT_datasets_experiment.csv")

## Import processed Analyses ----------------------------------------------------------------

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")

# Load Swish DGE *GENCODE*-annotation-based data
Swish_DGE_combined <- read_csv("Resources/Swish_DGE_combined.csv")

# Load edgeR DTE *GENCODE*-annotation-based data
edgeR_DTE_combined <- read_csv("Resources/edgeR_DTE_combined.csv")

# Load edgeR DTE *NMDRHT*-annotation-based data
edgeR_DTE_NMDRHT_combined <- read_csv("Resources/edgeR_DTE_NMDRHT_combined.csv")

# Load IsoformSwitchAnalyzeR (ISAR) DTU *GENCODE*-annotation-based data
# Note: only reduced datasets available!
ISAR_DTU_combined <- read_csv("Resources/ISAR_DTU_combined.csv")

# Load IsoformSwitchAnalyzeR (ISAR) DTU *NMDRHT*-annotation-based data
# Note: only reduced datasets available!
ISAR_DTU_NMDRHT_combined <- read_csv("Resources/ISAR_DTU_NMDRHT_combined.csv")

#
##
# Generate NMDRHT annotation --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run NMDRHT annotation generation

# source("UPF1_NMDRHT_Annotation.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load simplified gencode version 42 annotation separately
gtf_gencode_df_short <- read_csv("Resources/GENCODE/gencode.v42.gtf_df_short.csv")

# Load shortend NMDRHT information
NMDRHT.v1.1_tbl <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_short.csv")

#
##
# Supplement NMDRHT annotation with Ribo-Seq data --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run NMDRHT annotation supplementation with Ribo-Seq data
# source("UPF1_NMDRHT_ORFs.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load partially analyzed ORFquant and Ribo-TISH data
HCT_N_AID_UPF1_ORFquant <- read_csv("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant.csv")
HCT_N_AID_UPF1_RiboTISH_combined <- read_csv("Resources/RiboTISH/HCT_N_AID_UPF1_RiboTISH.csv")
HCT_N_AID_UPF1_ORFquant_RiboTish <- read_csv("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish.csv")

# Load updated NMDRHT information (contains Ribo-Seq ORFs and transcript properties)
NMDRHT.v1.2_tbl_length_GC_mfe <- as_tibble(read_csv(file.path("Resources/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe.csv")))

#
##
# Gene-level data preparation --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to perform gene-level data preparatgion
# source("UPF1_NMDRHT_Gene_Level.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load rRNA mapping QC multiQC output
rRNA_QC  <-  read_csv(file.path("Resources/QC", "rRNA_QC.csv"))

# Load rRNA mapping QC multiQC output
UPF1_Salmon_QC  <-  read_csv(file.path("Resources/QC", "UPF1_Salmon_QC.csv"))

# Read data - if necessary
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

#
##
# Metabolic mRNA labeling data (4SU) --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run metabolic labeling analysis via bakR
# source("UPF1_NMDRHT_bakR.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load 4SU kdeg data - Hybrid model from bakR (https://simonlabcode.github.io/bakR/index.html)
Mechs_UPF1_combined <- 
  read_csv(file=file.path("Resources/bakR/Mechs_UPF1_combined.csv"))

#
##
# Translation Efficiency --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run translation analysis
# source("UPF1_NMDRHT_Translation.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

#
##
# Transcript-level data preparation --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to perform gene-level data preparatgion
# source("UPF1_NMDRHT_Transcript_Level.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load Main transcript NMDRHT table
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

#
##
# Metabolic mRNA labeling data (4SU) transcript-level --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to run metabolic labeling analysis via bakR
# source("UPF1_NMDRHT_EZbakR_NMDRHT.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

# Load EZbakR combined results
UPF1_NMDRHT_EZbakR_TEC_combined <- read_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))

#
##
# Protein-level data preparation --------------------------------------------
##
#

# *IMPORTANT*: running this sub-script requires that all required tools, data and analyses outputs are present in the respective folders - thats why it is commented out by default

# Please uncomment the next command to perform gene-level data preparatgion
# source("UPF1_NMDRHT_Protein_Level.R")

# Comment: alternatively, load or import checkpoint data or final outputs within the script above

#
##
###
# Revision_1 (Rev_1) Plots and Tables ----------------------------------------------------------------
###
##
#

###
## Rev_1 - PCR analyses -----------------------------------------------------------
###

# Analyses and plots for (q)PCRs
source("UPF1_NMDRHT_PCR_analysis.R")

###
## Rev_1 - Figure 1 -----------------------------------------------------------
###

# Analyses and plots for Figure 1 + Supplements
source("UPF1_NMDRHT_Figure1.R")

###
## Rev_1 - Figure 2 -----------------------------------------------------------
###

# Analyses and plots for Figure 2 + Supplements
source("UPF1_NMDRHT_Figure2.R")

###
## Rev_1 - Figure 3 -----------------------------------------------------------
###

# Analyses and plots for Figure 3 + Supplements
source("UPF1_NMDRHT_Figure3.R")

###
## Rev_1 - Figure 4 -----------------------------------------------------------
###

# Analyses and plots for Figure 4 + Supplements
source("UPF1_NMDRHT_Figure4.R")

###
## Rev_1 - Figure 5 -----------------------------------------------------------
###

# Analyses and plots for Figure 5 + Supplements
source("UPF1_NMDRHT_Figure5.R")

###
## Rev_1 - Figure 6 -----------------------------------------------------------
###

# Analyses and plots for Figure 6 + Supplements
source("UPF1_NMDRHT_Figure6.R")

###
## Rev_1 - Figure 7 -----------------------------------------------------------
###

# Analyses and plots for Figure 7 + Supplements
source("UPF1_NMDRHT_Figure7.R")

###
## Rev_1 - Supplemental Tables -----------------------------------------------------------
###

# Generate Supplemental Tables
source("UPF1_NMDRHT_SupplementalTables.R")

###
## Rev_1 - Revision-specific Analysis -----------------------------------------------------------
###

# Generate Supplemental Tables
source("UPF1_NMDRHT_RevisionAnalyses.R")

# SessionInfo to file
writeLines(capture.output(sessionInfo()), paste0("UPF1_NMDRHT_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))