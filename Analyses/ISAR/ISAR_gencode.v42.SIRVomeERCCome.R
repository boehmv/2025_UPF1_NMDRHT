#!/usr/bin/env Rscript


# Title: ISAR_gencode.v42.SIRVomeERCCome
# Objective: Standard IsoformSwitchAnalyseR (ISAR) pipeline giving differential transcript usage (DTU) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(IsoformSwitchAnalyzeR)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Load simplified gencode annotation
gtf_gencode_df_short <- read_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

# Get list of accepted biotypes
accepted_biotypes <- gtf_gencode_df_short %>% distinct(gene_type) %>% pull(gene_type)

# Get design table
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)
myDesign <- data.frame(
    sampleID = samples$sample,
    condition = samples$condition
)

# condition string is the second, gets converted to character
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Compile sampleVector from experiment file
myFiles <- file.path(mydir, "Salmon", samples$sample, "quant.sf")
names(myFiles) <- samples$sample

# Import quantifications
salmonQuant <- importIsoformExpression(
    sampleVector = myFiles,
    addIsofomIdAsColumn = TRUE
)

# Get comparison data frame
myComparison <- data.frame(
    condition_1 = "control",
    condition_2 = cond[-1]
)

# Generate switchlist
SwitchList <- importRdata(
		isoformCountMatrix   = salmonQuant$counts,
		isoformRepExpression = salmonQuant$abundance,
		designMatrix         = myDesign,
		detectUnwantedEffects = FALSE,
		addAnnotatedORFs     = TRUE,
		onlyConsiderFullORF = FALSE,
		removeNonConvensionalChr = TRUE,
		isoformExonAnnoation = file.path(ref_dir, "Gencode",  "gencode.v42.SIRVomeERCCome.annotation.gtf"),
		comparisonsToMake= myComparison,	    
		showProgress = TRUE,
		ignoreAfterBar = TRUE,
    		ignoreAfterSpace = TRUE,
    		ignoreAfterPeriod = FALSE,	# Set to FALSE for Gencode
    		removeTECgenes = TRUE		# If set to TRUE, spike_ins need "gene_name" and "gene_type"
)

# What is in the Switchlist
SwitchList

setdff <- setdiff

# Filtering
SwitchList_filt <- preFilter(
  SwitchList,
  geneExpressionCutoff = 1, # FPMK threshold
  isoformExpressionCutoff = 0, # FPMK threshold
  acceptedGeneBiotype = accepted_biotypes,
  IFcutoff=0.01,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  quiet=FALSE
)

# Testing for Isoform Switches via DEXSeq
SwitchList_filt_Analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchList_filt,
    reduceToSwitchingGenes=TRUE
)

SwitchList_filt_Analyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_filt_Analyzed,
  quiet=FALSE
)

# Set WD
setwd((file.path(mydir, "ISAR")))

# Summarize switching features
extractSwitchSummary(SwitchList_filt_Analyzed)

# Predicting Switch Consequences
SwitchList_filt_Analyzed <- analyzeSwitchConsequences(
    SwitchList_filt_Analyzed,
    consequencesToAnalyze = c('NMD_status'), 
    dIFcutoff = 0.1, 
    showProgress=TRUE
)

# Summarize switching features without consequences
extractSwitchSummary(SwitchList_filt_Analyzed, dIFcutoff = 0.1, filterForConsequences = FALSE)

# Summarize switching features with consequences
extractSwitchSummary(SwitchList_filt_Analyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)

write.csv(SwitchList_filt_Analyzed$isoformFeatures, file="SwitchList_filt_Analyzed.csv")
saveRDS(SwitchList_filt_Analyzed, file = "SwitchList_filt_Analyzed.rds")

writeLines(capture.output(sessionInfo()), paste0(mydir, "/ISAR/ISAR_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
