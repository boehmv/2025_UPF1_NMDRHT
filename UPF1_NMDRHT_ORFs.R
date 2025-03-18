#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_ORFs
# Objective: Supplement NMDRHT annotation with ORFs identified from Ribo-Seq data (0h & 12h IAA on HCT116 N-AID-UPF1) with ORFquant and RiboTISH
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# *NOTE* information up-front: 50-nt rule is counted from the first nucleotide after the stop codon to the first nucleotide of the intron
# ORFquant and factR do this differently
##

# Setup -------------------------------------------------------------------

library("devtools")

# Install ORFquant and Ribo-seQC if necessary
# # ORFquant 1.0.2 - 481ec99
# install_github(repo = "lcalviell/ORFquant",
#                ref = "481ec99",
#                force = TRUE)


library("RiboseQC")
library("ORFquant")

#
##
###
# NMDRHT -----------------------------------------------------------------
###
##
#

# Load shortend NMDRHT information
NMDRHT.v1.1_tbl <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_short.csv")

# load annotation information
NMDRHT.v1.1_GTF <- rtracklayer::readGFF("Resources/NMDRHT/NMDRHT.v1.1.sort.gtf",
                                              filter = list(type=c("gene",
                                                                   "transcript")))
# Convert to tibble
NMDRHT.v1.1_GTF_tbl <- dplyr::as_tibble(NMDRHT.v1.1_GTF)

## Annotation --------------------------------------------------------------
# Prepare annotation file

ORFquant::prepare_annotation_files(annotation_directory = "Resources/ORFquant",
                         twobit_file = "Resources/ORFquant/GRCh38.primary_assembly.genome.2bit",
                         gtf_file = "Resources/NMDRHT/NMDRHT.v1.1.sort.gtf",
                         scientific_name = "Homo.sapiens",
                         annotation_name = "NMDRHT.v1.1",
                         export_bed_tables_TxDb = F,
                         forge_BSgenome = T,
                         create_TxDb = T)

# Load annotation
ORFquant::load_annotation("Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot")

## Run RiboseQC  -------------------------------------------------------
RiboseQC::RiboseQC_analysis(annotation_file="Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot",
                  bam_files = c("/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM_NMDRegHumanTxome/HCT_N_AID_UPF1_0h_IAA_merged.bam",
                                "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/BAM_NMDRegHumanTxome/HCT_N_AID_UPF1_12h_IAA_merged.bam"),
                  dest_names = c("/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/NMDRegHumanTxome/HCT_N_AID_UPF1_0h_IAA",
                                 "/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd/NMDRegHumanTxome/HCT_N_AID_UPF1_12h_IAA"),
                  report_file = "Resources/ORFquant/RiboseQC_merged.html",
                  sample_names = c("HCT_N_AID_UPF1_0h",
                                   "HCT_N_AID_UPF1_12h"),
                  write_tmp_files = F)

##
## Run ORFquant  -------------------------------------------------------
##

##
### control (HCT_N_AID_UPF1_0h_IAA)   -------------------------------------------------------
##

run_ORFquant(for_ORFquant_file = "Resources/ORFquant/HCT_N_AID_UPF1_0h_IAA_for_SaTAnn",
             annotation_file = "Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot",
             n_cores = 12)

# Generate plots
plot_ORFquant_results(for_ORFquant_file = "Resources/ORFquant/HCT_N_AID_UPF1_0h_IAA_for_SaTAnn",
                      ORFquant_output_file = "Resources/ORFquant/HCT_N_AID_UPF1_0h_IAA_for_SaTAnn_final_ORFquant_results",
                      annotation_file = "Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot",
                      output_plots_path = "Resources/ORFquant/")

# Load Control Data for analysis (if required)
load("Resources/ORFquant/HCT_N_AID_UPF1_0h_IAA_for_SaTAnn_final_ORFquant_results")

library(tidyverse)

# Obtain transcript- and gene-specific information
ORFquant_results$ORFs_tx -> ORFs_tx
ORFquant_results$ORFs_gen -> ORFs_gen

# Extract ORF region from gene-specific information
# Match by ORF_id_tr
# Fix missing stop codon
ORF_gen_df <- dplyr::as_tibble(ORFs_gen) %>% 
  mutate(ORF_id_tr = ORFs_gen@ranges@NAMES) %>% 
  group_by(ORF_id_tr) %>% 
  mutate(start = min(start),
         end = max(end)) %>% 
  mutate(ORF_region = paste0(seqnames,":",start,"-",end)) %>% 
  dplyr::distinct(ORF_id_tr, ORF_region)

# Extract other parameters from transcript-specific information
ORF_tx_df <- dplyr::tibble(gene_id = ORFs_tx$gene_id,
                           gene_biotype = ORFs_tx$gene_biotype,
                           gene_name =  ORFs_tx$gene_name,
                           transcript_id = ORFs_tx$transcript_id,
                           transcript_biotype = ORFs_tx$transcript_biotype,
                           ORF_id_tr = ORFs_tx$ORF_id_tr,
                           dplyr::as_tibble(ORFs_tx$region),
                           ave_pct_fr = ORFs_tx$ave_pct_fr,
                           pct_fr = ORFs_tx$pct_fr,
                           ave_pct_fr_st = ORFs_tx$ave_pct_fr_st,
                           pct_fr_st = ORFs_tx$pct_fr_st,
                           pval = ORFs_tx$pval,
                           pval_uniq = ORFs_tx$pval_uniq,
                           P_sites_raw = ORFs_tx$P_sites_raw,
                           P_sites_raw_uniq = ORFs_tx$P_sites_raw_uniq,
                           P_sites_raw_uniq_mm = ORFs_tx$P_sites_raw_uniq_mm,
                           P_sites = ORFs_tx$P_sites,
                           ORF_pct_P_sites = ORFs_tx$ORF_pct_P_sites,
                           ORF_pct_P_sites_pN = ORFs_tx$ORF_pct_P_sites_pN,
                           compatible_biotype = ORFs_tx$compatible_biotype,
                           NC_protein_isoform = ORFs_tx$NC_protein_isoform,
                           ORF_category_Tx = ORFs_tx$ORF_category_Tx,
                           ORF_category_Tx_compatible = ORFs_tx$ORF_category_Tx_compatible,
                           ORF_category_Gen = ORFs_tx$ORF_category_Gen,
                           NMD_candidate = ORFs_tx$NMD_candidate,
                           Distance_to_lastExEx = ORFs_tx$Distance_to_lastExEx,
                           ORFs_pM = ORFs_tx$ORFs_pM,
                           longest_ORF = as.character(ORFs_tx$longest_ORF),
                           Protein = as.character(ORFs_tx$Protein))

# Join both tibbles
HCT_N_AID_UPF1_0h_ORF_tx_df <- ORF_tx_df %>% left_join(ORF_gen_df)

#### Expand ORFquant ---------------------------------------------------------
# ORFquant reduces the output to less transcripts, but reports "compatible" transcript per ORF
# Since we are interested in all transcripts, expand the results to all compatible transcripts!

# Extract metadata from transcript-information
ORFquant_results_mcols_ORFs_tx <- mcols(ORFquant_results$ORFs_tx, use.names = FALSE)

# Get list of ORF_id_tr -to- compatible transcripts matching
ORFquant_results_mcols_ORFs_tx_compatible_tx_List <- ORFquant_results_mcols_ORFs_tx@listData$compatible_with
names(ORFquant_results_mcols_ORFs_tx_compatible_tx_List) <- ORFquant_results_mcols_ORFs_tx@listData$ORF_id_tr
ORFquant_results_mcols_ORFs_tx_compatible_tx <- as_tibble(ORFquant_results_mcols_ORFs_tx_compatible_tx_List) %>% 
  dplyr::rename("compatible_tx" = "value")

# Get distance to last exon-exon-junction for each compatible transcript
ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List <- ORFquant_results_mcols_ORFs_tx@listData$Distance_to_lastExEx_compatible_txs
names(ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List) <- ORFquant_results_mcols_ORFs_tx@listData$ORF_id_tr
ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs <- as_tibble(ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List) %>% 
  dplyr::rename("Distance_to_lastExEx_compatible" = "value")

# Combine both information sources
ORFquant_results_mcols_combined <- ORFquant_results_mcols_ORFs_tx_compatible_tx %>% 
  mutate(Distance_to_lastExEx_compatible = ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs$Distance_to_lastExEx_compatible) %>% 
  dplyr::select(-group) %>% 
  dplyr::rename("ORF_id_tr" = "group_name")

# Obtain expanded tibble with all compatible-transcript information
HCT_N_AID_UPF1_0h_ORF_tx_df_expand <- HCT_N_AID_UPF1_0h_ORF_tx_df %>% 
  dplyr::select(-c(gene_id,
                   gene_biotype,
                   gene_name,
                   transcript_id,
                   transcript_biotype,
                   seqnames,
                   start,
                   end,
                   width,
                   compatible_biotype,
                   NMD_candidate,
                   Distance_to_lastExEx,
                   longest_ORF)) %>% 
  left_join(ORFquant_results_mcols_combined) %>% 
  # subtract 4 nucleotides from distance - reason: ORFquant counts the distance from the last nt of the CDS (excluding the stop codon)
  # more conservative: count from first nt after the stop codon!
  mutate(Distance_to_lastExEx_compatible = as.numeric(Distance_to_lastExEx_compatible)-4) %>% 
  mutate(NMD_status_ORFquant = case_when(Distance_to_lastExEx_compatible > 50 ~ TRUE,
                                         Distance_to_lastExEx_compatible <= 50 ~ FALSE),
         .after=Distance_to_lastExEx_compatible) %>% 
  tidyr::separate(compatible_tx, into=c("transcript_id", NA, NA), remove = FALSE) %>% 
  left_join(NMDRHT.v1.1_GTF_tbl) %>% 
  relocate(gene_id,
           transcript_id,
           gene_name,
           transcript_name,
           gene_biotype,
           transcript_biotype) 

##
### UPF1 depletion (HCT_N_AID_UPF1_12h_IAA)   -------------------------------------------------------
##

run_ORFquant(for_ORFquant_file = "Resources/ORFquant/HCT_N_AID_UPF1_12h_IAA_for_SaTAnn",
             annotation_file = "Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot",
             n_cores = 12)

# Generate plots
plot_ORFquant_results(for_ORFquant_file = "Resources/ORFquant/NMDRegHumanTxome/HCT_N_AID_UPF1_12h_IAA_for_SaTAnn",
                      ORFquant_output_file = "Resources/ORFquant/HCT_N_AID_UPF1_12h_IAA_for_SaTAnn_final_ORFquant_results",
                      annotation_file = "Resources/ORFquant/NMDRHT.v1.1.sort.gtf_Rannot",
                      output_plots_path = "Resources/ORFquant/")

#### Create combined 0h&12h HTML report ------------------------------------------------------
create_ORFquant_html_report(input_files = c("Resources/ORFquant/HCT_N_AID_UPF1_0h_IAA_for_SaTAnn_ORFquant_plots_RData",
                                            "Resources/ORFquant/HCT_N_AID_UPF1_12h_IAA_for_SaTAnn_ORFquant_plots_RData"),
                            input_sample_names = c("HCT_N_AID_UPF1_0h",
                                                   "HCT_N_AID_UPF1_12h"), 
                            output_file= "Resources/ORFquant/HCT_N_AID_UPF1_ORFquant_report_merged.html")

# Load Control Data for analysis (if required)
load("Resources/ORFquant/HCT_N_AID_UPF1_12h_IAA_for_SaTAnn_final_ORFquant_results")

# Obtain transcript- and gene-specific information
ORFquant_results$ORFs_tx -> ORFs_tx
ORFquant_results$ORFs_gen -> ORFs_gen

# Extract ORF region from gene-specific information
# Match by ORF_id_tr
ORF_gen_df <- dplyr::as_tibble(ORFs_gen) %>% 
  mutate(ORF_id_tr = ORFs_gen@ranges@NAMES) %>% 
  group_by(ORF_id_tr) %>% 
  mutate(start = min(start),
         end = max(end)) %>% 
  mutate(ORF_region = paste0(seqnames,":",start,"-",end)) %>% 
  dplyr::distinct(ORF_id_tr, ORF_region)

# Extract other parameters from transcript-specific information
ORF_tx_df <- dplyr::tibble(gene_id = ORFs_tx$gene_id,
                    gene_biotype = ORFs_tx$gene_biotype,
                    gene_name =  ORFs_tx$gene_name,
                    transcript_id = ORFs_tx$transcript_id,
                    transcript_biotype = ORFs_tx$transcript_biotype,
                    ORF_id_tr = ORFs_tx$ORF_id_tr,
                    dplyr::as_tibble(ORFs_tx$region),
                    ave_pct_fr = ORFs_tx$ave_pct_fr,
                    pct_fr = ORFs_tx$pct_fr,
                    ave_pct_fr_st = ORFs_tx$ave_pct_fr_st,
                    pct_fr_st = ORFs_tx$pct_fr_st,
                    pval = ORFs_tx$pval,
                    pval_uniq = ORFs_tx$pval_uniq,
                    P_sites_raw = ORFs_tx$P_sites_raw,
                    P_sites_raw_uniq = ORFs_tx$P_sites_raw_uniq,
                    P_sites_raw_uniq_mm = ORFs_tx$P_sites_raw_uniq_mm,
                    P_sites = ORFs_tx$P_sites,
                    ORF_pct_P_sites = ORFs_tx$ORF_pct_P_sites,
                    ORF_pct_P_sites_pN = ORFs_tx$ORF_pct_P_sites_pN,
                    compatible_biotype = ORFs_tx$compatible_biotype,
                    NC_protein_isoform = ORFs_tx$NC_protein_isoform,
                    ORF_category_Tx = ORFs_tx$ORF_category_Tx,
                    ORF_category_Tx_compatible = ORFs_tx$ORF_category_Tx_compatible,
                    ORF_category_Gen = ORFs_tx$ORF_category_Gen,
                    NMD_candidate = ORFs_tx$NMD_candidate,
                    Distance_to_lastExEx = ORFs_tx$Distance_to_lastExEx,
                    ORFs_pM = ORFs_tx$ORFs_pM,
                    longest_ORF = as.character(ORFs_tx$longest_ORF),
                    Protein = as.character(ORFs_tx$Protein))

# Join both tibbles
HCT_N_AID_UPF1_12h_ORF_tx_df <- ORF_tx_df %>% left_join(ORF_gen_df)

#### Expand ORFquant ---------------------------------------------------------
# ORFquant reduces the output to less transcripts, but reports "compatible" transcript per ORF
# Since we are interested in all transcripts, expand the results to all compatible transcripts!

# Extract metadata from transcript-information
ORFquant_results_mcols_ORFs_tx <- mcols(ORFquant_results$ORFs_tx, use.names = FALSE)

# Get list of ORF_id_tr -to- compatible transcripts matching
ORFquant_results_mcols_ORFs_tx_compatible_tx_List <- ORFquant_results_mcols_ORFs_tx@listData$compatible_with
names(ORFquant_results_mcols_ORFs_tx_compatible_tx_List) <- ORFquant_results_mcols_ORFs_tx@listData$ORF_id_tr
ORFquant_results_mcols_ORFs_tx_compatible_tx <- as_tibble(ORFquant_results_mcols_ORFs_tx_compatible_tx_List) %>% 
  dplyr::rename("compatible_tx" = "value")

# Get distance to last exon-exon-junction for each compatible transcript
ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List <- ORFquant_results_mcols_ORFs_tx@listData$Distance_to_lastExEx_compatible_txs
names(ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List) <- ORFquant_results_mcols_ORFs_tx@listData$ORF_id_tr
ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs <- as_tibble(ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs_List) %>% 
  dplyr::rename("Distance_to_lastExEx_compatible" = "value")

# Combine both information sources
ORFquant_results_mcols_combined <- ORFquant_results_mcols_ORFs_tx_compatible_tx %>% 
  mutate(Distance_to_lastExEx_compatible = ORFquant_results_mcols_ORFs_tx_Distance_to_lastExEx_compatible_txs$Distance_to_lastExEx_compatible) %>% 
  dplyr::select(-group) %>% 
  dplyr::rename("ORF_id_tr" = "group_name")

# Obtain expanded tibble with all compatible-transcript information
HCT_N_AID_UPF1_12h_ORF_tx_df_expand <- HCT_N_AID_UPF1_12h_ORF_tx_df %>% 
  dplyr::select(-c(gene_id,
                   gene_biotype,
                   gene_name,
                   transcript_id,
                   transcript_biotype,
                   seqnames,
                   start,
                   end,
                   width,
                   compatible_biotype,
                   NMD_candidate,
                   Distance_to_lastExEx,
                   longest_ORF)) %>% 
  left_join(ORFquant_results_mcols_combined) %>% 
  # subtract 4 nucleotides from distance - reason: ORFquant counts the distance from the last nt of the CDS (excluding the stop codon)
  # more conservative: count from first nt after the stop codon!
  mutate(Distance_to_lastExEx_compatible = as.numeric(Distance_to_lastExEx_compatible)-4) %>%
  mutate(NMD_status_ORFquant = case_when(Distance_to_lastExEx_compatible > 50 ~ TRUE,
                                         Distance_to_lastExEx_compatible <= 50 ~ FALSE),
         .after=Distance_to_lastExEx_compatible) %>% 
  tidyr::separate(compatible_tx, into=c("transcript_id", NA, NA), remove = FALSE) %>% 
  left_join(NMDRHT.v1.1_GTF_tbl) %>% 
  relocate(gene_id,
           transcript_id,
           gene_name,
           transcript_name,
           gene_biotype,
           transcript_biotype) 

##
### Combine 0h&12h ----------------------------------------------------------
##

# Comnbine 0h and 12h IAA ORFquant data
HCT_N_AID_UPF1_ORFquant <- HCT_N_AID_UPF1_0h_ORF_tx_df_expand %>% 
  mutate(condition = "HCT_N_AID_UPF1_0h") %>% 
  bind_rows(HCT_N_AID_UPF1_12h_ORF_tx_df_expand %>% 
              mutate(condition = "HCT_N_AID_UPF1_12h")) %>% 
  mutate(ORF_length = str_length(Protein), .after=Protein)

##
#### save information ----------------------------------------------------------
##
save(HCT_N_AID_UPF1_ORFquant,
     file = paste0("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant.rds"))

HCT_N_AID_UPF1_ORFquant %>% write_csv("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant.csv")

##
# Ribo-TISH ---------------------------------------------------------------
##

# version 0.2.7

# Run Ribo-TISH predict on replicates - only canonical start codons!
## Was done in independent script (CRTP_QTI_V003_MODIFIED.sh)
# Use it to independently confirm ORFquant predicted ORFs
# Import prediction files (tab-separated)
# Fix genomic posiitons to be compatible with ORFquant (plus stand = start-1 & end-3; minus strand = start-4)

## Control (HCT_N_AID_UPF1_0h) ---------------------------------------------------------------
HCT_N_AID_UPF1_0h_RiboTISH_NMDRHT <- read_delim("Resources/RiboTISH/control_riboseq_only.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE) %>% 
  dplyr::rename("gene_id" = "Gid",
                "transcript_id" = "Tid",
                "gene_name" = "Symbol",
                "gene_biotype" = "GeneType",
                "ORF_region" = "GenomePos") %>% 
  separate(ORF_region, into=c("chr", "coords", "strand"), sep=":") %>% 
  mutate(ORF_region = paste0(chr, ":", coords), .before=chr) %>% 
  separate(coords, into=c("start_coord","end_coord")) %>% 
  mutate(ORF_region_ORFquantComp = case_when(strand == "+" ~ paste0(chr, ":", as.numeric(start_coord)+1, "-", as.numeric(end_coord)-3),
                                             strand == "-" ~ paste0(chr, ":", as.numeric(start_coord)+4, "-", as.numeric(end_coord))), .before=ORF_region) %>% 
  mutate(condition="HCT_N_AID_UPF1_0h")

## UPF1 depletion (HCT_N_AID_UPF1_12h) ---------------------------------------------------------------
HCT_N_AID_UPF1_12h_RiboTISH_NMDRHT <- read_delim("Resources/RiboTISH/UPF1_riboseq_only.txt", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE) %>% 
  dplyr::rename("gene_id" = "Gid",
                "transcript_id" = "Tid",
                "gene_name" = "Symbol",
                "gene_biotype" = "GeneType",
                "ORF_region" = "GenomePos") %>% 
  separate(ORF_region, into=c("chr", "coords", "strand"), sep=":") %>% 
  mutate(ORF_region = paste0(chr, ":", coords), .before=chr) %>% 
  separate(coords, into=c("start_coord","end_coord")) %>% 
  mutate(ORF_region_ORFquantComp = case_when(strand == "+" ~ paste0(chr, ":", as.numeric(start_coord)+1, "-", as.numeric(end_coord)-3),
                                             strand == "-" ~ paste0(chr, ":", as.numeric(start_coord)+4, "-", as.numeric(end_coord))), .before=ORF_region) %>% 
  mutate(condition="HCT_N_AID_UPF1_12h")

## Combine both
HCT_N_AID_UPF1_RiboTISH_combined <- HCT_N_AID_UPF1_0h_RiboTISH_NMDRHT %>% 
  bind_rows(HCT_N_AID_UPF1_12h_RiboTISH_NMDRHT) %>% 
  left_join(NMDRHT.v1.1_GTF_tbl %>% dplyr::select(-c(ORF_type, NMD_status, stop_to_lastEJ, num_of_downEJs, UTR3_length)))

##
#### save information ----------------------------------------------------------
##
save(HCT_N_AID_UPF1_RiboTISH_combined,
     file = paste0("Resources/RiboTISH/HCT_N_AID_UPF1_RiboTISH.rds"))

HCT_N_AID_UPF1_RiboTISH_combined %>% write_csv("Resources/RiboTISH/HCT_N_AID_UPF1_RiboTISH.csv")

## Compare ORFquant & Ribo-TISH ---------------------------------------------------------------
# Join ORFquant and Ribo-TISH
# If ORFquant ORF was never found in Ribo-TISH -> flag it as RiboTISH=FALSE
HCT_N_AID_UPF1_ORFquant_RiboTish <- HCT_N_AID_UPF1_ORFquant %>% 
  left_join(HCT_N_AID_UPF1_RiboTISH_combined %>% dplyr::select(transcript_id, gene_name,
                                                               ORF_region_ORFquantComp,
                                                               strand,
                                                               TisType,
                                                               RiboPStatus,
                                                               FrameQvalue,
                                                               InFrameCount,
                                                               condition) %>% 
              dplyr::rename("ORF_region" = "ORF_region_ORFquantComp")) %>% 
  relocate(TisType,
           RiboPStatus,
           FrameQvalue,
           InFrameCount,.before = ave_pct_fr) %>% 
  group_by(ORF_id_tr) %>% 
  mutate(RiboTISH_perc = sum(!is.na(RiboPStatus))/n(),.before = TisType) %>% 
  ungroup() %>% 
  mutate(RiboTISH = case_when(RiboTISH_perc != 0 ~ TRUE,
                              TRUE ~ FALSE),.before = RiboTISH_perc) 

# Save information for plots
HCT_N_AID_UPF1_ORFquant_RiboTish %>% write_csv("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish.csv")

## Stats -----------------------------------------------------------

# How many ORFquant ORFs are also found in Ribo-TISH (per condition)
HCT_N_AID_UPF1_ORFquant_RiboTish %>% 
  distinct(ORF_region,condition, RiboTISH) %>% 
  dplyr::count(condition,RiboTISH) %>% 
  arrange(desc(n)) 

##
# NMD reason --------------------------------------------------------------
##

# Load NMD-sensitivity analysis on GENCODE canonical transcripts (from UPF1_NMDRHT_Annotation.R)
# Those were used to predict ORFs - which are treated by ORFquant as "annotated"
gencode.v42.annotation.canonical_NMDprediction.out <- read_csv(("Resources/GENCODE/gencode.v42.annotation.canonical_NMDprediction.csv"))

# Filter for NMD-predicted canonical transcripts
# 378 canonical transcripts are predicted to be "natural" NMD targets!
gencode.v42.annotation.canonical_NMDprediction.out_NMDplus <- gencode.v42.annotation.canonical_NMDprediction.out %>% 
  filter(is_NMD == TRUE)

##
# Define NMD reason
##
# 1) If NMD status is FALSE = no NMD reason ~ "none"
# 2) If NMD is TRUE and ORF is "annotated" (coming from ORFanage and GENCODE-canonical) & ORF type is exact match *and* 
#    is already NMD = TRUE in GENCODE ~ "annotated_NMD"
# 3) If NMD is TRUE and ORF is "annotated" & ORF type is exact match -> must be 3'UTR alternative splicing event ~ "AS_NMD_UTR3"
# 4) If NMD is TRUE and ORF is "annotated" & ORF type is start match -> means different stop codon - probably AS ~ "AS_NMD"
# 5) If NMD is TRUE and ORF is NOT "annotated" ~ use ORF_category from ORFquant as reason

HCT_N_AID_UPF1_ORFquant_RiboTish_NMD <- HCT_N_AID_UPF1_ORFquant_RiboTish %>% 
  mutate(NMD_reason = case_when(NMD_status_ORFquant == FALSE ~ "none",
                                NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "exact_match" & GENCODE_transcript_id %in% gencode.v42.annotation.canonical_NMDprediction.out_NMDplus$transcript_id ~ "annotated_NMD",
                                NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "exact_match" ~ "AS_NMD_UTR3",
                                NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "start_match" ~ "AS_NMD",
                                NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "start_other" ~ "ORFanage_altStart",
                                NMD_status_ORFquant == TRUE ~ ORF_category_Tx,
  ),
  .after = NMD_status_ORFquant) %>% 
  # modify ORFanage-based parameters
  dplyr::rename("orfanage_ORF_type"="ORF_type",
                "orfanage_start_codon_type"="start_codon_type",
                "orfanage_stop_codon_type"="stop_codon_type",
                "orfanage_ORF_id"="ORF_id",
                "orfanage_NMD_status"="NMD_status",
                "orfanage_stop_to_lastEJ"="stop_to_lastEJ",
                "orfanage_num_of_downEJs"="num_of_downEJs",
                "orfanage_UTR3_length"="UTR3_length") %>% 
  relocate(orfanage_status,
           orfanage_duplicity,
           orfanage_template,
           orfanage_template_source,
           orfanage_ORF_type,
           orfanage_start_codon_type,
           orfanage_stop_codon_type,
           orfanage_ORF_id,
           orfanage_NMD_status,
           orfanage_stop_to_lastEJ,
           orfanage_num_of_downEJs,
           orfanage_UTR3_length,
           .after = last_col())

# Check for ORF_annotated without NMD_reason
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD %>% 
  filter(NMD_status_ORFquant == TRUE) %>% 
  dplyr::count(NMD_reason)

# Export final tibble
write_csv(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD,
          file="Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish_NMD.csv")

# read file in if neccesary
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD <- read_csv(file="Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish_NMD.csv")

## Stats -------------------------------------------------------------------

# How many transcripts with detected ORFs per condition
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD %>% 
  distinct(transcript_id, condition, RiboTISH) %>% 
  dplyr::count(condition, RiboTISH)

#
##
###
# NMDRHT ORF annotation ---------------------------------------------------
###
##
#

# Objective: Identify the best high-quality (found in both ORFquant & Ribo-TISH) ORF per NMDRHT transcript
# Favor NMD-activating (50-nt rule == TRUE) ORFs in "mixed" situations (where both NMD-sensitive and -insensitive ORFs are identified)


## All ORFquant ORFs -------------------------------------------------------

##
## *NOTE* - ORFquant-only ORFs are filtered out from most later analyses - only those that are also found in Ribo-TISH are kept! (see below; proteomics will include them again!)
##

# Calculate transcript-wise ORF_percentage (ORF length normalized)
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_tx_ORF_pct_all <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD %>% 
  mutate(P_sites_pN = P_sites/ORF_length, .after = P_sites) %>% 
  group_by(transcript_id, condition) %>% 
  mutate(tx_P_sites_sum_pN = sum(P_sites_pN), .after = P_sites_pN) %>% 
  ungroup() %>% 
  mutate(tx_ORF_pct_P_sites_pN = (P_sites_pN/tx_P_sites_sum_pN)*100, .after = tx_P_sites_sum_pN)

# Summarize tx_ORF_pct_P_sites_pN depending on NMD status for each transcript
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_tx_ORF_pct_all %>% 
  group_by(transcript_id, condition) %>% 
  mutate(NMD_tx_ORF_pct_P_sites_pN = sum(tx_ORF_pct_P_sites_pN[NMD_status_ORFquant == TRUE]),
         Non_NMD_tx_ORF_pct_P_sites_pN = sum(tx_ORF_pct_P_sites_pN[NMD_status_ORFquant == FALSE]),
         .after = tx_ORF_pct_P_sites_pN) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(ORFquant_NMD_tx_status = case_when(NMD_tx_ORF_pct_P_sites_pN == 100 ~ "NMD",
                                            NMD_tx_ORF_pct_P_sites_pN == 0 ~ "coding",
                                            TRUE ~ "mixed")) %>% 
  ungroup()

# Check for overlap between 0h and 12h
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all %>% 
  dplyr::select(gene_id,transcript_id,gene_name,transcript_name, condition, ORF_id_tr, NMD_status_ORFquant, ORFquant_NMD_tx_status, NMD_reason, tx_P_sites_sum_pN, P_sites_pN) %>%
  group_by(condition,transcript_id, ORF_id_tr) %>% 
  slice_max(tx_P_sites_sum_pN, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  pivot_wider(names_from = condition,
              values_from =c(ORFquant_NMD_tx_status, tx_P_sites_sum_pN, P_sites_pN),
              names_glue = "{condition}_{.value}") %>% 
  # pseudocount of 0.001 for not detected ORFs to allow log2 calculation later on
  mutate(HCT_N_AID_UPF1_0h_valid = case_when(is.na(HCT_N_AID_UPF1_0h_P_sites_pN) ~ FALSE,
                                             TRUE ~ TRUE),
         HCT_N_AID_UPF1_12h_valid = case_when(is.na(HCT_N_AID_UPF1_12h_P_sites_pN) ~ FALSE,
                                             TRUE ~ TRUE)) %>% 
  replace_na(list(HCT_N_AID_UPF1_0h_P_sites_pN = 0.001, HCT_N_AID_UPF1_12h_P_sites_pN = 0.001)) %>% 
  mutate(log2FC_P_sites_pN = log2(HCT_N_AID_UPF1_12h_P_sites_pN/HCT_N_AID_UPF1_0h_P_sites_pN)) %>% 
  mutate(ORFquant_NMD_tx_status_combined = case_when(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status == HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     is.na(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status) ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     is.na(HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status) ~ HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status,
                                                     HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN >= HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     TRUE ~ HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status)) %>% 
  mutate(ORFquant_NMD_tx_status_source = case_when(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status == HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status ~ "ORFquant_HCT_N_AID_UPF1_0h&12h",
                                                   is.na(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status) ~ "ORFquant_HCT_N_AID_UPF1_12h",
                                                   is.na(HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status) ~ "ORFquant_HCT_N_AID_UPF1_0h",
                                                   HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN >= HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN ~ "ORFquant_HCT_N_AID_UPF1_12h",
                                                   TRUE ~ "ORFquant_HCT_N_AID_UPF1_0h"))

# Save data for e.g. usage in Proteomics
save(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique,
     file = paste0("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique.rds"))

#
## HQ overlapping ORFs -------------------------------------------------------
#

# Calculate **ONLY HQ = Ribo-TISH-supported** transcript-wise ORF_percentage (ORF length normalized)
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_tx_ORF_pct_HQ <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD %>% 
  filter(RiboTISH == TRUE) %>% 
  mutate(P_sites_pN = P_sites/ORF_length, .after = P_sites) %>% 
  group_by(transcript_id, condition) %>% 
  mutate(tx_P_sites_sum_pN = sum(P_sites_pN), .after = P_sites_pN) %>% 
  ungroup() %>% 
  mutate(tx_ORF_pct_P_sites_pN = (P_sites_pN/tx_P_sites_sum_pN)*100, .after = tx_P_sites_sum_pN)

# Summarize tx_ORF_pct_P_sites_pN depending on NMD status for each transcript
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_tx_ORF_pct_HQ %>% 
  group_by(transcript_id, condition) %>% 
  mutate(NMD_tx_ORF_pct_P_sites_pN = sum(tx_ORF_pct_P_sites_pN[NMD_status_ORFquant == TRUE]),
         Non_NMD_tx_ORF_pct_P_sites_pN = sum(tx_ORF_pct_P_sites_pN[NMD_status_ORFquant == FALSE]),
         .after = tx_ORF_pct_P_sites_pN) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(ORFquant_NMD_tx_status = case_when(NMD_tx_ORF_pct_P_sites_pN == 100 ~ "NMD",
                                            NMD_tx_ORF_pct_P_sites_pN == 0 ~ "coding",
                                            TRUE ~ "mixed")) %>% 
  ungroup()

# Check for overlap between 0h and 12h
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_unique <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF %>% 
  dplyr::select(gene_id,transcript_id,gene_name,transcript_name, condition, ORFquant_NMD_tx_status, tx_P_sites_sum_pN) %>%
  group_by(condition,transcript_id) %>% 
  slice_max(tx_P_sites_sum_pN, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  pivot_wider(names_from = condition,
              values_from =c(ORFquant_NMD_tx_status, tx_P_sites_sum_pN),
              names_glue = "{condition}_{.value}") %>% 
  # pseudocount of 0.001 for not detected ORFs to allow log2 calculation later on
  replace_na(list(HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN = 0.001, HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN = 0.001))

# Determine combined (0h & 12h) NMD_status - if unsure, use highest tx_P_sites_sum_pN value
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined <-  HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_unique %>% 
  mutate(ORFquant_NMD_tx_status_combined = case_when(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status == HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     is.na(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status) ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     is.na(HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status) ~ HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status,
                                                     HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN >= HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN ~ HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status,
                                                     TRUE ~ HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status)) %>% 
  mutate(ORFquant_NMD_tx_status_source = case_when(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status == HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status ~ "ORFquant_HCT_N_AID_UPF1_0h&12h",
                                                   is.na(HCT_N_AID_UPF1_0h_ORFquant_NMD_tx_status) ~ "ORFquant_HCT_N_AID_UPF1_12h",
                                                   is.na(HCT_N_AID_UPF1_12h_ORFquant_NMD_tx_status) ~ "ORFquant_HCT_N_AID_UPF1_0h",
                                                   HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN >= HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN ~ "ORFquant_HCT_N_AID_UPF1_12h",
                                                   TRUE ~ "ORFquant_HCT_N_AID_UPF1_0h"))

# How many transcripts with which status?
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined_number <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined %>% 
  dplyr::count(ORFquant_NMD_tx_status_combined,ORFquant_NMD_tx_status_source)

# Recover most relevant ORF statistics
HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined_txORF <- HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined %>% 
  filter(ORFquant_NMD_tx_status_source %in% c("ORFquant_HCT_N_AID_UPF1_0h&12h", "ORFquant_HCT_N_AID_UPF1_12h")) %>% 
  left_join(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF %>% 
              filter(condition=="HCT_N_AID_UPF1_12h") %>% 
              dplyr::select(gene_id,transcript_id,gene_name,transcript_name, ORF_id_tr, ORF_region, strand, P_sites, Distance_to_lastExEx_compatible, NMD_status_ORFquant, NMD_reason, tx_ORF_pct_P_sites_pN, NMD_tx_ORF_pct_P_sites_pN) %>% 
              mutate(NMD_status_ORFquant_num = case_when(NMD_status_ORFquant == TRUE ~ 1,
                                                         NMD_status_ORFquant == FALSE ~ 0)) %>% 
              arrange(NMD_status_ORFquant_num,tx_ORF_pct_P_sites_pN) %>% 
              group_by(transcript_id) %>% 
              slice_max(tibble(NMD_status_ORFquant_num, tx_ORF_pct_P_sites_pN), n = 1, with_ties = FALSE) %>% 
              ungroup()) %>% 
  bind_rows(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined %>% 
              filter(ORFquant_NMD_tx_status_source %in% c("ORFquant_HCT_N_AID_UPF1_0h")) %>% 
              left_join(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF %>% 
                          filter(condition=="HCT_N_AID_UPF1_0h") %>% 
                          dplyr::select(gene_id,transcript_id,gene_name,transcript_name, ORF_id_tr, ORF_region, strand, P_sites, Distance_to_lastExEx_compatible, NMD_status_ORFquant, NMD_reason, tx_ORF_pct_P_sites_pN, NMD_tx_ORF_pct_P_sites_pN) %>% 
                          mutate(NMD_status_ORFquant_num = case_when(NMD_status_ORFquant == TRUE ~ 1,
                                                                     NMD_status_ORFquant == FALSE ~ 0)) %>% 
                          arrange(NMD_status_ORFquant_num,tx_ORF_pct_P_sites_pN) %>% 
                          group_by(transcript_id) %>% 
                          slice_max(tibble(NMD_status_ORFquant_num, tx_ORF_pct_P_sites_pN), n = 1, with_ties = FALSE) %>% 
                          ungroup()))

# Combine with initial transcriptome
NMDRHT.v1.1_tbl_ORFquant <- NMDRHT.v1.1_tbl %>% 
  left_join(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_combined_txORF %>% 
              dplyr::select(gene_id,transcript_id,gene_name,transcript_name,ORF_id_tr, ORF_region, strand, P_sites, Distance_to_lastExEx_compatible, ORFquant_NMD_tx_status_combined, ORFquant_NMD_tx_status_source, NMD_reason, tx_ORF_pct_P_sites_pN, NMD_tx_ORF_pct_P_sites_pN, HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN, HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN )) %>% 
  mutate(ORFquant_NMD_tx_status_combined = case_when(!is.na(ORFquant_NMD_tx_status_combined) ~ ORFquant_NMD_tx_status_combined,
                                                     NMD_status == TRUE ~ "predicted_NMD",
                                                     NMD_status == FALSE ~ "predicted_coding",
                                                     TRUE ~ "lncRNA")) %>% 
  mutate(log2FC_tx_ORF_pct_P_sites_pN = log2(HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN/HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN)) %>% 
  dplyr::select(-c(HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN,
                   HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN)) %>% 
  mutate(ORFquant_NMD_tx_status_source = case_when(!is.na(ORFquant_NMD_tx_status_source) ~ ORFquant_NMD_tx_status_source,
                                                   !is.na(ORFquant_NMD_tx_status_combined) ~ "ORFanage",
                                                   TRUE ~ "lncRNA")) %>% 
  mutate(NMD_reason = case_when(!is.na(NMD_reason) ~ NMD_reason,
                                ORFquant_NMD_tx_status_combined == "predicted_NMD" & ORF_type == "exact_match" ~ "AS_NMD_UTR3",
                                ORFquant_NMD_tx_status_combined == "predicted_NMD" & ORF_type == "start_match" ~ "AS_NMD",
                                ORFquant_NMD_tx_status_combined == "predicted_coding" ~ "none",
                                ORFquant_NMD_tx_status_combined == "lncRNA" ~ "lncRNA",
                                TRUE ~ "other"))

# export
write_csv(NMDRHT.v1.1_tbl_ORFquant,
          file="Resources/NMDRHT/NMDRHT.v1.1_tbl_ORFquant.csv")

# read data - if necessary
NMDRHT.v1.1_tbl_ORFquant <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_tbl_ORFquant.csv")


## Stats -------------------------------------------------------------------

# How many transcripts with HQ-ORF support from Ribo-Seq?
NMDRHT.v1.1_tbl_ORFquant %>% 
  filter(ORFquant_NMD_tx_status_combined %in% c("coding",
                                                "mixed",
                                                "NMD")) %>% 
  dplyr::count()

# How many transcripts with ORFanage prediction?
NMDRHT.v1.1_tbl_ORFquant %>% 
  filter(ORFquant_NMD_tx_status_combined %in% c("predicted_coding",
                                                "predicted_NMD")) %>% 
  dplyr::count()

#
##
###
# Match ORFs & transcript-properties --------------------------------------
###
##
#

# Objective: Using the best ORF per transcript (favoring potentially NMD-activating ORFs) - generate ORF-annotated GTF file
# Two stages: 1) intersect tx-wise ORF information with NMDRHT annotation -> obtain Ribo-Seq-supported ORFs
#             2) Add ORFanage-based ORFs (not found in Ribo-Seq, but still contain predicted ORF information)

library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(tidygenomics)

# Read in NMDRHT ORFquant file
NMDRHT.v1.1_tbl_ORFquant <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_tbl_ORFquant.csv")

# Prepare NMDRHT ORFquant data for matching
# Fix missing stop codon positions
NMDRHT.v1.1_tbl_ORFquant_forMatch <- NMDRHT.v1.1_tbl_ORFquant %>% 
  filter(!is.na(ORF_region))  %>% 
  separate(col = ORF_region, into = c("chr", "start", "end")) %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  # mutate(start = as.numeric(start),
  #        end = as.numeric(end)) %>% 
  mutate(start = case_when(strand == "+" ~ as.numeric(start),
                           strand == "-" ~ as.numeric(start)-3),
         end = case_when(strand == "+" ~ as.numeric(end)+3,
                         strand == "-" ~ as.numeric(end)))

# read in GTF file
NMDRHT_GTF <- rtracklayer::readGFF("Resources/NMDRHT/NMDRHT.v1.1.sort.gtf")

# Prepare NMDRHT GTF data for matching
NMDRHT_GTF_forMatch <- NMDRHT_GTF %>% 
  filter(type == "exon") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqid") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect with tidygenomics
NMDRHT_intersect <- genome_intersect(NMDRHT.v1.1_tbl_ORFquant_forMatch, NMDRHT_GTF_forMatch, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
NMDRHT_intersect_filt <- NMDRHT_intersect %>% 
  filter(transcript_id.x == transcript_id.y & strand.x == strand.y) %>% 
  dplyr::select(-c(transcript_id.y, strand.y)) %>% 
  dplyr::rename("transcript_id" = "transcript_id.x",
                "strand" = "strand.x") %>% 
  mutate(source="ORFquant",
         type="CDS",
         score=NA,
         phase=NA) %>% 
  relocate(chr,
           source,
           type,
           start,
           end,
           score,
           strand,
           phase,
           transcript_id
  ) %>% 
  left_join(NMDRHT.v1.1_tbl_ORFquant %>% 
              filter(!is.na(ORF_region)))

# Prepare NMDRHT-ORFquant supplemented annotation
NMDRHT_GTF_ORFquant_supplemented <- NMDRHT_GTF %>% 
  filter(transcript_id %in% NMDRHT_intersect_filt$transcript_id) %>% 
  filter(type != "CDS") %>% 
  dplyr::rename("chr" = "seqid") %>% 
  mutate(UIC_total_support = as.numeric(UIC_total_support),
         NMD_status = as.logical(NMD_status),
         stop_to_lastEJ = as.numeric(stop_to_lastEJ),
         num_of_downEJs = as.numeric(num_of_downEJs),
         UTR3_length = as.numeric(UTR3_length)) %>% 
  bind_rows(NMDRHT_intersect_filt) %>%  
  bind_rows(NMDRHT_GTF %>% 
              filter(!transcript_id %in% NMDRHT_intersect_filt$transcript_id) %>% 
              dplyr::rename("chr" = "seqid") %>% 
              mutate(UIC_total_support = as.numeric(UIC_total_support),
                     NMD_status = as.logical(NMD_status),
                     stop_to_lastEJ = as.numeric(stop_to_lastEJ),
                     num_of_downEJs = as.numeric(num_of_downEJs),
                     UTR3_length = as.numeric(UTR3_length))) %>% 
  arrange(transcript_id) %>% 
  relocate(GENCODE_gene_id,
           GENCODE_transcript_id,
           GENCODE_gene_name,
           GENCODE_transcript_name,
           GENCODE_gene_type,
           GENCODE_transcript_type,
           GENCODE_transcript_support_level,
           ensembl_canonical,
           GTF_transcript_id,
           NMDRHT_locus_id,
           .after=transcript_biotype) %>% 
  dplyr::rename("ORFanage_status"="orfanage_status",
                "ORFanage_duplicity"="orfanage_duplicity",
                "ORFanage_template"="orfanage_template",
                "ORFanage_template_source"="orfanage_template_source",
                "ORFanage_ORF_type"="ORF_type",
                "ORFanage_start_codon_type"="start_codon_type",
                "ORFanage_stop_codon_type"="stop_codon_type",
                "ORFanage_ORF_id"="ORF_id",
                "ORFanage_NMD_status"="NMD_status",
                "ORFanage_stop_to_lastEJ"="stop_to_lastEJ",
                "ORFanage_num_of_downEJs"="num_of_downEJs",
                "ORFanage_UTR3_length"="UTR3_length") %>% 
  dplyr::rename("ORFquant_ORF_id_tr"="ORF_id_tr",
                "ORFquant_ORF_region"="ORF_region",
                "ORFquant_P_sites"="P_sites",
                "ORFquant_Distance_to_lastExEx"="Distance_to_lastExEx_compatible",
                "ORFquant_NMD_reason"="NMD_reason",
                "ORFquant_tx_ORF_pct_P_sites_pN"="tx_ORF_pct_P_sites_pN",
                "ORFquant_NMD_tx_ORF_pct_P_sites_pN"="NMD_tx_ORF_pct_P_sites_pN",
                "ORFquant_log2FC_tx_ORF_pct_P_sites_pN"="log2FC_tx_ORF_pct_P_sites_pN") %>% 
  relocate(ORFanage_status,
           ORFanage_duplicity,
           ORFanage_template,
           ORFanage_template_source,
           ORFanage_ORF_type,
           ORFanage_start_codon_type,
           ORFanage_stop_codon_type,
           ORFanage_ORF_id,
           ORFanage_NMD_status,
           ORFanage_stop_to_lastEJ,
           ORFanage_num_of_downEJs,
           ORFanage_UTR3_length,
           .after = last_col())

# Check total number of transcripts
NMDRHT_GTF_ORFquant_supplemented %>% distinct(transcript_id) %>% dplyr::count()

##
# Generate new NMDRHT v1.2 tibble (new version 1.2 = contains Ribo-Seq-derived information)
##
NMDRHT.v1.2_tbl <- NMDRHT_GTF_ORFquant_supplemented %>% 
  filter(type == "transcript") %>% 
  mutate(coordinates = paste0(chr,":",start,"-",end,":",strand),.after=NMDRHT_locus_id) %>% 
  dplyr::select(-c(chr,source,type,start,end,score,strand,phase,starts_with(c("ORFquant")))) %>% 
  left_join(NMDRHT.v1.1_tbl_ORFquant %>% 
              dplyr::select(transcript_id,
                            ORFquant_NMD_tx_status_combined,
                            ORFquant_NMD_tx_status_source,
                            NMD_reason) %>% 
              dplyr::rename("NMD_tx_status"="ORFquant_NMD_tx_status_combined",
                            "NMD_tx_source"="ORFquant_NMD_tx_status_source",
                            "NMD_tx_reason"="NMD_reason")) %>% 
  left_join(NMDRHT_GTF_ORFquant_supplemented %>% 
              filter(type == "CDS") %>%
              distinct(transcript_id, .keep_all=TRUE) %>% 
              dplyr::select(transcript_id, starts_with(c("ORFquant")))) %>% 
  relocate(ORFanage_status,
           ORFanage_duplicity,
           ORFanage_template,
           ORFanage_template_source,
           ORFanage_ORF_type,
           ORFanage_start_codon_type,
           ORFanage_stop_codon_type,
           ORFanage_ORF_id,
           ORFanage_NMD_status,
           ORFanage_stop_to_lastEJ,
           ORFanage_num_of_downEJs,
           ORFanage_UTR3_length,
           .after = last_col())

##
#### save information ----------------------------------------------------------
##
NMDRHT.v1.2_tbl %>% write_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl.csv")

NMDRHT.v1.2_tbl <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl.csv")

# Generate final ORF-supplemented GTF - collapse all information on transcript_level and remove unneccesary information
NMDRHT_GTF_ORFquant_supplemented_final <- NMDRHT_GTF_ORFquant_supplemented %>% 
  dplyr::select(chr,source,type,start,end,score,strand,phase, gene_id, transcript_id) %>% 
  left_join(NMDRHT.v1.2_tbl) %>% 
  mutate(across(.cols = c(GENCODE_gene_id,
                          GENCODE_transcript_id,
                          GENCODE_gene_name,
                          GENCODE_transcript_name,
                          GENCODE_gene_type,
                          GENCODE_transcript_type,
                          GENCODE_transcript_support_level,
                          NMDRHT_locus_id,
                          coordinates,
                          gene_id_conflict,
                          structural_category,
                          structural_category_simple,
                          subcategory,
                          class_code,
                          class_code_simple,
                          UIC_total_support,
                          UIC_LR_support,
                          UIC_SR_support,
                          structural_category,
                          subcategory,
                          ensembl_canonical,
                          within_CAGE_peak,
                          within_polyA_site,
                          polyA_motif_found,
                          study,
                          GTF_transcript_id,
                          cell_line,
                          method,
                          sequencing,
                          annotation,
                          NMD_tx_status,
                          NMD_tx_source,
                          NMD_tx_reason,
                          starts_with(c("ORFquant", "ORFanage"))
  ),
  .fns = ~case_when(type %in% c("exon", "CDS") ~ NA,
                    TRUE ~ (.))))

# Export annotation
# *Note* v1.2 due to ORF annotation
rtracklayer::export(NMDRHT_GTF_ORFquant_supplemented_final,
                    "Resources/NMDRHT/NMDRHT.v1.2.gtf")

# Run sort and index
command_sort <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.2.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), sep="", collapse = "")
output_sort <- system(command_sort, intern = TRUE)

command_index <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                         "index ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), sep="", collapse = "")
output_index <- system(command_index, intern = TRUE)

#
##
# FactR -------------------------------------------------------------------
##
#
# Use to predict NMD sensitivity of NMDRHT_v1.2 transcripts

library(factR) 

# Import
NMDRHT_v1.2.sort.gtf_for_factR <- factR::importGTF("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf") 

# Predict NMD
NMDRHT_v1.2.sort.gtf_for_factR_NMDprediction.out <- predictNMD(NMDRHT_v1.2.sort.gtf_for_factR)  

# Join, rename and relocate with NMDRHT.v1.2_tbl
NMDRHT.v1.2_tbl_NMD <- NMDRHT.v1.2_tbl %>% 
  left_join(NMDRHT_v1.2.sort.gtf_for_factR_NMDprediction.out %>% 
              dplyr::rename("transcript_id" = "transcript",
                            "UTR3_length" = "3'UTR_length",
                            "NMD_50nt_rule" = "is_NMD")) %>% 
  relocate("NMD_50nt_rule",
           "stop_to_lastEJ",
           "num_of_downEJs",
           "UTR3_length", .after=NMD_tx_reason)


# How many NMD-sensitive (50-nt rule) transcripts
NMDRHT.v1.2_tbl_NMD %>% 
  dplyr::count(NMD_tx_status,NMD_50nt_rule) %>% 
  arrange(desc(n))

# Check for problematic candidates (factR and ORFquant do not agree?)
# should be empty
NMDRHT.v1.2_tbl_NMD %>% 
  filter(NMD_50nt_rule == FALSE & NMD_tx_status %in% c("NMD", "mixed")) 

# export
write_csv(NMDRHT.v1.2_tbl_NMD,
          file="Resources/NMDRHT/NMDRHT.v1.2_tbl_NMD.csv")

NMDRHT.v1.2_tbl_NMD <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl_NMD.csv")

## factR - CDS sequence ---------------------------------------------------------

# Codes from factR
NMDRHT_v1.2.sort.gtf_for_factR_cleaned <- .extractCDSchecks(NMDRHT_v1.2.sort.gtf_for_factR, Hsapiens)

NMDRHT_v1.2.sort.gtf_for_factR_cleaned_AA <- .getSequence(NMDRHT_v1.2.sort.gtf_for_factR_cleaned, Hsapiens) %>% 
  filter(noATG == FALSE) %>% 
  mutate(x = stringr::str_replace(x, '\\*', ''))

NMDRHT_v1.2.sort.gtf_for_factR_cleaned_AA_fasta <- c(rbind(paste0(">", NMDRHT_v1.2.sort.gtf_for_factR_cleaned_AA$id), NMDRHT_v1.2.sort.gtf_for_factR_cleaned_AA$x))

write(x = NMDRHT_v1.2.sort.gtf_for_factR_cleaned_AA_fasta, file = "Resources/NMDRHT/NMDRHT.v1.2.fa")  

###
## Get transcript and 3'UTR lengths and sequences  ----------------------------------------------------------
###

# Import GTF file
NMDRHT.v1.2.sort.gtf <- rtracklayer::import(file.path("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), format="gtf")

# Generate TxDb from annotation
NMDRHT.v1.2_txdb1 <- makeTxDbFromGRanges(NMDRHT.v1.2.sort.gtf, drop.stop.codons=FALSE)

# Get transcript lengths with 5'UTR, CDS and 3'UTR lengths individually
NMDRHT_tx_length <- transcriptLengths(NMDRHT.v1.2_txdb1, with.cds_len=TRUE,
                                      with.utr5_len=TRUE, with.utr3_len=TRUE)

# Remove tx_id and gene_id columns as they are redundant and not needed
NMDRHT_tx_length_mod <- tibble::as_tibble(NMDRHT_tx_length) %>% 
  dplyr::select(-c(tx_id, gene_id))

# Combine the NMDRegHumanTxome.v1.1_tbl table (from Main analysis) with transcripts lengths
NMDRHT.v1.2_tbl_length <- left_join(NMDRHT.v1.2_tbl, 
                                              NMDRHT_tx_length_mod, by = c("transcript_id" = "tx_name"))

# Get genome fasta file - Required to obtain nucleotide sequence
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

## Extract the exon ranges grouped by transcript from 'txdb':
NMDRHT_transcripts <- exonsBy(NMDRHT.v1.2_txdb1, by="tx", use.names=TRUE)

# Get 5'UTRs from NMDRHT annotation and generate relevant data frame
NMDRHT_utr5 <- fiveUTRsByTranscript(NMDRHT.v1.2_txdb1, use.names=TRUE)

# Get CDS from NMDRHT annotation and generate relevant data frame
NMDRHT_cds <- cdsBy(NMDRHT.v1.2_txdb1, by="tx", use.names=TRUE)

# Get 3'UTRs from NMDRHT annotation and generate relevant data frame
NMDRHT_utr3 <- threeUTRsByTranscript(NMDRHT.v1.2_txdb1, use.names=TRUE)

# Get sequences of transcripts, CDS, 5'UTRs and 3'UTRs
NMDRHT_transcripts_seqs <- extractTranscriptSeqs(Hsapiens, NMDRHT_transcripts)
NMDRHT_utr5_seqs <- extractTranscriptSeqs(Hsapiens, NMDRHT_utr5)
NMDRHT_cds_seqs <- extractTranscriptSeqs(Hsapiens, NMDRHT_cds)
NMDRHT_utr3_seqs <- extractTranscriptSeqs(Hsapiens, NMDRHT_utr3)

# Only first 200 nucleotides of 3' UTR
NMDRHT_utr3_seqs_all <- extractTranscriptSeqs(Hsapiens, NMDRHT_utr3) 

NMDRHT_utr3_seqs_over200 <- NMDRHT_utr3_seqs_all[width(NMDRHT_utr3_seqs_all) >= 200]
NMDRHT_utr3_seqs_under200 <- NMDRHT_utr3_seqs_all[width(NMDRHT_utr3_seqs_all) < 200]

NMDRHT_utr3_seqs_200 <- subseq(NMDRHT_utr3_seqs_over200, start=1,width=200)

NMDRHT_utr3_seqs_200_final <- c(NMDRHT_utr3_seqs_200,
                                NMDRHT_utr3_seqs_under200)

####
### GC-content  --------------------------------------------------------
####

# Calculate GC content for transcripts, CDS, 5'UTRs and 3'UTRss
NMDRHT_transcript_GC <- as.numeric(Biostrings::letterFrequency(x = NMDRHT_transcripts_seqs, letters = "GC", as.prob = TRUE))
names(NMDRHT_transcript_GC) <- names(NMDRHT_transcripts_seqs)
NMDRHT_transcript_GC <- enframe(NMDRHT_transcript_GC, name = "tx_id", value = "GC_tx")

NMDRHT_utr5_GC <- as.numeric(Biostrings::letterFrequency(x = NMDRHT_utr5_seqs, letters = "GC", as.prob = TRUE))
names(NMDRHT_utr5_GC) <- names(NMDRHT_utr5_seqs)
NMDRHT_utr5_GC <- enframe(NMDRHT_utr5_GC, name = "tx_id", value = "GC_utr5")

NMDRHT_cds_GC <- as.numeric(Biostrings::letterFrequency(x = NMDRHT_cds_seqs, letters = "GC", as.prob = TRUE))
names(NMDRHT_cds_GC) <- names(NMDRHT_cds_seqs)
NMDRHT_cds_GC <- enframe(NMDRHT_cds_GC, name = "tx_id", value = "GC_cds")

NMDRHT_utr3_GC <- as.numeric(Biostrings::letterFrequency(x = NMDRHT_utr3_seqs, letters = "GC", as.prob = TRUE))
names(NMDRHT_utr3_GC) <- names(NMDRHT_utr3_seqs)
NMDRHT_utr3_GC <- enframe(NMDRHT_utr3_GC, name = "tx_id", value = "GC_utr3")

NMDRHT_utr3_GC_200nt <- as.numeric(Biostrings::letterFrequency(x = NMDRHT_utr3_seqs_200_final, letters = "GC", as.prob = TRUE))
names(NMDRHT_utr3_GC_200nt) <- names(NMDRHT_utr3_seqs_200_final)
NMDRHT_utr3_GC_200nt <- enframe(NMDRHT_utr3_GC_200nt, name = "tx_id", value = "GC_utr3_200nt")


# Save as csv
NMDRHT_utr3_GC_200nt %>% write_csv("Resources/NMDRHT/NMDRHT_utr3_GC_200nt.csv")

# Combine the NMDRegHumanTxome.v1.1_tbl_MasterTable_length table with GC contents
NMDRHT.v1.2_tbl_length_GC <- left_join(NMDRHT.v1.2_tbl_length, 
                                                             NMDRHT_transcript_GC,
                                                             by = c("transcript_id" = "tx_id")) %>% 
  left_join(NMDRHT_utr5_GC,
            by = c("transcript_id" = "tx_id")) %>% 
  left_join(NMDRHT_cds_GC,
            by = c("transcript_id" = "tx_id")) %>% 
  left_join(NMDRHT_utr3_GC,
            by = c("transcript_id" = "tx_id"))

#
##
# Structure prediction ----------------------------------------------------
##
#

# Write 3'UTRs as fasta file
writeXStringSet(NMDRHT_utr3_seqs, file.path("Resources/NMDRHT/NMDRHT_ORFquant_NMD.v1.2.utr3_seqs.fa"), format="fasta")

##TEST
# Write 3'UTRs as fasta file
# writeXStringSet(gencode_utr3_seqs["ENST00000641515.2"], file.path(ref_dir, "Gencode", "gencode.v42.SIRVomeERCCome.utr3_seqs_TEST.fa"), format="fasta")

# Run RNAfold (note: this took about 8 hours for the full gencode 3UTR annotation!!!)
command <- paste(c("RNAfold ", "--noPS --jobs='0'", " --infile='", file.path("/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT_ORFquant_NMD.v1.2.utr3_seqs.fa"), "'"), sep="", collapse = "")
output <- system(command, intern = TRUE)

# Integrate output into data frame

UTR3fold <- tibble("transcript_id" = output[seq(1, length(output), 3)],
                   "sequence" = output[seq(2, length(output), 3)],
                   "fold" = output[seq(3, length(output), 3)]) %>% 
  separate("transcript_id", c(NA,"transcript_id"), ">")

UTR3fold$fold2 <- gsub("\\.(?=[^.]*\\.)", "", UTR3fold$fold, perl=TRUE)

UTR3fold_final <- UTR3fold %>% 
  mutate(utr3_mfe = parse_number(fold2)) %>% 
  dplyr::select(-fold2) %>% 
  mutate(utr3_length_str = str_length(sequence)) %>% 
  mutate(utr3_mfe_nt = utr3_mfe/utr3_length_str)

# Save output to file (to not lose it!)
write_csv(UTR3fold_final, file=file.path("Resources/NMDRHT/NMDRHT_ORFquant_NMD.v1.2.utr3_fold.csv"))

# In case needed, load UTR3fold file again
UTR3fold_final <- read_csv("Resources/NMDRHT/NMDRHT_ORFquant_NMD.v1.2.utr3_fold.csv") 

# Combine with other characteristics
NMDRHT.v1.2_tbl_length_GC_mfe <-  NMDRHT.v1.2_tbl_length_GC %>% 
  left_join(UTR3fold_final %>% dplyr::select(transcript_id,utr3_mfe, utr3_length_str,utr3_mfe_nt))

# Safe final data frame 
NMDRHT.v1.2_tbl_length_GC_mfe %>% write_csv(file.path("Resources/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe.csv"))

# Read in if necessary
NMDRHT.v1.2_tbl_length_GC_mfe <- read_csv(file.path("Resources/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe.csv"))

## Stats -------------------------------------------------------------------

# How many transcripts with HQ-ORF support from Ribo-Seq?
NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  filter(NMD_tx_status %in% c("coding",
                              "mixed",
                              "NMD")) %>% 
  dplyr::count()

# How many transcripts with ORFanage prediction?
NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  filter(NMD_tx_status %in% c("predicted_coding",
                              "predicted_NMD")) %>% 
  dplyr::count()

# How many ORFs with which NMD_tx_status?
NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  dplyr::count(NMD_tx_status) %>% 
  mutate(n_per = n/sum(n))

# How many ORFs with which NMD_tx_status?
NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  filter(NMD_tx_status %in% c("coding",
                              "mixed",
                              "NMD")) %>% 
  dplyr::count(NMD_tx_status) %>% 
  mutate(n_per = n/sum(n))

#
##
# Check splice junc - Splam -----------------------------------------------
##
#

# https://github.com/Kuanhao-Chao/splam
# PMID: 39285451
# Run on Gehring-Server2024 with CUDA
# Standard pipeline


## Import splam junction_score ---------------------------------------------
NMDRHT_v1.2_splam_junction_score <- read_delim("Resources/splam/NMDRHT/junction_score.bed", 
                                               delim = "\t", escape_double = FALSE, 
                                               col_names = c("chr", "start", "end", "junc_name", "intron_num", "stramd", "donor_score", "acceptor_score", "transcript_id"), trim_ws = TRUE)  %>%
  separate_rows(transcript_id, sep = ",")

# Calculate average donor and acceptor score for each transcript
NMDRHT_v1.2_splam_junction_score_avg <- NMDRHT_v1.2_splam_junction_score %>% 
  group_by(transcript_id) %>% 
  mutate(avg_donor_score = mean(donor_score, na.rm = TRUE),
         avg_acceptor_score = mean(acceptor_score, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct(transcript_id, avg_donor_score, avg_acceptor_score)

# Combine with NMDRHT_Global_Tx_tbl (import from main)
NMDRHT.v1.2_tbl_length_GC_mfe_splam <- NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  left_join(NMDRHT_v1.2_splam_junction_score_avg)

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>% write_csv(file.path("Resources/splam/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe_splam.csv"))

# Check for missing transcripts - two explanations: 1) single exon transcripts (n=3267); 2) transcripts from overlapping genes (n=128) [manually detected; must be a problem when assigning the introns back to transcripts]
NMDRHT.v1.2_tbl_length_GC_mfe_splam %>% 
  filter(is.na(avg_donor_score) & nexon == 1) %>% 
  dplyr::count()

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>% 
  filter(is.na(avg_donor_score) & nexon > 1) %>% 
  dplyr::count()

##
# Match ORFs for complete ORFquant (used for proteomics analysis) ----------------------------------------
##

# Prepare NMDRHT ORFquant data for matching
# Fix missing stop codon positions
HCT_N_AID_UPF1_ORFquant_forMatch <- HCT_N_AID_UPF1_ORFquant %>% 
  filter(!is.na(ORF_region))  %>% 
  separate(col = ORF_region, into = c("chr", "start", "end")) %>% 
  dplyr::select(ORF_id_tr, transcript_id, chr, start, end, strand) %>% 
  distinct() %>% 
  # mutate(start = as.numeric(start),
  #        end = as.numeric(end)) %>% 
  mutate(start = case_when(strand == "+" ~ as.numeric(start),
                           strand == "-" ~ as.numeric(start)-3),
         end = case_when(strand == "+" ~ as.numeric(end)+3,
                         strand == "-" ~ as.numeric(end)))

# read in GTF file
NMDRHT_GTF <- rtracklayer::readGFF("Resources/NMDRHT/NMDRHT.v1.1.sort.gtf")

# Prepare NMDRHT GTF data for matching
NMDRHT_GTF_forMatch <- NMDRHT_GTF %>% 
  filter(type == "exon") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqid") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect
HCT_N_AID_UPF1_ORFquant_intersect <- genome_intersect(NMDRHT_GTF_forMatch, HCT_N_AID_UPF1_ORFquant_forMatch, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
HCT_N_AID_UPF1_ORFquant_intersect_filt <- HCT_N_AID_UPF1_ORFquant_intersect %>% 
  filter(transcript_id.x == transcript_id.y & strand.x == strand.y) %>% 
  dplyr::select(-c(transcript_id.y, strand.y)) %>% 
  dplyr::rename("transcript_id" = "transcript_id.x",
                "strand" = "strand.x") %>% 
  mutate(source="ORFquant",
         type="CDS",
         score=NA,
         phase=NA) %>% 
  relocate(chr,
           source,
           type,
           start,
           end,
           score,
           strand,
           phase,
           transcript_id
  ) %>% 
  left_join(HCT_N_AID_UPF1_ORFquant %>% 
              dplyr::select(-c(source, type, start, end, score, strand, phase)) %>% 
              filter(!is.na(ORF_region))) %>% 
  mutate(transcript_id_ORF = paste0(transcript_id,"_",ORF_id_tr))

# Create table with transcript_id_original and new ORF-containing IDs
HCT_N_AID_UPF1_ORFquant_intersect_filt_distinct <- HCT_N_AID_UPF1_ORFquant_intersect_filt %>% 
  distinct(transcript_id_ORF, transcript_id)

# Prepare NMDRHT-ORFquant supplemented annotation
HCT_N_AID_UPF1_ORFquant_supplemented <- NMDRHT_GTF %>% 
  filter(transcript_id %in% HCT_N_AID_UPF1_ORFquant_intersect_filt$transcript_id) %>% 
  left_join(HCT_N_AID_UPF1_ORFquant_intersect_filt_distinct,
            relationship = "many-to-many",
            by=c("transcript_id"="transcript_id")) %>% 
  filter(type != "CDS") %>% 
  dplyr::rename("chr" = "seqid") %>% 
  # mutate(UIC_total_support = as.numeric(UIC_total_support),
  #        NMD_status = as.logical(NMD_status),
  #        stop_to_lastEJ = as.numeric(stop_to_lastEJ),
  #        num_of_downEJs = as.numeric(num_of_downEJs),
  #        UTR3_length = as.numeric(UTR3_length)) %>% 
  bind_rows(HCT_N_AID_UPF1_ORFquant_intersect_filt) %>%  
  arrange(transcript_id, ORF_id_tr) %>% 
  relocate(GENCODE_gene_id,
           GENCODE_transcript_id,
           GENCODE_gene_name,
           GENCODE_transcript_name,
           GENCODE_gene_type,
           GENCODE_transcript_type,
           GENCODE_transcript_support_level,
           ensembl_canonical,
           GTF_transcript_id,
           NMDRHT_locus_id,
           .after=transcript_biotype) %>% 
  dplyr::rename("transcript_id_NMDRHT" = "transcript_id",
                "transcript_id" = "transcript_id_ORF",
                "ORFanage_status"="orfanage_status",
                "ORFanage_duplicity"="orfanage_duplicity",
                "ORFanage_template"="orfanage_template",
                "ORFanage_template_source"="orfanage_template_source",
                "ORFanage_ORF_type"="ORF_type",
                "ORFanage_start_codon_type"="start_codon_type",
                "ORFanage_stop_codon_type"="stop_codon_type",
                "ORFanage_ORF_id"="ORF_id",
                "ORFanage_NMD_status"="NMD_status",
                "ORFanage_stop_to_lastEJ"="stop_to_lastEJ",
                "ORFanage_num_of_downEJs"="num_of_downEJs",
                "ORFanage_UTR3_length"="UTR3_length") %>% 
  dplyr::rename("ORFquant_ORF_id_tr"="ORF_id_tr",
                "ORFquant_ORF_region"="ORF_region",
                "ORFquant_P_sites"="P_sites",
                "ORFquant_Distance_to_lastExEx"="Distance_to_lastExEx_compatible") %>% 
  relocate(ORFanage_status,
           ORFanage_duplicity,
           ORFanage_template,
           ORFanage_template_source,
           ORFanage_ORF_type,
           ORFanage_start_codon_type,
           ORFanage_stop_codon_type,
           ORFanage_ORF_id,
           ORFanage_NMD_status,
           ORFanage_stop_to_lastEJ,
           ORFanage_num_of_downEJs,
           ORFanage_UTR3_length,
           .after = last_col())

HCT_N_AID_UPF1_ORFquant_supplemented %>%
  filter(type == "transcript") %>% 
  distinct(transcript_id, ORFquant_ORF_id_tr) %>%
  dplyr::count()

# Check total number of transcripts
HCT_N_AID_UPF1_ORFquant_supplemented %>% distinct(ORFquant_ORF_id_tr, transcript_id) %>% dplyr::count()

# Generate new tibble
HCT_N_AID_UPF1_ORFquant_supplemented_tbl <- HCT_N_AID_UPF1_ORFquant_supplemented %>% 
  filter(type == "transcript") %>% 
  mutate(coordinates = paste0(chr,":",start,"-",end,":",strand),.after=NMDRHT_locus_id) %>% 
  dplyr::select(-c(chr,source,type,start,end,score,strand,phase,starts_with(c("ORFquant")))) %>% 
  left_join(HCT_N_AID_UPF1_ORFquant_supplemented %>% 
              filter(type == "CDS") %>%
              distinct(transcript_id, .keep_all=TRUE) %>% 
              dplyr::select(transcript_id, starts_with(c("ORFquant")))) %>% 
  relocate(ORFanage_status,
           ORFanage_duplicity,
           ORFanage_template,
           ORFanage_template_source,
           ORFanage_ORF_type,
           ORFanage_start_codon_type,
           ORFanage_stop_codon_type,
           ORFanage_ORF_id,
           ORFanage_NMD_status,
           ORFanage_stop_to_lastEJ,
           ORFanage_num_of_downEJs,
           ORFanage_UTR3_length,
           .after = last_col())

# Generate final ORF-supplemented GTF - collapse all information on transcript_level and remove 
HCT_N_AID_UPF1_ORFquant_supplemented_final <- HCT_N_AID_UPF1_ORFquant_supplemented %>% 
  dplyr::select(chr,source,type,start,end,score,strand,phase, gene_id, gene_name, transcript_id, ORFquant_ORF_id_tr, transcript_id_NMDRHT) %>% 
  mutate(name = paste0(gene_name, ";", transcript_id)) %>% 
  relocate(name, .before = gene_id)

# Export annotation
rtracklayer::export.gff(HCT_N_AID_UPF1_ORFquant_supplemented_final,
                    "Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_supplemented_final.gtf")

# Run sort and index
command_sort <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_supplemented_final.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_supplemented_final.sort.gtf"), sep="", collapse = "")
output_sort <- system(command_sort, intern = TRUE)

command_index <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                         "index ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_supplemented_final.sort.gtf"), sep="", collapse = "")
output_index <- system(command_index, intern = TRUE)

# Import GTF file
HCT_N_AID_UPF1_ORFquant_supplemented_gtf <- rtracklayer::import(file.path("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_supplemented_final.sort.gtf"), format="gtf")

as_tibble(HCT_N_AID_UPF1_ORFquant_supplemented_gtf) %>% 
  filter(is.na(gene_id))

# Generate TxDb from annotation
HCT_N_AID_UPF1_ORFquant_txdb1 <- makeTxDbFromGRanges(HCT_N_AID_UPF1_ORFquant_supplemented_gtf, drop.stop.codons=FALSE)

# Get genome fasta file - Required to obtain nucleotide sequence
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38

# Get CDS from NMDRHT annotation and generate relevant data frame
HCT_N_AID_UPF1_ORFquant_cds <- cdsBy(HCT_N_AID_UPF1_ORFquant_txdb1, by="tx", use.names=TRUE)

save(HCT_N_AID_UPF1_ORFquant_cds,
     file = paste0("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant_cds.rds"))
