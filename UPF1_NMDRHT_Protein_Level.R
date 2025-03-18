#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Protein_Level
# Objective: Perform protien-level data preparation and initial analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)
library(GenomicFeatures)

##
# NMDRHT protein sequence --------------------------------------------------------
##

# Load required data
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

# Import GTF file
NMDRHT.v1.2.sort.gtf <- rtracklayer::import(file.path("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), format="gtf")

# Generate TxDb from annotation
NMDRHT.v1.2_txdb1 <- GenomicFeatures::makeTxDbFromGRanges(NMDRHT.v1.2.sort.gtf, drop.stop.codons=FALSE)

# Get genome fasta file - Required to obtain nucleotide sequence
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Get CDS from NMDRHT annotation and generate relevant data frame
NMDRHT_cds <- cdsBy(NMDRHT.v1.2_txdb1, by="tx", use.names=TRUE)

# Get sequences of transcripts, CDS
NMDRHT_utr5_seqs <- extractTranscriptSeqs(Hsapiens, NMDRHT_utr5)

## A sanity check:
stopifnot(identical(width(NMDRHT_cds_seqs), unname(sum(width(NMDRHT_cds)))))

# Translate CDS nucleotide sequence to amino acids
NMDRHT_cds_seqs_AA <- translate(NMDRHT_cds_seqs, if.fuzzy.codon="solve")

# create dataframe and join with NMDRHT information
NMDRHT_cds_seqs_AA_df <- as.data.frame(NMDRHT_cds_seqs_AA) %>% 
  rownames_to_column(var="transcript_id") %>% 
  dplyr::rename("sequence" = "x") %>% 
  # remove stop codon from width calculation
  mutate(width = str_length(sequence)-1) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% 
              dplyr::select(gene_id, transcript_id, gene_name, transcript_name, starts_with("GENCODE"), starts_with("NMD"), starts_with("ORFquant"), starts_with("ORFanage"), DTE_cluster))

# Save data as csv
save(NMDRHT_cds_seqs_AA_df,
     file = paste0("Resources/NMDRHT/NMDRHT_cds_seqs_AA_df.rds"))

save(NMDRHT_cds,
     file = paste0("Resources/NMDRHT/NMDRHT_cds.rds"))

# Check for how many stop codons are present
NMDRHT_cds_seqs_AA_df_2 <- NMDRHT_cds_seqs_AA_df %>% 
  mutate(stop_n = str_count(sequence, pattern = "\\*"), .after=sequence) 

# Two transcripts with seemingly no stop codon?
NMDRHT_cds_seqs_AA_df_2 %>% dplyr::count(stop_n)

# Load GENCODE canonical protein information - generated in "UPF1_NMMDRHT_Annotation.R"
load("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds_AA_df.rds")

# Match NMDRHT with GENCODE v42 **canonical** protein sequence - by gene_id!
NMDRHT_cds_seqs_AA_df_GENCODE_match <- NMDRHT_cds_seqs_AA_df %>% 
  filter(NMD_tx_reason == "AS_NMD") %>% 
  left_join(gencode.v42.annotation_canonical_txDb_cds_AA_df %>% 
              dplyr::select(sequence, width, gene_id, gene_id) %>% 
              dplyr::rename("GENCODE_sequence" = "sequence",
                            "GENCODE_width" = "width")) %>% 
  mutate(rel_proteinWidth_AS_NMD = width/GENCODE_width*100, .after=width)

## Save information to rds -------------------------------------------
save(NMDRHT_cds_seqs_AA_df_GENCODE_match,
     file = paste0("Resources/NMDRHT/NMDRHT_cds_seqs_AA_df_GENCODE_match.rds"))

##
# AS-NMD protein domains ---------------------------------------------------------------
##

# Objective: determine difference in encoded protein domains between AS-NMD and canonical transcripts

##
## Pfam --------------------------------------------------
##

# Comment: was run locally with InterProScan - Pfam output imported

# GENCODE data
GENCODE_canonical_interPro_Pfam <- read_delim("Resources/GENCODE/InterProScan/GENCODE_canonical.tsv", 
                                              delim = "\t", escape_double = FALSE, 
                                              col_names = c("transcript_id",
                                                            "seq_md5",
                                                            "seq_length",
                                                            "analysis",
                                                            "signature_id",
                                                            "signature_desc",
                                                            "start_pos",
                                                            "stop_pos",
                                                            "score",
                                                            "status",
                                                            "run_date",
                                                            "interpro_id",
                                                            "interpro_desc",
                                                            "GO_annot",
                                                            "pathways"), trim_ws = TRUE) %>% 
  dplyr::select(-c(GO_annot, pathways)) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "transcript") %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name)) %>% 
  dplyr::relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  group_by(transcript_id) %>% 
  mutate(n_domains_GENCODE = dplyr::n()) %>% 
  ungroup()

# NDMRHT data
NMDRHT_interPro_Pfam <- read_delim("Resources/NMDRHT/InterProScan/NMDRHT.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   col_names = c("transcript_id",
                                                 "seq_md5",
                                                 "seq_length",
                                                 "analysis",
                                                 "signature_id",
                                                 "signature_desc",
                                                 "start_pos",
                                                 "stop_pos",
                                                 "score",
                                                 "status",
                                                 "run_date",
                                                 "interpro_id",
                                                 "interpro_desc",
                                                 "GO_annot",
                                                 "pathways"), trim_ws = TRUE) %>% 
  dplyr::select(-c(GO_annot, pathways)) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name, NMD_tx_status, NMD_tx_reason, DTE_cluster)) %>% 
  dplyr::relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  group_by(transcript_id) %>% 
  mutate(n_domains_NMDRHT = dplyr::n()) %>% 
  ungroup()

##
## MobiDBLite --------------------------------------------------
##

# Comment: was run locally with InterProScan - MobiDBLite output imported

# GENCODE data
GENCODE_canonical_interPro_MobiDBLite <- read_delim("Resources/GENCODE/InterProScan/GENCODE_canonical_MobiDBLite.tsv", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    col_names = c("transcript_id",
                                                                  "seq_md5",
                                                                  "seq_length",
                                                                  "analysis",
                                                                  "signature_id",
                                                                  "signature_desc",
                                                                  "start_pos",
                                                                  "stop_pos",
                                                                  "score",
                                                                  "status",
                                                                  "run_date",
                                                                  "interpro_id",
                                                                  "interpro_desc",
                                                                  "GO_annot",
                                                                  "pathways"), trim_ws = TRUE) %>% 
  dplyr::select(-c(GO_annot, pathways)) %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "transcript") %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name)) %>% 
  dplyr::relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  group_by(transcript_id) %>% 
  mutate(n_unstruct_GENCODE = n()) %>% 
  ungroup()

# NDMRHT data
NMDRHT_interPro_MobiDBLite <- read_delim("Resources/NMDRHT/InterProScan/NMDRHT_MobiDBLite.tsv", 
                                         delim = "\t", escape_double = FALSE, 
                                         col_names = c("transcript_id",
                                                       "seq_md5",
                                                       "seq_length",
                                                       "analysis",
                                                       "signature_id",
                                                       "signature_desc",
                                                       "start_pos",
                                                       "stop_pos",
                                                       "score",
                                                       "status",
                                                       "run_date",
                                                       "interpro_id",
                                                       "interpro_desc",
                                                       "GO_annot",
                                                       "pathways"), trim_ws = TRUE) %>% 
  dplyr::select(-c(GO_annot, pathways)) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, gene_name, transcript_id, transcript_name, NMD_tx_status, NMD_tx_reason, DTE_cluster)) %>% 
  dplyr::relocate(gene_id, gene_name, transcript_id, transcript_name) %>% 
  group_by(transcript_id) %>% 
  mutate(n_unstruct_NMDRHT = n()) %>% 
  ungroup()

## Combine -----------------------------------------------------------------

NMDRHT_interPro_Pfam_MobiDBLite_combined <- NMDRHT_cds_seqs_AA_df_GENCODE_match %>% 
  dplyr::select(gene_id, gene_name, transcript_id, transcript_name, NMD_tx_status, NMD_tx_reason, DTE_cluster, width, GENCODE_width, rel_proteinWidth_AS_NMD) %>% 
  dplyr::rename("NMDRHT_width" = "width") %>% 
  left_join(NMDRHT_interPro_Pfam  %>% distinct(transcript_id, n_domains_NMDRHT)) %>% 
  left_join(GENCODE_canonical_interPro_Pfam  %>% distinct(gene_id, n_domains_GENCODE)) %>% 
  left_join(NMDRHT_interPro_MobiDBLite%>% distinct(transcript_id, n_unstruct_NMDRHT)) %>% 
  left_join(GENCODE_canonical_interPro_MobiDBLite %>% distinct(gene_id, n_unstruct_GENCODE)) %>% 
  filter(!is.na(rel_proteinWidth_AS_NMD)) %>% 
  replace(is.na(.), 0) %>% 
  as_tibble()  %>% 
  mutate(PFAM_diff = n_domains_NMDRHT-n_domains_GENCODE, .after=n_domains_GENCODE) %>% 
  mutate(unstr_diff = n_unstruct_NMDRHT-n_unstruct_GENCODE, .after=n_unstruct_GENCODE)

# Save as csv
NMDRHT_interPro_Pfam_MobiDBLite_combined %>% write_csv("Resources/NMDRHT/InterProScan/NMDRHT_interPro_Pfam_MobiDBLite_combined.csv")

#
##
###
# Proteomics --------------------------------------------------------------
###
##
#

#
###
## Potential reasons for bad detection of C-terminal peptides -------------------------------------
###
#

# Positive charge at the peptide C-terminus is advantageous for mass spectrometry analysis:
## PMID: 16159109
## https://doi.org/10.1021/ja9542193

## PMID: 34626679
### In comparison to N-termini, the analysis of protein C-termini (C-terminomics) is
### a more challenging task due to several factors.
### First of all, the ionization efficiency of C-terminal peptides derived from tryptic
### digestion is reduced compared to internal or N-terminal peptides, as, in
### the absence of missed cleavage sites, only a single protonatable amino
### group located at the peptides' amino terminus is present.

## PMID: 36959352
### Alternative splicing, which is pervasive at the transcript level, was
### previously largely undetected at the proteomics level due to the low
### degree of peptide coverage in most shotgun MS experiments 

#
###
## Whole cell proteomics (WCP) ---------------------------------------------------
###
#

##
### Peptide-Level -> Data preparation (2025-01-27) --------------------------------------------------------
##

# Scope: identify peptides that allow specific detection of NMD-annotated protein isoforms

# Overview: Use three datasources: 
# (1) in silico translated NMDRHT v1.2 - from 36898 transcripts
# (2) in silico translated ORFquant-detected ORFs (avoid overlaps with NMDRHT - i.e. focus on additional ORFs)
# (3) in silico translated Ensembl-Canonical ORFs (used to filter out potential false-positives)
# 
# Create peptidome with "cleaver" using PeptideRanger-based approach an two cleavage rule sets
# I - canonical trypsin cleavage pattern (rules: https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html)
# II - simple cleavage pattern "^M|K|R" (i.e. cut initial Methionine or after K or after R)

# From NMDRHT_ORFs.R - get ORFquant info
load("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant.rds")

# From NMDRHT_ORFs.R - get CDS info
load("Resources/NMDRHT/NMDRHT_cds_seqs_AA_df.rds")

###
#### All NMDRHT v1.2 peptides ---------------------------------------------------------
###

# Note: includes only the selected ORFs for 36898 transcripts (i.e. excluding "lncRNA"-annotated transcripts or multiple ORFs per tx)
# Problem: in "mixed"-annotated transcripts -> lack the non-NMD ORF or other potential NMD-inducing ORFs 
# Solution: (see below; next section) include all ORFquant-detected ORFs as well

NMDRHT_translation_tbl <- NMDRHT_cds_seqs_AA_df %>% 
  dplyr::select(gene_id, gene_name, transcript_id, sequence, NMD_tx_status, NMD_tx_source, NMD_tx_reason, NMD_50nt_rule, ORFquant_ORF_id_tr, ORFanage_ORF_id) %>% 
  distinct() %>% 
  dplyr::rename("Protein" =  "sequence") %>% 
  # define row number to allow precise matching after peptidome generation
  mutate(tx_row_id = row_number()) %>% 
  # detect number of stop codons (*) -> should be 1 for most transcripts (> 1 ~ selenocysteine-containing proteins?)
  mutate(stop_n = str_count(Protein, pattern = "\\*"), .after=Protein) 

##
# Determine stats
##

# Distinct transcripts
NMDRHT_translation_tbl %>% distinct(transcript_id) %>% dplyr::count()
# Distinct protein sequences
NMDRHT_translation_tbl %>% distinct(Protein) %>% dplyr::count()

# Create named vector with ORF id and sequence
NMDRHT_translation_tbl_namedVector <- NMDRHT_translation_tbl %>% 
  dplyr::select(transcript_id, Protein) %>% 
  tibble::deframe()

# Convert to AAStringSet
NMDRHT_translation_tbl_AAstring <- Biostrings::AAStringSet(NMDRHT_translation_tbl_namedVector)

# Use the "cleaver" package to perform in silico trypsin cleavage - function from "UPF1_NMDRHT_functions.R"
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
NMDRHT_translation_tbl_peptidome <- create_peptidome_VB(NMDRHT_translation_tbl_AAstring,
                                                        missed_cleavages = c(0,2), 
                                                        aa_range = c(7,52)) %>% 
  left_join(NMDRHT_translation_tbl) %>% 
  relocate(gene_id, gene_name, transcript_id)

# Use the "cleaver" package to perform in silico trypsin cleavage - SIMPLE RULES: cut_sites <- "^M|K|R" - function from "UPF1_NMDRHT_functions.R"
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
NMDRHT_translation_tbl_peptidome_simple <- create_peptidome_VB_simple(NMDRHT_translation_tbl_AAstring,
                                                                      missed_cleavages = c(0,2), 
                                                                      aa_range = c(7,52)) %>% 
  left_join(NMDRHT_translation_tbl) %>% 
  relocate(gene_id, gene_name, transcript_id)

# Join both peptidomes
NMDRHT_translation_tbl_peptidome_complete <- NMDRHT_translation_tbl_peptidome %>% 
  mutate(cleavageRules = "canonical") %>% 
  bind_rows(NMDRHT_translation_tbl_peptidome_simple %>% 
              filter(!sequence %in% NMDRHT_translation_tbl_peptidome$sequence) %>% 
              mutate(cleavageRules = "simple"))

# How many distinct peptide sequences per cleavage Rule?
NMDRHT_translation_tbl_peptidome_complete %>% 
  distinct(sequence, .keep_all = TRUE) %>% 
  dplyr::count()

###
#### All ORFquant (minus NMDRHT v1.2) peptides ---------------------------------------------------------
###

# Note: includes *all* ORFquant-detected ORFs (also those not supported by Ribo-TISH)
# Reason: ORFquant seemed to also detect rather small ORFs -> potentially better to detect more possible ORFs
# Problem: avoid overlap with already NMDRHT-derived ORFs
# Solution: filter those and retain only informative, additional ORFs

# Load NMD-sensitivity analysis on GENCODE canonical transcripts (from UPF1_NMDRHT_Annotation.R)
# Those were used to predict ORFs - which are treated by ORFquant as "annotated"
gencode.v42.annotation.canonical_NMDprediction.out <- read_csv(("Resources/GENCODE/gencode.v42.annotation.canonical_NMDprediction.csv"))

# Filter for NMD-predicted canonical transcripts
# 378 canonical transcripts are predicted to be "natural" NMD targets!
gencode.v42.annotation.canonical_NMDprediction.out_NMDplus <- gencode.v42.annotation.canonical_NMDprediction.out %>% 
  filter(is_NMD == TRUE)

# First, select all relevant ORFquant columns and retain all ORFs
ORFquant_translation_tbl <- HCT_N_AID_UPF1_ORFquant %>% 
  dplyr::select(gene_id,
                gene_name, 
                transcript_id,
                transcript_name,
                GENCODE_transcript_id,
                Protein,
                ORF_id_tr,
                ORF_category_Tx,
                ORF_type,
                NMD_status_ORFquant) %>% 
  dplyr::rename("ORFquant_ORF_id_tr" = "ORF_id_tr") %>% 
  mutate(ORFquant_NMD_reason = case_when(NMD_status_ORFquant == FALSE ~ "none",
                                         NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "exact_match" & GENCODE_transcript_id %in% gencode.v42.annotation.canonical_NMDprediction.out_NMDplus$transcript_id ~ "annotated_NMD",
                                         NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "exact_match" ~ "AS_NMD_UTR3",
                                         NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "start_match" ~ "AS_NMD",
                                         NMD_status_ORFquant == TRUE & ORF_category_Tx == "ORF_annotated" & ORF_type == "start_other" ~ "ORFanage_altStart",
                                         NMD_status_ORFquant == TRUE ~ ORF_category_Tx,
  )) %>% 
  distinct() %>% 
  # Add stop codon (*) as last character to sequence -> use later on to identify terminal peptides
  mutate(Protein = paste0(Protein,"*"))

# Check for overlap with NMDRHT table (see above) and identical NMD_tx_reason
ORFquant_translation_tbl_overlap_check <- ORFquant_translation_tbl %>% 
  # filter(ORFquant_ORF_id_tr %in% NMDRHT_translation_tbl$ORFquant_ORF_id_tr) %>% 
  dplyr::select(transcript_id, ORFquant_ORF_id_tr, NMD_status_ORFquant, ORFquant_NMD_reason) %>% 
  full_join(NMDRHT_translation_tbl %>% dplyr::select(transcript_id, ORFquant_ORF_id_tr, NMD_tx_source, NMD_tx_reason)) %>% 
  mutate(same_reason = case_when(ORFquant_NMD_reason == NMD_tx_reason ~ TRUE, 
                                 TRUE ~ FALSE))

# Simple count reveals non-identical NMD reason only for ORFanage-predicted (NMD_tx_source == ORFanage & same_reason == FALSE)
# or for ORFquant ORFs not selected for final NMDRHT transcriptome (NMD_tx_source == NA & same_reason == FALSE)
ORFquant_translation_tbl_overlap_check %>% 
  dplyr::count(NMD_tx_source, same_reason)

# Filter ORFquant ORFs for those not present in NMDRHT
# use combination of transcript_id and ORF_id for filtering
ORFquant_translation_tbl_filt <- ORFquant_translation_tbl %>% 
  filter(ORFquant_ORF_id_tr %in% (ORFquant_translation_tbl_overlap_check %>% 
                                    filter(is.na(NMD_tx_source)) %>% 
                                    pull(ORFquant_ORF_id_tr)) &
           transcript_id %in% (ORFquant_translation_tbl_overlap_check %>% 
                                 filter(is.na(NMD_tx_source)) %>% 
                                 pull(transcript_id))) %>% 
  # define row number to allow precise matching after peptidome generation
  mutate(tx_row_id = row_number()) %>% 
  # detect number of stop codons (*) -> should be 1 for most transcripts (> 1 ~ selenocysteine-containing proteins?)
  mutate(stop_n = str_count(Protein, pattern = "\\*"), .after=Protein) 

##
# Determine stats
##

# Distinct transcripts
ORFquant_translation_tbl_filt %>% distinct(transcript_id) %>% dplyr::count()
# Distinct protein sequences
ORFquant_translation_tbl_filt %>% distinct(Protein) %>% dplyr::count()

# Create named vector with ORF id and sequence
ORFquant_translation_tbl_filt_namedVector <- ORFquant_translation_tbl_filt %>% 
  dplyr::select(transcript_id, Protein) %>% 
  tibble::deframe()

# Convert to AAStringSet
ORFquant_translation_tbl_filt_AAstring <- Biostrings::AAStringSet(ORFquant_translation_tbl_filt_namedVector)

# Use the "cleaver" package to perform in silico trypsin cleavage
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
ORFquant_translation_tbl_filt_peptidome <- create_peptidome_VB(ORFquant_translation_tbl_filt_AAstring,
                                                               missed_cleavages = c(0,2), 
                                                               aa_range = c(7,52)) %>% 
  left_join(ORFquant_translation_tbl_filt) %>% 
  relocate(gene_id, gene_name, transcript_id)

# Use the "cleaver" package to perform in silico trypsin cleavage - SIMPLE RULES: cut_sites <- "^M|K|R"
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
ORFquant_translation_tbl_filt_peptidome_simple <- create_peptidome_VB_simple(ORFquant_translation_tbl_filt_AAstring,
                                                                             missed_cleavages = c(0,2), 
                                                                             aa_range = c(7,52)) %>% 
  left_join(ORFquant_translation_tbl_filt) %>% 
  relocate(gene_id, gene_name, transcript_id)

# Join both peptidomes
ORFquant_translation_tbl_filt_peptidome_complete <- ORFquant_translation_tbl_filt_peptidome %>% 
  mutate(cleavageRules = "canonical") %>% 
  bind_rows(ORFquant_translation_tbl_filt_peptidome_simple %>% 
              filter(!sequence %in% ORFquant_translation_tbl_filt_peptidome$sequence) %>% 
              mutate(cleavageRules = "simple"))

# How many distinct peptide sequences per cleavage Rule?
ORFquant_translation_tbl_filt_peptidome_complete %>% 
  distinct(sequence, .keep_all = TRUE) %>% 
  dplyr::count()

###
#### GENCODE canonical sequences ---------------------------------------------
###

# Note: includes *all* Ensembl-canonical in silico translated proteins 
# Reason: serves as source to filter out false-positives. Peptides found here are very very unlikely to be NMD-specific
# Problem: many (206) translated proteins with more than 1 stop codon -> potentially selenocysteine-containing proteins?
# Solution: simply filtering those out does not help and rather leads to false-positives. At the moment, retain those proteins

load("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds_AA_df_forExport.rds")

# Add row ID for joining later on
gencode.v42.annotation_txDb_cds_AA <- gencode.v42.annotation_canonical_txDb_cds_AA_df_forExport %>% 
  # Add stop codon (*) as last character to sequence -> use later on to identify terminal peptides
  mutate(Protein = paste0(Protein,"*")) %>% 
  # define row number to allow precise matching after peptidome generation
  mutate(tx_row_id = row_number()) %>% 
  # detect number of stop codons (*) -> should be 1 for most transcripts (> 1 ~ selenocysteine-containing proteins?)
  mutate(stop_n = str_count(Protein, pattern = "\\*"), .after=Protein) 

##
# Determine stats
##

# Distinct transcripts
gencode.v42.annotation_txDb_cds_AA %>% distinct(transcript_id) %>% dplyr::count()
# Distinct protein sequences
gencode.v42.annotation_txDb_cds_AA %>% distinct(Protein) %>% dplyr::count()

# Create named vector with ORF id and sequence
gencode_canonical_tbl_namedVector <- gencode.v42.annotation_txDb_cds_AA %>% 
  # filter(gene_id %in% AS_NMD_ORFquant_tbl$gene_id) %>%
  dplyr::select(transcript_id, Protein) %>% 
  tibble::deframe()

# Convert to AAStringSet
gencode_canonical_tbl_AAstring <- Biostrings::AAStringSet(gencode_canonical_tbl_namedVector)

# Use the "cleaver" package to perform in silico trypsin cleavage
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
gencode_canonical_tbl_peptidome <- create_peptidome_VB(gencode_canonical_tbl_AAstring,
                                                       missed_cleavages = c(0,2), 
                                                       aa_range = c(7,52)) %>% 
  left_join(gencode.v42.annotation_txDb_cds_AA) %>% 
  relocate(gene_id, gene_name, transcript_id) %>% 
  mutate(NMD_tx_reason = "none")

# Use the "cleaver" package to perform in silico trypsin cleavage - SIMPLE RULES: cut_sites <- "^M|K|R"
# Filter for at least 7 aa long peptides (shortest in proteomics peptide list)
gencode_canonical_tbl_peptidome_simple <- create_peptidome_VB_simple(gencode_canonical_tbl_AAstring,
                                                                     missed_cleavages = c(0,2), 
                                                                     aa_range = c(7,52)) %>% 
  left_join(gencode.v42.annotation_txDb_cds_AA) %>% 
  relocate(gene_id, gene_name, transcript_id) %>% 
  mutate(NMD_tx_reason = "none")

# Join both peptidomes
gencode_canonical_tbl_peptidome_complete <- gencode_canonical_tbl_peptidome %>% 
  mutate(cleavageRules = "canonical") %>% 
  bind_rows(gencode_canonical_tbl_peptidome_simple %>% 
              filter(!sequence %in% gencode_canonical_tbl_peptidome$sequence) %>% 
              mutate(cleavageRules = "simple"))

# How many distinct peptide sequences per cleavage Rule?
gencode_canonical_tbl_peptidome_complete %>% 
  distinct(sequence, .keep_all = TRUE) %>% 
  dplyr::count()

###
#### Export ORFs for Spectronaut Analysis  -----------------------------------
###

## First for NMDRHT v1.2 -> select ORFquant-annotated ORFs per transcript
NMDRHT_translation_tbl_distinct_ORFquant_forMS <- NMDRHT_translation_tbl %>% 
  filter(!is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFquant_ORF_id_tr, transcript_id, Protein)

# Prepare for FASTA export
NMDRHT_translation_tbl_distinct_ORFquant_forMS_fasta <- c(rbind(paste0(">", 
                                                                       NMDRHT_translation_tbl_distinct_ORFquant_forMS$gene_name, 
                                                                       "_",
                                                                       NMDRHT_translation_tbl_distinct_ORFquant_forMS$ORFquant_ORF_id_tr, 
                                                                       "_",
                                                                       NMDRHT_translation_tbl_distinct_ORFquant_forMS$transcript_id), 
                                                                NMDRHT_translation_tbl_distinct_ORFquant_forMS$Protein))

head(NMDRHT_translation_tbl_distinct_ORFquant_forMS_fasta)

## Second for NMDRHT v1.2 -> select ORFanage-predicted ORFs per transcript
NMDRHT_translation_tbl_distinct_ORFanage_forMS <- NMDRHT_translation_tbl %>% 
  filter(is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFanage_ORF_id, transcript_id, Protein) %>% 
  separate(ORFanage_ORF_id, c("start", "end")) %>% 
  mutate(ORFanage_ORF_id_exp = paste0("ORFanage", "_", start, "_", end))

# Prepare for FASTA export
NMDRHT_translation_tbl_distinct_ORFanage_forMS_fasta <- c(rbind(paste0(">", 
                                                                       NMDRHT_translation_tbl_distinct_ORFanage_forMS$gene_name, 
                                                                       "_",
                                                                       NMDRHT_translation_tbl_distinct_ORFanage_forMS$ORFanage_ORF_id_exp, 
                                                                       "_",
                                                                       NMDRHT_translation_tbl_distinct_ORFanage_forMS$transcript_id), 
                                                                NMDRHT_translation_tbl_distinct_ORFanage_forMS$Protein))

head(NMDRHT_translation_tbl_distinct_ORFanage_forMS_fasta)

## As last source -> select (extra) ORFquant ORFs
ORFquant_translation_tbl_filt_forMS <- ORFquant_translation_tbl_filt %>% 
  filter(!is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFquant_ORF_id_tr, transcript_id, Protein) 

# Prepare for FASTA export
ORFquant_translation_tbl_filt_forMS_fasta <- c(rbind(paste0(">", 
                                                            ORFquant_translation_tbl_filt_forMS$gene_name, 
                                                            "_",
                                                            ORFquant_translation_tbl_filt_forMS$ORFquant_ORF_id_tr, 
                                                            "_",
                                                            ORFquant_translation_tbl_filt_forMS$transcript_id), 
                                                     ORFquant_translation_tbl_filt_forMS$Protein))

head(ORFquant_translation_tbl_filt_forMS_fasta)

# Combine all three FASTA data
combined_forMS_fasta <- c(NMDRHT_translation_tbl_distinct_ORFquant_forMS_fasta,
                          NMDRHT_translation_tbl_distinct_ORFanage_forMS_fasta,
                          ORFquant_translation_tbl_filt_forMS_fasta)

write(x = combined_forMS_fasta, file = "Resources/Proteomics/HCT_N_AID_UPF1_combined_forMS_fasta_20250123.fa")

# combine as data frame
combined_forMS <- NMDRHT_translation_tbl_distinct_ORFquant_forMS %>% 
  mutate(source = "NMDRHT") %>% 
  bind_rows(NMDRHT_translation_tbl_distinct_ORFanage_forMS %>% 
              mutate(source = "ORFanage")) %>% 
  bind_rows(ORFquant_translation_tbl_filt_forMS %>% 
              mutate(source = "ORFquant_extra")) %>% 
  mutate(stop_n = str_count(Protein, pattern = "\\*"), .after=Protein) %>% 
  mutate(id = case_when(source == "ORFanage" ~ paste0(gene_name, 
                                                      "_",
                                                      ORFanage_ORF_id_exp, 
                                                      "_",
                                                      transcript_id),
                        source != "ORFanage" ~ paste0(gene_name, 
                                                      "_",
                                                      ORFquant_ORF_id_tr, 
                                                      "_",
                                                      transcript_id)))

combined_forMS %>% 
  dplyr::count(stop_n)

combined_forMS %>% 
  distinct(Protein) %>% 
  dplyr::count()

###
#### w/o Stop - Export ORFs for Spectronaut Analysis  -----------------------------------
###

## First for NMDRHT v1.2 -> select ORFquant-annotated ORFs per transcript
NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop <- NMDRHT_translation_tbl %>% 
  filter(!is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFquant_ORF_id_tr, transcript_id, Protein) %>% 
  mutate(Protein = stringr::str_replace(Protein, '\\*', ''))

# Prepare for FASTA export
NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop_fasta <- c(rbind(paste0(">", 
                                                                              NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop$gene_name, 
                                                                              "_",
                                                                              NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop$ORFquant_ORF_id_tr, 
                                                                              "_",
                                                                              NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop$transcript_id), 
                                                                       NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop$Protein))

head(NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop_fasta)

## Second for NMDRHT v1.2 -> select ORFanage-predicted ORFs per transcript
NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop <- NMDRHT_translation_tbl %>% 
  filter(is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFanage_ORF_id, transcript_id, Protein) %>% 
  separate(ORFanage_ORF_id, c("start", "end")) %>% 
  mutate(ORFanage_ORF_id_exp = paste0("ORFanage", "_", start, "_", end)) %>% 
  mutate(Protein = stringr::str_replace(Protein, '\\*', ''))

# Prepare for FASTA export
NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop_fasta <- c(rbind(paste0(">", 
                                                                              NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop$gene_name, 
                                                                              "_",
                                                                              NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop$ORFanage_ORF_id_exp, 
                                                                              "_",
                                                                              NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop$transcript_id), 
                                                                       NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop$Protein))

head(NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop_fasta)

## As last source -> select (extra) ORFquant ORFs
ORFquant_translation_tbl_filt_forMS_woStop <- ORFquant_translation_tbl_filt %>% 
  filter(!is.na(ORFquant_ORF_id_tr)) %>% 
  distinct(gene_name, ORFquant_ORF_id_tr, transcript_id, Protein) %>% 
  mutate(Protein = stringr::str_replace(Protein, '\\*', ''))

# Prepare for FASTA export
ORFquant_translation_tbl_filt_forMS_woStop_fasta <- c(rbind(paste0(">", 
                                                                   ORFquant_translation_tbl_filt_forMS_woStop$gene_name, 
                                                                   "_",
                                                                   ORFquant_translation_tbl_filt_forMS_woStop$ORFquant_ORF_id_tr, 
                                                                   "_",
                                                                   ORFquant_translation_tbl_filt_forMS_woStop$transcript_id), 
                                                            ORFquant_translation_tbl_filt_forMS_woStop$Protein))

head(ORFquant_translation_tbl_filt_forMS_woStop_fasta)

# Combine all three FASTA data
combined_forMS_fasta_woStop <- c(NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop_fasta,
                                 NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop_fasta,
                                 ORFquant_translation_tbl_filt_forMS_woStop_fasta)

write(x = combined_forMS_fasta_woStop, file = "Resources/Proteomics/HCT_N_AID_UPF1_combined_forMS_fasta_woStop_20250123.fa")

# combine as data frame
combined_forMS_woStop <- NMDRHT_translation_tbl_distinct_ORFquant_forMS_woStop %>% 
  bind_rows(NMDRHT_translation_tbl_distinct_ORFanage_forMS_woStop) %>% 
  bind_rows(ORFquant_translation_tbl_filt_forMS_woStop) %>% 
  mutate(stop_n = str_count(Protein, pattern = "\\*"), .after=Protein) 

combined_forMS_woStop %>% 
  dplyr::count(stop_n)

combined_forMS_woStop %>% 
  distinct(Protein) %>% 
  dplyr::count()

###
#### Join and classify  -----------------------------------
###

# Combine NMDRHT-, ORFquant- and Ensembl-Canonical-derived peptidomes
# harmonize certain columns
combined_tbl_peptidome_complete <- NMDRHT_translation_tbl_peptidome_complete %>% 
  mutate(peptide_set = "NMDRHT") %>% 
  bind_rows(ORFquant_translation_tbl_filt_peptidome_complete %>% 
              dplyr::select(-c(transcript_name, GENCODE_transcript_id, ORF_category_Tx, ORF_type)) %>% 
              mutate(NMD_tx_status = case_when(NMD_status_ORFquant == TRUE ~ "NMD",
                                               NMD_status_ORFquant == FALSE ~ "coding",
                                               TRUE ~ "lncRNA")) %>% 
              mutate(NMD_tx_source = "ORFquant_only") %>% 
              dplyr::rename("NMD_tx_reason"="ORFquant_NMD_reason",
                            "NMD_50nt_rule" = "NMD_status_ORFquant") %>% 
              mutate(peptide_set = "ORFquant")) %>% 
  bind_rows(gencode_canonical_tbl_peptidome_complete %>% 
              dplyr::select(-c(transcript_name)) %>% 
              mutate(NMD_tx_status = "coding") %>% 
              mutate(NMD_tx_source = "canonical") %>% 
              mutate(NMD_50nt_rule = FALSE) %>% 
              mutate(peptide_set = "canonical"))

# Classify each peptide and determine whether NMD reason can be defined uniquely
# Rationale: each peptide that is exclusively found in NMD-annotated proteins and *NOT* in canonical or coding -> designate as NMD
# Problem: NMD_reason can be ambigiuous depending on the transcript-context
# Solution: assign "multiple" if NMD_reason cannot be uniquely assigned

combined_tbl_peptidome_complete_class <- combined_tbl_peptidome_complete %>% 
  group_by(sequence) %>% 
  mutate(n_NMD = sum(NMD_50nt_rule == TRUE),
         n_noNMD = sum(NMD_50nt_rule == FALSE)) %>% 
  mutate(n_unique_NMD_tx_reason = length(unique(NMD_tx_reason))) %>% 
  ungroup() %>% 
  mutate(peptide_class = case_when(n_NMD > 0 & n_noNMD > 0 ~ "shared",
                                   n_NMD > 0 & n_noNMD == 0  ~ "NMD",
                                   n_NMD == 0 & n_noNMD > 0 ~ "coding")) %>% 
  mutate(NMD_peptide_status = case_when(n_unique_NMD_tx_reason == 1 ~ "unique",
                                        TRUE ~ "multiple")) %>% 
  mutate(NMD_peptide_reason = case_when(NMD_peptide_status == "unique" ~ NMD_tx_reason,
                                        NMD_peptide_status != "unique" ~ "multiple"))

# Summary stats
combined_tbl_peptidome_complete_class %>% 
  distinct(sequence) %>% 
  dplyr::count()

###
#### Junction-peptides -------------------------------
###

##### NMDRHT ------------------------------------------------------------------

# load NMDRHT cds information (GRangesList)
load("Resources/NMDRHT/NMDRHT_cds.rds")

# Determine exon-exon junction positions relative to ORF
NMDRHT_cds_df <- as_tibble(NMDRHT_cds) %>% 
  group_by(group_name) %>% 
  # define rank based on strand
  mutate(cds_rank = case_when(strand == "+" ~ rank(start),
                              strand == "-" ~ rank(-start))) %>% 
  # convert from genomic to protein level coordinates
  mutate(width_aa = width/3) %>%
  mutate(cumulative_width_aa=cumsum(width_aa)) %>% 
  ungroup() 

# Match peptidome with cds information
NMDRHT_peptidome_complete_class_junction <- combined_tbl_peptidome_complete_class %>% 
  filter(peptide_set == "NMDRHT") %>% 
  # match only by transcript_id
  left_join(NMDRHT_cds_df %>% 
              dplyr::select(-c(start,end)) %>% 
              dplyr::rename("transcript_id" = "group_name")) %>% 
  # check if reading frame is in phase 0, 1, 2 (round number, .3 or .6)
  mutate(cum_width_aa_round = round(cumulative_width_aa,0)) %>% 
  # calculate distance between peptide start or end with exon-exon junction (EEJ)
  mutate(dist_start_EEJ = start - cumulative_width_aa,
         dist_end_EEJ = end - cumulative_width_aa) %>% 
  # define rules to detect junction-spanning peptides
  # different rules depending on phase of translation - note corner case for phase 1 or 2
  mutate(peptide_overJunction = case_when(dist_start_EEJ < 0 & dist_end_EEJ <= 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ >= 0 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ >= 1 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ < 1 & dist_end_EEJ > 0 ~ TRUE,
                                          dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          TRUE ~ FALSE)) %>% 
  # determine absolute distance of peptide from EEJ
  mutate(abs_diffSum = abs(dist_start_EEJ) + abs(dist_end_EEJ)) %>% 
  group_by(transcript_id, sequence, tx_row_id, start, end, ORFquant_ORF_id_tr, peptide_set) %>% 
  # select only nearest EEJ to keep for each peptide
  slice_min(abs_diffSum, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  # Determine exon number the peptide comes from (X_Y for spanning peptides)
  mutate(peptide_inExon = case_when(peptide_overJunction == TRUE ~ paste0(exon_rank,"_",exon_rank+1),
                                    peptide_overJunction == FALSE & dist_start_EEJ < 0 & dist_end_EEJ < 0 ~ as.character(exon_rank),
                                    peptide_overJunction == FALSE & dist_start_EEJ > 0 & dist_end_EEJ > 0 ~ as.character(exon_rank+1))) %>% 
  dplyr::rename("dist_ref_exon" = "exon_rank") %>% 
  relocate(dist_ref_exon,
           dist_start_EEJ, 
           dist_end_EEJ,
           peptide_inExon, 
           peptide_overJunction, .after=size) %>% 
  dplyr::select(-c(cds_rank,
                   width_aa, 
                   cumulative_width_aa,
                   cum_width_aa_round, 
                   abs_diffSum, 
                   group,
                   seqnames,
                   width,
                   strand,
                   cds_id,
                   cds_name))

##### ORFquant ------------------------------------------------------------------

# load ORFquant cds information (GRangesList)
load("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant_cds.rds")

# Determine exon-exon junction positions relative to ORF
HCT_N_AID_UPF1_ORFquant_cds_df <- as_tibble(HCT_N_AID_UPF1_ORFquant_cds) %>% 
  group_by(group_name) %>% 
  # define rank based on strand
  mutate(cds_rank = case_when(strand == "+" ~ rank(start),
                              strand == "-" ~ rank(-start))) %>% 
  # convert from genomic to protein level coordinates
  mutate(width_aa = width/3) %>%
  mutate(cumulative_width_aa=cumsum(width_aa)) %>% 
  ungroup() 

# Match peptidome with cds information
ORFquant_peptidome_complete_class_junction <- combined_tbl_peptidome_complete_class %>% 
  filter(peptide_set == "ORFquant") %>% 
  # match only by transcript_id
  left_join(HCT_N_AID_UPF1_ORFquant_cds_df %>% 
              separate(group_name, c("transcript_id", "ORFquant_ORF_id_tr"),
                       sep = "_",
                       extra = "merge") %>% 
              dplyr::select(-c(start,end)) %>% 
              dplyr::rename("transcript_id" = "transcript_id")) %>% 
  # check if reading frame is in phase 0, 1, 2 (round number, .3 or .6)
  mutate(cum_width_aa_round = round(cumulative_width_aa,0)) %>% 
  # calculate distance between peptide start or end with exon-exon junction (EEJ)
  mutate(dist_start_EEJ = start - cumulative_width_aa,
         dist_end_EEJ = end - cumulative_width_aa) %>% 
  # define rules to detect junction-spanning peptides
  # different rules depending on phase of translation - note corner case for phase 1 or 2
  mutate(peptide_overJunction = case_when(dist_start_EEJ < 0 & dist_end_EEJ <= 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ >= 0 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ >= 1 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ < 1 & dist_end_EEJ > 0 ~ TRUE,
                                          dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          TRUE ~ FALSE)) %>% 
  # determine absolute distance of peptide from EEJ
  mutate(abs_diffSum = abs(dist_start_EEJ) + abs(dist_end_EEJ)) %>% 
  group_by(transcript_id, sequence, tx_row_id, start, end, ORFquant_ORF_id_tr, peptide_set) %>% 
  # select only nearest EEJ to keep for each peptide
  slice_min(abs_diffSum, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  # Determine exon number the peptide comes from (X_Y for spanning peptides)
  mutate(peptide_inExon = case_when(peptide_overJunction == TRUE ~ paste0(exon_rank,"_",exon_rank+1),
                                    peptide_overJunction == FALSE & dist_start_EEJ < 0 & dist_end_EEJ < 0 ~ as.character(exon_rank),
                                    peptide_overJunction == FALSE & dist_start_EEJ > 0 & dist_end_EEJ > 0 ~ as.character(exon_rank+1))) %>% 
  dplyr::rename("dist_ref_exon" = "exon_rank") %>% 
  relocate(dist_ref_exon,
           dist_start_EEJ, 
           dist_end_EEJ,
           peptide_inExon, 
           peptide_overJunction, .after=size) %>% 
  dplyr::select(-c(cds_rank,
                   width_aa, 
                   cumulative_width_aa,
                   cum_width_aa_round, 
                   abs_diffSum, 
                   group,
                   seqnames,
                   width,
                   strand,
                   cds_id,
                   cds_name))

##### GENCODE canonical ------------------------------------------------------------------

# load gencode.v42 ensembl-canonical cds information (GRangesList)
load("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds.rds")

# Determine exon-exon junction positions relative to ORF
gencode.v42.annotation_txDb_cds_df <- as_tibble(gencode.v42.annotation_canonical_txDb_cds) %>% 
  group_by(group_name) %>% 
  # define rank based on strand
  mutate(cds_rank = case_when(strand == "+" ~ rank(start),
                              strand == "-" ~ rank(-start))) %>% 
  # convert from genomic to protein level coordinates
  mutate(width_aa = width/3) %>%
  mutate(cumulative_width_aa=cumsum(width_aa)) %>% 
  ungroup() 

# Match peptidome with cds information
gencode_canonical_peptidome_complete_class_junction <- combined_tbl_peptidome_complete_class %>% 
  filter(peptide_set == "canonical") %>% 
  # match only by transcript_id
  left_join(gencode.v42.annotation_txDb_cds_df %>% 
              dplyr::select(-c(start,end)) %>% 
              dplyr::rename("transcript_id" = "group_name")) %>% 
  # check if reading frame is in phase 0, 1, 2 (round number, .3 or .6)
  mutate(cum_width_aa_round = round(cumulative_width_aa,0)) %>% 
  # calculate distance between peptide start or end with exon-exon junction (EEJ)
  mutate(dist_start_EEJ = start - cumulative_width_aa,
         dist_end_EEJ = end - cumulative_width_aa) %>% 
  # define rules to detect junction-spanning peptides
  # different rules depending on phase of translation - note corner case for phase 1 or 2
  mutate(peptide_overJunction = case_when(dist_start_EEJ < 0 & dist_end_EEJ <= 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ >= 0 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round == cumulative_width_aa & dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ >= 1 & dist_end_EEJ > 0 ~ FALSE,
                                          cum_width_aa_round != cumulative_width_aa & dist_start_EEJ < 1 & dist_end_EEJ > 0 ~ TRUE,
                                          dist_start_EEJ < 0 & dist_end_EEJ > 0 ~ TRUE,
                                          TRUE ~ FALSE)) %>% 
  # determine absolute distance of peptide from EEJ
  mutate(abs_diffSum = abs(dist_start_EEJ) + abs(dist_end_EEJ)) %>% 
  group_by(transcript_id, sequence, tx_row_id, start, end, peptide_set) %>% 
  # select only nearest EEJ to keep for each peptide
  slice_min(abs_diffSum, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  # Determine exon number the peptide comes from (X_Y for spanning peptides)
  mutate(peptide_inExon = case_when(peptide_overJunction == TRUE ~ paste0(exon_rank,"_",exon_rank+1),
                                    peptide_overJunction == FALSE & dist_start_EEJ < 0 & dist_end_EEJ < 0 ~ as.character(exon_rank),
                                    peptide_overJunction == FALSE & dist_start_EEJ > 0 & dist_end_EEJ > 0 ~ as.character(exon_rank+1))) %>% 
  dplyr::rename("dist_ref_exon" = "exon_rank") %>% 
  relocate(dist_ref_exon,
           dist_start_EEJ, 
           dist_end_EEJ,
           peptide_inExon, 
           peptide_overJunction, .after=size) %>% 
  dplyr::select(-c(cds_rank,
                   width_aa, 
                   cumulative_width_aa,
                   cum_width_aa_round, 
                   abs_diffSum, 
                   group,
                   seqnames,
                   width,
                   strand,
                   cds_id,
                   cds_name))


##### Combine -----------------------------------------------------------------
combined_peptidome_complete_class_junction <- NMDRHT_peptidome_complete_class_junction %>% 
  bind_rows(ORFquant_peptidome_complete_class_junction) %>% 
  bind_rows(gencode_canonical_peptidome_complete_class_junction)

# Stats
combined_peptidome_complete_class_junction %>% 
  dplyr::count(peptide_set, peptide_overJunction)

# Save data
save(combined_peptidome_complete_class_junction,
     file = paste0("Resources/Proteomics/combined_peptidome_complete_class_junction.rds"))

# load("Resources/Proteomics/combined_peptidome_complete_class_junction.rds")

###
#### Distinct peptides -----------------------------------
###

# Determine duplicate sequences - obtain Unique peptides based on few important parameters
combined_tbl_peptidome_class_unique <- combined_peptidome_complete_class_junction %>% 
  group_by(sequence) %>% 
  mutate(duplicated = case_when(n() > 1 ~ TRUE,
                                n() == 1 ~ FALSE)) %>% 
  ungroup() %>% 
  distinct(sequence, .keep_all = TRUE) 

# Stats
combined_tbl_peptidome_class_unique %>% 
  dplyr::count()

###
#### Peptide detectability  -------------------------------
###

# Use PeptideRanger - peptide_predictions

# Perform on distinct sequences - no need to process duplications

# Perform peptide detectability analysis using the Random forest classifier trained on the ProteomicsDB database
combined_tbl_peptidome_class_unique_predictions <- PeptideRanger::peptide_predictions(peptides = combined_tbl_peptidome_class_unique$sequence,
                                                                       prediction_model = RFmodel_ProteomicsDB,
                                                                       missed_cleavages = combined_tbl_peptidome_class_unique$missed_cleavages)

# Combine information
combined_tbl_peptidome_class_unique_detectability <- combined_tbl_peptidome_class_unique %>% 
  mutate(RF_score = combined_tbl_peptidome_class_unique_predictions$RF_score)

# Save data
save(combined_tbl_peptidome_class_unique_detectability,
     file = paste0("Resources/Proteomics/combined_tbl_peptidome_class_unique_detectability.rds"))

# load("Resources/Proteomics/combined_tbl_peptidome_class_unique_detectability.rds")

# Export 
combined_tbl_peptidome_class_unique_detectability %>% write_csv("Resources/Proteomics/peptide_classes.csv")



###
### SN Detected peptides (27.01.25) -----------------------------------
###

# Code adapted from Oliver Popp (MDC)


#### Data import & tidy-up -------------------------------------------------------------

# import of the converted peptide-level data
WCP_peptide_data_20250127_SN <- read_delim("Resources/Proteomics/2025_01_27_Reanalysis/UPF1_NMDRHT_WCP_peptides_SN.txt")

# same Q values as in protein-based analysis used
# Check length of amino acids - should be between 7-52
c(min(nchar(WCP_peptide_data_20250127_SN$Sequence)), max(nchar(WCP_peptide_data_20250127_SN$Sequence)))

# Check dimensions 
dim(WCP_peptide_data_20250127_SN)

# remove unnecessary strings and rename experiments for better sorting
names(WCP_peptide_data_20250127_SN) <- gsub("Elmo_|Fozzie_|[0-9]{8,8}_|Landthaler_|HS_|HSdia_|PaWu_|UPF1_|OP_|A[0-9]+_|B[0-9]+_",
                                            "", names(WCP_peptide_data_20250127_SN))
names(WCP_peptide_data_20250127_SN) <- sub("6h", "06h", names(WCP_peptide_data_20250127_SN))
names(WCP_peptide_data_20250127_SN) <- sub("hDMSO", "h_DMSO", names(WCP_peptide_data_20250127_SN))
names(WCP_peptide_data_20250127_SN) <- sub("h_Depl", "h_IAA", names(WCP_peptide_data_20250127_SN))
names(WCP_peptide_data_20250127_SN) <- sub("hDepl", "h_IAA", names(WCP_peptide_data_20250127_SN))
names(WCP_peptide_data_20250127_SN)

# Generate data frame
WCP_peptide_data_20250127_SN_df <- WCP_peptide_data_20250127_SN %>% dplyr::select(starts_with("Intensity."))
WCP_peptide_data_20250127_SN_df <- WCP_peptide_data_20250127_SN_df[,order(names(WCP_peptide_data_20250127_SN_df))]
names(WCP_peptide_data_20250127_SN_df)

# combine relevant columns from original dataframe with adapted columns
WCP_peptide_data_20250127_SN_final <- cbind(WCP_peptide_data_20250127_SN[,1:3], WCP_peptide_data_20250127_SN_df)
names(WCP_peptide_data_20250127_SN_final)

# the output can be used for peptide maps etc.
write_tsv(WCP_peptide_data_20250127_SN_final, "Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250125_SN_final_peptides.txt")

# Any duplicated sequences? Should be false
any(duplicated(WCP_peptide_data_20250127_SN_final$Sequence))

# sanity checks using a matrix version of the peptides data -- log2 transformation
WCP_peptide_data_20250127_SN_matrix <- WCP_peptide_data_20250127_SN_final %>% 
  dplyr::select(Sequence, starts_with("Intensity")) %>% column_to_rownames("Sequence")

WCP_peptide_data_20250127_SN_matrix <- as.matrix(WCP_peptide_data_20250127_SN_matrix)
head(WCP_peptide_data_20250127_SN_matrix)

# log-transform
WCP_peptide_data_20250127_SN_matrix <- mylog(WCP_peptide_data_20250127_SN_matrix)
# myboxplot(WCP_peptide_data_20250127_SN_matrix) # own function

# remove outlier (12h_Depl_A) from matrix and dataframe
WCP_peptide_data_20250127_SN_matrix <- WCP_peptide_data_20250127_SN_matrix[,-grep("12h_IAA_A", colnames(WCP_peptide_data_20250127_SN_matrix))]
# myboxplot(WCP_peptide_data_20250127_SN_matrix)
WCP_peptide_data_20250127_SN_final <- WCP_peptide_data_20250127_SN_final[,-grep("12h_IAA_A", names(WCP_peptide_data_20250127_SN_final))]

### sanity checks -- renaming of columns -- merging of linear-space and log2 data
colnames(WCP_peptide_data_20250127_SN_matrix) <- sub("Intensity.", "log2_Intensity.", colnames(WCP_peptide_data_20250127_SN_matrix))
WCP_peptide_data_20250127_SN_matrix_df <- as.data.frame(WCP_peptide_data_20250127_SN_matrix) %>% rownames_to_column("Sequence")
head(WCP_peptide_data_20250127_SN_matrix_df)
names(WCP_peptide_data_20250127_SN_matrix_df);names(WCP_peptide_data_20250127_SN_final)
WCP_peptide_data_20250127_SN_final_log2 <- merge(WCP_peptide_data_20250127_SN_final, WCP_peptide_data_20250127_SN_matrix_df, by = "Sequence")
dim(WCP_peptide_data_20250127_SN_final);dim(WCP_peptide_data_20250127_SN_matrix_df);dim(WCP_peptide_data_20250127_SN_final_log2)
head(WCP_peptide_data_20250127_SN_final_log2)

#### Merge with classes -------------------------------------------------------------

# merge peptide data with peptide classes (from above)
WCP_peptide_data_20250127_SN_final_log2_combined <- merge(WCP_peptide_data_20250127_SN_final_log2, 
                                                          combined_tbl_peptidome_class_unique_detectability, 
                                                          by.x = "Sequence", by.y = "sequence", all.x = TRUE, all.y = FALSE)
head(WCP_peptide_data_20250127_SN_final_log2_combined)

# sanity check
table(WCP_peptide_data_20250127_SN_final_log2_combined$peptide_class)
length(which(is.na(WCP_peptide_data_20250127_SN_final_log2_combined$peptide_class))) 
# 1144/162949 (~ 1%) sequences are not assigned to a peptide class

# check which peptides have no classification 
# sanity check: should not contain NMDRHT-derived ORFs
WCP_peptide_data_20250127_SN_final_log2_combined %>% 
  filter(is.na(peptide_class)) %>% 
  mutate(source = case_when(str_detect(Protein.Group, "zzCont") ~ "zzCont",
                            str_detect(Protein.Group, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  dplyr::count(source)

dim(WCP_peptide_data_20250127_SN_final_log2_combined)
any(duplicated(WCP_peptide_data_20250127_SN_final_log2_combined$Sequence))
head(WCP_peptide_data_20250127_SN_final_log2_combined)

# stats
WCP_peptide_data_20250127_SN_final_log2_combined %>% 
  dplyr::count(peptide_class,peptide_overJunction)

# save as csv
WCP_peptide_data_20250127_SN_final_log2_combined %>% write_csv("Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250127_SN_final_log2_combined.csv")

# Filter for NMD-informative pepetides
WCP_peptide_data_20250127_SN_final_log2_combined_NMD <- WCP_peptide_data_20250127_SN_final_log2_combined %>% 
  filter(peptide_class == "NMD")

# Check NMD peptide reason
WCP_peptide_data_20250127_SN_final_log2_combined_NMD %>% 
  dplyr::count(NMD_peptide_reason)

#### sample-level ------------------------------------------------------------

# Pivot longer
WCP_peptide_data_20250127_SN_final_log2_combined_expand <- WCP_peptide_data_20250127_SN_final_log2_combined %>% 
  # clean_names() %>% 
  mutate(source = case_when(str_detect(Protein.Group, "zzCont") ~ "zzCont",
                            str_detect(Protein.Group, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  pivot_longer(cols = c(starts_with("Intensity."), starts_with("log2_Intensity.")),
               names_to = c(".value", "sample"),
               names_sep = "\\.") %>% 
  mutate(condition = case_when(str_detect(sample, "12h_DMSO") ~ "12h_DMSO",
                               str_detect(sample, "12h_IAA") ~ "12h_IAA",
                               str_detect(sample, "24h_DMSO") ~ "24h_DMSO",
                               str_detect(sample, "24h_IAA") ~ "24h_IAA",
                               str_detect(sample, "06h_DMSO") ~ "06h_DMSO",
                               str_detect(sample, "06h_IAA") ~ "06h_IAA",
                               str_detect(sample, "15h_DMSO") ~ "15h_DMSO",
                               str_detect(sample, "15h_IAA") ~ "15h_IAA",
                               str_detect(sample, "18h_DMSO") ~ "18h_DMSO",
                               str_detect(sample, "18h_IAA") ~ "18h_IAA"))    

# save as csv
WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% write_csv("Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250127_SN_final_log2_combined_expand.csv")

WCP_peptide_data_20250127_SN_final_log2_combined_expand_NMD <- WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  get_summary_stats(log2_Intensity)

###
### DiaNN Detected peptides (06.02.25) -----------------------------------
###

# Code adapted from Oliver Popp (MDC)

#### Data import & tidy-up -------------------------------------------------------------

# import of the converted peptide-level data
WCP_peptide_data_20250206_DiaNN <- read_delim("Resources/Proteomics/2025_02_06_Reanalysis/UPF1_NMDRHT_WCP_peptides_DiaNN.txt")

# same Q values as in protein-based analysis used
# Check length of amino acids
c(min(nchar(WCP_peptide_data_20250206_DiaNN$Sequence)), max(nchar(WCP_peptide_data_20250206_DiaNN$Sequence)))

# Check dimensions 
dim(WCP_peptide_data_20250206_DiaNN)

# remove unnecessary strings and rename experiments for better sorting
names(WCP_peptide_data_20250206_DiaNN) <- gsub("Elmo_|Fozzie_|[0-9]{8,8}_|Landthaler_|HS_|HSdia_|PaWu_|UPF1_|OP_|A[0-9]+_|B[0-9]+_",
                                               "", names(WCP_peptide_data_20250206_DiaNN))
names(WCP_peptide_data_20250206_DiaNN) <- sub("6h", "06h", names(WCP_peptide_data_20250206_DiaNN))
names(WCP_peptide_data_20250206_DiaNN) <- sub("hDMSO", "h_DMSO", names(WCP_peptide_data_20250206_DiaNN))
names(WCP_peptide_data_20250206_DiaNN) <- sub("h_Depl", "h_IAA", names(WCP_peptide_data_20250206_DiaNN))
names(WCP_peptide_data_20250206_DiaNN) <- sub("hDepl", "h_IAA", names(WCP_peptide_data_20250206_DiaNN))
names(WCP_peptide_data_20250206_DiaNN)

# Generate data frame
WCP_peptide_data_20250206_DiaNN_df <- WCP_peptide_data_20250206_DiaNN %>% dplyr::select(starts_with("Intensity."))
WCP_peptide_data_20250206_DiaNN_df <- WCP_peptide_data_20250206_DiaNN_df[,order(names(WCP_peptide_data_20250206_DiaNN_df))]
names(WCP_peptide_data_20250206_DiaNN_df)

# combine relevant columns from original dataframe with adapted columns
WCP_peptide_data_20250206_DiaNN_final <- cbind(WCP_peptide_data_20250206_DiaNN[,1:3], WCP_peptide_data_20250206_DiaNN_df)
names(WCP_peptide_data_20250206_DiaNN_final)

# the output can be used for peptide maps etc.
write_tsv(WCP_peptide_data_20250206_DiaNN_final, "Resources/Proteomics/2025_02_06_Reanalysis/WCP_peptide_data_20250206_DiaNN_final_peptides.txt")

# Any duplicated sequences? Should be false
any(duplicated(WCP_peptide_data_20250206_DiaNN_final$Sequence))

# sanity checks using a matrix version of the peptides data -- log2 transformation
WCP_peptide_data_20250206_DiaNN_matrix <- WCP_peptide_data_20250206_DiaNN_final %>% 
  dplyr::select(Sequence, starts_with("Intensity")) %>% column_to_rownames("Sequence")

WCP_peptide_data_20250206_DiaNN_matrix <- as.matrix(WCP_peptide_data_20250206_DiaNN_matrix)
head(WCP_peptide_data_20250206_DiaNN_matrix)

# log-transform
WCP_peptide_data_20250206_DiaNN_matrix <- mylog(WCP_peptide_data_20250206_DiaNN_matrix)
# myboxplot(WCP_peptide_data_20250206_DiaNN_matrix) # own function

# remove outlier (12h_Depl_A) from matrix and dataframe
WCP_peptide_data_20250206_DiaNN_matrix <- WCP_peptide_data_20250206_DiaNN_matrix[,-grep("12h_IAA_A", colnames(WCP_peptide_data_20250206_DiaNN_matrix))]
# myboxplot(WCP_peptide_data_20250206_DiaNN_matrix)
WCP_peptide_data_20250206_DiaNN_final <- WCP_peptide_data_20250206_DiaNN_final[,-grep("12h_IAA_A", names(WCP_peptide_data_20250206_DiaNN_final))]

### sanity checks -- renaming of columns -- merging of linear-space and log2 data
colnames(WCP_peptide_data_20250206_DiaNN_matrix) <- sub("Intensity.", "log2_Intensity.", colnames(WCP_peptide_data_20250206_DiaNN_matrix))
WCP_peptide_data_20250206_DiaNN_matrix_df <- as.data.frame(WCP_peptide_data_20250206_DiaNN_matrix) %>% rownames_to_column("Sequence")
head(WCP_peptide_data_20250206_DiaNN_matrix_df)
names(WCP_peptide_data_20250206_DiaNN_matrix_df);names(WCP_peptide_data_20250206_DiaNN_final)
WCP_peptide_data_20250206_DiaNN_final_log2 <- merge(WCP_peptide_data_20250206_DiaNN_final, WCP_peptide_data_20250206_DiaNN_matrix_df, by = "Sequence")
dim(WCP_peptide_data_20250206_DiaNN_final);dim(WCP_peptide_data_20250206_DiaNN_matrix_df);dim(WCP_peptide_data_20250206_DiaNN_final_log2)
head(WCP_peptide_data_20250206_DiaNN_final_log2)

#### Merge with classes -------------------------------------------------------------

# merge peptide data with peptide classes (from above)
WCP_peptide_data_20250206_DiaNN_final_log2_combined <- merge(WCP_peptide_data_20250206_DiaNN_final_log2, 
                                                             combined_tbl_peptidome_class_unique_detectability, 
                                                             by.x = "Sequence", by.y = "sequence", all.x = TRUE, all.y = FALSE)
head(WCP_peptide_data_20250206_DiaNN_final_log2_combined)

# sanity check
table(WCP_peptide_data_20250206_DiaNN_final_log2_combined$peptide_class)
length(which(is.na(WCP_peptide_data_20250206_DiaNN_final_log2_combined$peptide_class))) 
# 1333/145181 (~ 1%) sequences are not assigned to a peptide class

# check which peptides have no classification 
# sanity check: should not contain NMDRHT-derived ORFs
WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% 
  filter(is.na(peptide_class)) %>% 
  mutate(source = case_when(str_detect(Protein.Group, "zzCont") ~ "zzCont",
                            str_detect(Protein.Group, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  dplyr::count(source)

dim(WCP_peptide_data_20250206_DiaNN_final_log2_combined)
any(duplicated(WCP_peptide_data_20250206_DiaNN_final_log2_combined$Sequence))
head(WCP_peptide_data_20250206_DiaNN_final_log2_combined)

# stats
WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% 
  dplyr::count(peptide_class,peptide_overJunction)

WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD <- WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>% 
  dplyr::count(NMD_peptide_reason)

# Save
WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% write_csv("Resources/Proteomics/2025_02_06_Reanalysis/WCP_peptide_data_20250206_DiaNN_final_log2_combined.csv")

#### sample-level ------------------------------------------------------------

# Pivot longer
WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand <- WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% 
  # clean_names() %>% 
  mutate(source = case_when(str_detect(Protein.Group, "zzCont") ~ "zzCont",
                            str_detect(Protein.Group, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  pivot_longer(cols = c(starts_with("Intensity."), starts_with("log2_Intensity.")),
               names_to = c(".value", "sample"),
               names_sep = "\\.") %>% 
  mutate(condition = case_when(str_detect(sample, "12h_DMSO") ~ "12h_DMSO",
                               str_detect(sample, "12h_IAA") ~ "12h_IAA",
                               str_detect(sample, "24h_DMSO") ~ "24h_DMSO",
                               str_detect(sample, "24h_IAA") ~ "24h_IAA",
                               str_detect(sample, "06h_DMSO") ~ "06h_DMSO",
                               str_detect(sample, "06h_IAA") ~ "06h_IAA",
                               str_detect(sample, "15h_DMSO") ~ "15h_DMSO",
                               str_detect(sample, "15h_IAA") ~ "15h_IAA",
                               str_detect(sample, "18h_DMSO") ~ "18h_DMSO",
                               str_detect(sample, "18h_IAA") ~ "18h_IAA"))    

# save as csv
WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand %>% write_csv("Resources/Proteomics/2025_02_06_Reanalysis/WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand.csv")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand_NMD <- WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand %>% 
  get_summary_stats(log2_Intensity)

### SN - Add Gene and tx-level info ------------------------------------------------------

##### Data sources ------------------------------------------------------------

# Read gene data - if necessary
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

# Read in - if necessary
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

# Load ORFquant P-site quantification data
load("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique.rds")

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")

# Load edgeR DTE *NMDRHT*-annotation-based data
edgeR_DTE_NMDRHT_combined <- read_csv("Resources/edgeR_DTE_NMDRHT_combined.csv")

# Load IsoformSwitchAnalyzeR (ISAR) DTU *NMDRHT*-annotation-based data
ISAR_DTU_NMDRHT_combined <- read_csv("Resources/ISAR_DTU_NMDRHT_combined.csv")

load("Resources/Proteomics/combined_peptidome_complete_class_junction.rds")

##### All Peptide-Tx-combinations ------------------------------------------------------------
# Join with non-distinct peptide table -> obtain all possible transcripts that this peptide might originate from
WCP_peptide_data_20250127_SN_final_log2_combined_all <- merge(WCP_peptide_data_20250127_SN_final_log2, 
                                                              combined_peptidome_complete_class_junction, 
                                                              by.x = "Sequence", by.y = "sequence", all.x = TRUE, all.y = FALSE)

##### Combine data ------------------------------------------------------------

# Join multiple data sources
WCP_peptide_data_20250127_SN_all_GeneTxORF <- WCP_peptide_data_20250127_SN_final_log2_combined_all %>% 
  left_join(NMDRHT.v1.2_MainTable %>% 
              dplyr::select(gene_id,
                            gene_name, 
                            transcript_id,
                            transcript_name, 
                            NMD_tx_status, 
                            NMD_tx_reason, 
                            NMD_50nt_rule,
                            DTE_cluster,
                            NMD_bin_tx)) %>% 
  # add 24h IAA DTE
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(transcript_id, logFC, FDR, logCPM) %>% 
              dplyr::rename("tx_log2FC" = "logFC",
                            "tx_FDR" = "FDR",
                            "tx_logCPM" = "logCPM")) %>% 
  # add 24h IAA DTU
  left_join(ISAR_DTU_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(isoform_id, dIF, isoform_switch_q_value) %>% 
              dplyr::rename("transcript_id" = "isoform_id",
                            "tx_ISAR_dIF" = "dIF",
                            "tx_ISAR_qval" = "isoform_switch_q_value")) %>% 
  # add ORF data
  left_join(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique %>% 
              dplyr::select(transcript_id,
                            ORF_id_tr,
                            HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN,
                            HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN,
                            HCT_N_AID_UPF1_0h_P_sites_pN,
                            HCT_N_AID_UPF1_12h_P_sites_pN,
                            HCT_N_AID_UPF1_0h_valid,
                            HCT_N_AID_UPF1_12h_valid,
                            log2FC_P_sites_pN,
                            ORFquant_NMD_tx_status_combined,
                            NMD_reason,
                            NMD_status_ORFquant) %>% 
              distinct() %>% 
              dplyr::rename("ORFquant_ORF_id_tr" = "ORF_id_tr",
                            "ORFquant_NMD_reason" = "NMD_reason",
                            "ORFquant_50nt_rule" = "NMD_status_ORFquant")) %>% 
  mutate(NMD_ORF_tx_status = case_when(!is.na(NMD_tx_status) ~ NMD_tx_status,
                                       TRUE ~ ORFquant_NMD_tx_status_combined),
         NMD_ORF_tx_reason = case_when(!is.na(NMD_tx_reason) ~ NMD_tx_status,
                                       TRUE ~ ORFquant_NMD_reason),
         NMD_ORF_50nt_rule = case_when(!is.na(NMD_50nt_rule) ~ NMD_50nt_rule,
                                       TRUE ~ ORFquant_50nt_rule)) %>% 
  # classify based on transcript NMD reason
  group_by(Sequence) %>% 
  mutate(tx_n_NMD = sum(NMD_ORF_50nt_rule == TRUE),
         tx_n_noNMD = sum(NMD_ORF_50nt_rule == FALSE)) %>% 
  mutate(n_unique_NMD_ORF_tx_reason = length(unique(NMD_ORF_tx_reason))) %>% 
  ungroup() %>% 
  # add gene data
  left_join(GENCODE_v42_MainTable %>% 
              dplyr::select(gene_id,
                            gene_name, 
                            DGE_cluster,
                            NMD_bin,
                            L2FC_kdeg,
                            Mech_conclusion) %>% 
              dplyr::rename("NMD_bin_gene" = "NMD_bin")) %>% 
  # add 24h IAA DGE
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(gene_id, log2FoldChange, padj, baseMean) %>% 
              dplyr::rename("gene_log2FC" = "log2FoldChange",
                            "gene_padj" = "padj",
                            "gene_baseMean" = "baseMean")) %>% 
  # aggregate DGE and DTE info per protein ID and 50_nt_rule
  group_by(Sequence) %>% 
  # no significance cutoffs!
  mutate(tx_log2FC_NMD_Mean = mean(tx_log2FC,na.rm = TRUE),
         .after=tx_logCPM) %>% 
  mutate(tx_dIF_NMD_Sum = sum(tx_ISAR_dIF,na.rm = TRUE),
         .after=tx_ISAR_qval) %>% 
  mutate(gene_log2FC_Mean = mean(gene_log2FC,na.rm = TRUE),
         .after=gene_baseMean) %>% 
  ungroup() 

save(WCP_peptide_data_20250127_SN_all_GeneTxORF,
     file="Resources/Proteomics/WCP_peptide_data_20250127_SN_all_GeneTxORF.rds")

### SN ANOVA F-Test ---------------------------------------------------------

# prepare Spectronaut peptides for moderated ANOVA and F-testing -- valid value filtering
# For statistical analysis, peptide intensity values were log2-transformed. 
# Peptides were filtered to retain those with at least three valid intensity measurements within each experimental group: 06h_DMSO, 06h_Depl, 12h_DMSO, 12h_Depl, 15h_DMSO, 15h_Depl, 18h_DMSO, 18h_Depl, 24h_DMSO, and 24h_Depl. 
# This filtering step resulted in the retention of approximately 64% of all peptide sequences

m <- WCP_peptide_data_20250127_SN_matrix

# m is the matrix containing log2 data
dim(m)
head(m)
any(duplicated(rownames(m)))

# valid value filter using 35% valid values
# Reason: allows DMSO samples to be emtpy and 06h IAA treatment to be empty as well
# Otherwise we exclude biologically interesting peptides that are detectable only after 12h depletion
vv <- apply(m, 1, valid.n)
hist(vv)
w <- which(vv >= 0.353*ncol(m))
length(w)/nrow(m)

filtered_data <- m[w,]

all(rownames(m) %in% rownames(filtered_data))
all(rownames(filtered_data) %in% rownames(m))

vv2 <- apply(filtered_data, 1, valid.n)
hist(vv2)
# generate the imputed matrix with downshift imputation
mf <- Impute(filtered_data)
dim(mf) / dim(m)
# myboxplot(mf)

vv3 <- apply(mf, 1, valid.n)
hist(vv3)

colnames(mf) <- paste0("Imputed_", colnames(mf))
# sanity check with a few example groups -- all should have at least 3 valid values
hist(apply(mf[,grep("06h_Depl", colnames(mf))], 1, valid.n))
hist(apply(mf[,grep("15h_DMSO", colnames(mf))], 1, valid.n))
hist(apply(mf[,grep("24h_Depl", colnames(mf))], 1, valid.n))

head(mf)
# mydensityplot(mf, show.legend = F)

# generate dataframe from filtered and imputed matrix and experimental design dataframe
filtered_dataframe <- as.data.frame(mf) %>% rownames_to_column("Sequence")
head(filtered_dataframe)
ed <- data.frame(Column.Name=make.names(names(filtered_dataframe)), Experiment=rep("", length(names(filtered_dataframe))))
ed$Experiment[grep("log2", ed$Column.Name)] <- make.names(gsub("log2_Intensity_|_[A-Z]$", "", ed$Column.Name[grep("log2", ed$Column.Name)]))

write_tsv(ed, "Resources/Proteomics/expDesign.txt")

# perform moderated F-test (modF) and moderated ANOVA (modA)
# dataframe is "filtered_dataframe" -- experimental design is "ed"

columnsForTest <- names(filtered_dataframe) %in% ed$Column.Name[ed$Experiment!=""]
class.vector   <- ed$Experiment[ed$Experiment != ""]

names(filtered_dataframe)
library(limma) # package needed for moderated statistics

# run modF
modf_test <- modF(d = filtered_dataframe, columnsForTest = columnsForTest, class.vector = class.vector, id.col = "Sequence")
head(modf_test)
modf_test$Sequence <- filtered_dataframe$Sequence

dim(filtered_dataframe);dim(mf)

# run modA
grps <- make.names(gsub("log2_Intensity_|_[A-Z]$", "", colnames(mf)))
modA <- oneWayModAnova(mat = mf, group = grps, adjust.method = "BH", sort.by = "none", trend = TRUE, repeated = FALSE, subject = NULL)
names(modA) <- paste0("ANOVA_", names(modA))
names(modA) <- make.names(names(modA))
head(modA)
modA$Sequence <- rownames(mf)

# combine both in one dataframe
names(modf_test);names(modA)
stats <- merge(modf_test, modA, by = "Sequence")
head(stats)
all(stats$modF_id == stats$Sequence)

# merge stats with original dataframe --> new dataframe r for "results"
names(filtered_dataframe);names(stats)
stats <- merge(filtered_dataframe, stats, by = "Sequence")
names(WCP_peptide_data_20250127_SN_final_log2_combined)
WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA <- merge(WCP_peptide_data_20250127_SN_final_log2_combined, stats, by = "Sequence") # results data frame

names(WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA)
any(duplicated(names(WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA)))

# export complete results including peptide_class = NA
write_csv(WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA, "Resources/Proteomics/Results_modANOVA_modF.csv")

save(WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA,
     file="Resources/Proteomics/WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA.rds")



## Conservation ---------------------------------------------------------------  

load("Resources/Proteomics/WCP_peptide_data_20250127_SN_Selection_ANOVA.rds")

# load NMDRHT cds information (GRangesList)
load("Resources/NMDRHT/NMDRHT_cds.rds")

# load ORFquant cds information (GRangesList)
load("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant_cds.rds")

##
# NMDRHT first
##

WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT <- WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>% 
              dplyr::select(Sequence, transcript_id, start, end, peptide_set, NMD_peptide_status)) %>% 
  filter(peptide_set == "NMDRHT")

# Generate IRanges
WCP_peptide_data_20250127_SN_IRanges_NMDRHT <- IRanges(start=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT$start, 
                                                       end=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT$end, 
                                                       names=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT$transcript_id,
                                                       gene_seq=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT$gene_seq,
                                                       Sequence= WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT$Sequence)

# Generate GRanges
WCP_peptide_data_20250127_SN_GRanges_NMDRHT <- GenomicFeatures::proteinToGenome(WCP_peptide_data_20250127_SN_IRanges_NMDRHT, NMDRHT_cds)
# Progate metadata as outer mcols
mcols(WCP_peptide_data_20250127_SN_GRanges_NMDRHT) <- mcols(WCP_peptide_data_20250127_SN_IRanges_NMDRHT)

##
# ORFquant second
##

WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant <- WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>% 
              dplyr::select(Sequence, transcript_id, ORFquant_ORF_id_tr, start, end, peptide_set, NMD_peptide_status)) %>% 
  filter(peptide_set == "ORFquant")

# Generate IRanges
WCP_peptide_data_20250127_SN_IRanges_ORFquant <- IRanges(start=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$start, 
                                                         end=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$end, 
                                                         names=paste0(WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$transcript_id,
                                                                      "_",
                                                                      WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$ORFquant_ORF_id_tr),
                                                         gene_seq=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$gene_seq,
                                                         Sequence= WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant$Sequence)

# Generate GRanges
WCP_peptide_data_20250127_SN_GRanges_ORFquant <- GenomicFeatures::proteinToGenome(WCP_peptide_data_20250127_SN_IRanges_ORFquant, HCT_N_AID_UPF1_ORFquant_cds)
# Progate metadata as outer mcols
mcols(WCP_peptide_data_20250127_SN_GRanges_ORFquant) <- mcols(WCP_peptide_data_20250127_SN_IRanges_ORFquant)

##
# Combine
##

WCP_peptide_data_20250127_SN_GRanges_combined <- c(WCP_peptide_data_20250127_SN_GRanges_NMDRHT,
                                                   WCP_peptide_data_20250127_SN_GRanges_ORFquant)

# Generate data frame for final merging - important: include outer mcols (contains peptide sequence ...)
WCP_peptide_data_20250127_SN_GRanges_combined_df <- as.data.frame(WCP_peptide_data_20250127_SN_GRanges_combined, use.outer.mcols=TRUE)

# For conservation score calculation - unlist
WCP_peptide_data_20250127_SN_GRanges_combined_unlist <- unlist(WCP_peptide_data_20250127_SN_GRanges_combined)

##
# phastCons100way
##
library(phastCons100way.UCSC.hg38)
phast <- phastCons100way.UCSC.hg38
class(phast)

# Get phastCons100way scores on unlisted GRanges
WCP_peptide_data_20250127_SN_GRanges_combined_phastCons <- score(phast, 
                                                                 WCP_peptide_data_20250127_SN_GRanges_combined_unlist,
                                                                 pop="DP2")

# Join for final data frame
WCP_peptide_data_20250127_SN_GRanges_combined_df$phastScore <- WCP_peptide_data_20250127_SN_GRanges_combined_phastCons

save(WCP_peptide_data_20250127_SN_GRanges_combined_df,
     file="Resources/Proteomics/WCP_peptide_data_20250127_SN_GRanges_combined_df.rds")

## Conservation Single-peptides ---------------------------------------------------------------  

##
# NMDRHT first
##

WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single <- WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>% 
              dplyr::select(Sequence, transcript_id, start, end, peptide_set, NMD_peptide_status)) %>% 
  filter(peptide_set == "NMDRHT")

# Generate IRanges
WCP_peptide_data_20250127_SN_IRanges_NMDRHT_single <- IRanges(start=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single$start, 
                                                              end=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single$end, 
                                                              names=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single$transcript_id,
                                                              gene_seq=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single$gene_seq,
                                                              Sequence= WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_NMDRHT_single$Sequence)

# Generate GRanges
WCP_peptide_data_20250127_SN_GRanges_NMDRHT_single <- GenomicFeatures::proteinToGenome(WCP_peptide_data_20250127_SN_IRanges_NMDRHT_single, NMDRHT_cds)
# Progate metadata as outer mcols
mcols(WCP_peptide_data_20250127_SN_GRanges_NMDRHT_single) <- mcols(WCP_peptide_data_20250127_SN_IRanges_NMDRHT_single)

##
# ORFquant second
##

WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single <- WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>% 
              dplyr::select(Sequence, transcript_id, ORFquant_ORF_id_tr, start, end, peptide_set, NMD_peptide_status)) %>% 
  filter(peptide_set == "ORFquant")

# Generate IRanges
WCP_peptide_data_20250127_SN_IRanges_ORFquant_single <- IRanges(start=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$start, 
                                                                end=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$end, 
                                                                names=paste0(WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$transcript_id,
                                                                             "_",
                                                                             WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$ORFquant_ORF_id_tr),
                                                                gene_seq=WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$gene_seq,
                                                                Sequence= WCP_peptide_data_20250127_SN_Selection_ANOVA_forCons_ORFquant_single$Sequence)

# Generate GRanges
WCP_peptide_data_20250127_SN_GRanges_ORFquant_single <- GenomicFeatures::proteinToGenome(WCP_peptide_data_20250127_SN_IRanges_ORFquant_single, HCT_N_AID_UPF1_ORFquant_cds)
# Progate metadata as outer mcols
mcols(WCP_peptide_data_20250127_SN_GRanges_ORFquant_single) <- mcols(WCP_peptide_data_20250127_SN_IRanges_ORFquant_single)

##
# Combine
##

WCP_peptide_data_20250127_SN_GRanges_combined_single <- c(WCP_peptide_data_20250127_SN_GRanges_NMDRHT_single,
                                                          WCP_peptide_data_20250127_SN_GRanges_ORFquant_single)

# Generate data frame for final merging - important: include outer mcols (contains peptide sequence ...)
WCP_peptide_data_20250127_SN_GRanges_combined_df_single <- as.data.frame(WCP_peptide_data_20250127_SN_GRanges_combined_single, use.outer.mcols=TRUE)

# For conservation score calculation - unlist
WCP_peptide_data_20250127_SN_GRanges_combined_unlist_single <- unlist(WCP_peptide_data_20250127_SN_GRanges_combined_single)

##
# phastCons100way
##
library(phastCons100way.UCSC.hg38)
phast <- phastCons100way.UCSC.hg38
class(phast)

# Get phastCons100way scores on unlisted GRanges
WCP_peptide_data_20250127_SN_GRanges_combined_phastCons_single <- score(phast, 
                                                                        WCP_peptide_data_20250127_SN_GRanges_combined_unlist_single,
                                                                        pop="DP2")

# Join for final data frame
WCP_peptide_data_20250127_SN_GRanges_combined_df_single$phastScore <- WCP_peptide_data_20250127_SN_GRanges_combined_phastCons_single

save(WCP_peptide_data_20250127_SN_GRanges_combined_df_single,
     file="Resources/Proteomics/WCP_peptide_data_20250127_SN_GRanges_combined_df_single.rds")

##
# Protein-Group-Level -> Data preparation (2025-01-27) --------------------------------------------------------
##

# Load proteomics quantification (from 2025-01-17) -> clean names
WCP_PG_data_20250127 <- read_excel("Resources/Proteomics/2025_01_27_Reanalysis/UPF1_NMDRHT_WCP_protein_group_SN.xlsx") %>% 
  janitor::clean_names()

# How many protein groups
WCP_PG_data_20250127 %>% 
  dplyr::count()
## 9030 quantified proteins

# How many IDs per source 
WCP_PG_data_20250127 %>% 
  # define source of protein group information
  mutate(source = case_when(str_detect(pg_protein_groups, "zzCont") ~ "zzCont",
                            str_detect(pg_protein_groups, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  dplyr::count(source)
## 8317 NMDRHT 
## 658 UniProt
## 55 zzCont

# Which columns are present?
colnames(WCP_PG_data_20250127)

# Problem: some protein groups match multiple database entries (e.g. multiple IDs in UniProt and/or NMDRHT-derived ORFs)
# Approach: Unlist protein groups -> filter for NMDRHT or (if not found) UniProt IDs -> filter out zzCont (common contaminations)
WCP_PG_data_20250127_unlist <- WCP_PG_data_20250127 %>% 
  # Filter out potential contaminants - even if they could contain NMDRHT information!
  filter(str_detect(pg_protein_groups, "zzCont") == FALSE) %>% 
  # expand informative columns
  separate_rows(c("pg_protein_groups", "pg_genes", "pg_organisms", "pg_fasta_files", "pg_molecular_weight"), sep =";", convert = TRUE) %>% 
  # select only basal columns (no quantification yet)
  dplyr::select(id, genes, uniprot, starts_with("pg")) %>% 
  dplyr::select(-c(pg_protein_descriptions, pg_uni_prot_ids, pg_protein_names))  %>% 
  # define source of protein group information
  mutate(source = case_when(str_detect(pg_protein_groups, "zzCont") ~ "zzCont",
                            str_detect(pg_protein_groups, "NMDRHT") ~ "NMDRHT",
                            TRUE ~ "UniProt")) %>% 
  # exclude contaminations
  filter(source != "zzCont") %>% 
  # filter for NMDRHT-defined protein groups or only-UniProt-defined
  group_by(id) %>% 
  mutate(NMDRHT_count = sum(source == "NMDRHT")) %>% 
  ungroup() %>% 
  filter((NMDRHT_count > 0 & source == "NMDRHT") | (NMDRHT_count == 0)) %>% 
  distinct(id, pg_protein_groups, .keep_all = TRUE) %>% 
  group_by(id) %>% 
  mutate(pg_count = sum(!is.na(source))) %>% 
  ungroup() %>% 
  dplyr::select(-NMDRHT_count)

# How many distinct ids - relation to source
WCP_PG_data_20250127_unlist %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::count(source) 
# Numbers match with initial data

# Further tidy up the data and include quantifications
WCP_PG_data_20250127_unlist_quant <- WCP_PG_data_20250127_unlist %>% 
  # Fix problem with PAR_Y starting protein groups
  mutate(pg_protein_groups = str_remove(pg_protein_groups, "PAR_Y_")) %>% 
  # recover information from NMDRHT ORF description
  separate(pg_protein_groups, c("ORF_transcript_id", "start", "end", "transcript_id"),remove = FALSE, extra = "merge") %>% 
  mutate(transcript_id = case_when(source == "NMDRHT" ~ transcript_id,
                                   TRUE ~ NA)) %>% 
  mutate(ORFquant_ORF_id_tr = case_when(source == "NMDRHT" ~ paste0(ORF_transcript_id, "_", start, "_", end))) %>% 
  dplyr::rename("UniProtKB" = "uniprot") %>% 
  dplyr::mutate(UniProtKB = case_when(source == "UniProt" ~ UniProtKB,
                                      TRUE ~ NA)) %>% 
  dplyr::select(-c(ORF_transcript_id, start, end)) %>% 
  relocate(transcript_id, ORFquant_ORF_id_tr, .after = UniProtKB) %>% 
  left_join(WCP_PG_data_20250127 %>% dplyr::select(id,
                                                   starts_with("log_fc"),
                                                   starts_with("adj_p")))

# Sanity Check: How many distinct ids - relation to source
WCP_PG_data_20250127_unlist_quant %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::count(source) 

# How many protein-Tx-ORF matchings total
WCP_PG_data_20250127_unlist_quant %>% 
  dplyr::count(source)

# How many ORF ids per protein id
WCP_PG_data_20250127_unlist_quant %>% 
  # filter(source != "UniProt") %>% 
  distinct(id, ORFquant_ORF_id_tr, .keep_all = TRUE) %>% 
  dplyr::count(source)
# 10821 distinct NMDRHT id-to-ORFid matches

# Distribution of unique/multiple hits per ID
WCP_PG_data_20250127_unlist_quant %>% 
  dplyr::count(source, pg_count) %>% 
  ggplot() +
  geom_density(aes(x=pg_count,
                   fill=source),
               alpha=0.5)

# Verdict: many detected protein groups match to multiple ORFs (high number of pg_count)

###
#### Add gene and tx-level info ----------------------------------------------
###

##### Data sources ------------------------------------------------------------

# Read gene data - if necessary
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

# Read in - if necessary
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

# Load ORFquant P-site quantification data
load("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique.rds")

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")

# Load edgeR DTE *NMDRHT*-annotation-based data
edgeR_DTE_NMDRHT_combined <- read_csv("Resources/edgeR_DTE_NMDRHT_combined.csv")

# Load IsoformSwitchAnalyzeR (ISAR) DTU *NMDRHT*-annotation-based data
ISAR_DTU_NMDRHT_combined <- read_csv("Resources/ISAR_DTU_NMDRHT_combined.csv")


##### Combine data ------------------------------------------------------------

# Join multiple data sources
WCP_PG_data_20250127_GeneTxORF <- WCP_PG_data_20250127_unlist_quant %>% 
  # Filter out UniProt
  # filter(source != "UniProt") %>% 
  # add transcript data
  left_join(NMDRHT.v1.2_MainTable %>% 
              dplyr::select(gene_id,
                            gene_name, 
                            transcript_id,
                            transcript_name, 
                            NMD_tx_status, 
                            NMD_tx_reason, 
                            NMD_50nt_rule,
                            DTE_cluster,
                            NMD_bin_tx)) %>% 
  # add 24h IAA DTE
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(transcript_id, logFC, FDR, logCPM) %>% 
              dplyr::rename("tx_log2FC" = "logFC",
                            "tx_FDR" = "FDR",
                            "tx_logCPM" = "logCPM")) %>% 
  # add 24h IAA DTU
  left_join(ISAR_DTU_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(isoform_id, dIF, isoform_switch_q_value) %>% 
              dplyr::rename("transcript_id" = "isoform_id",
                            "tx_ISAR_dIF" = "dIF",
                            "tx_ISAR_qval" = "isoform_switch_q_value")) %>% 
  # add ORF data
  left_join(HCT_N_AID_UPF1_ORFquant_RiboTish_NMD_sumORF_all_unique %>% 
              dplyr::select(transcript_id,
                            ORF_id_tr,
                            HCT_N_AID_UPF1_0h_tx_P_sites_sum_pN,
                            HCT_N_AID_UPF1_12h_tx_P_sites_sum_pN,
                            HCT_N_AID_UPF1_0h_P_sites_pN,
                            HCT_N_AID_UPF1_12h_P_sites_pN,
                            HCT_N_AID_UPF1_0h_valid,
                            HCT_N_AID_UPF1_12h_valid,
                            log2FC_P_sites_pN,
                            ORFquant_NMD_tx_status_combined,
                            NMD_reason,
                            NMD_status_ORFquant) %>% 
              distinct() %>% 
              dplyr::rename("ORFquant_ORF_id_tr" = "ORF_id_tr",
                            "ORFquant_NMD_reason" = "NMD_reason",
                            "ORFquant_50nt_rule" = "NMD_status_ORFquant")) %>% 
  mutate(NMD_ORF_tx_status = case_when(!is.na(NMD_tx_status) ~ NMD_tx_status,
                                       TRUE ~ ORFquant_NMD_tx_status_combined),
         NMD_ORF_tx_reason = case_when(!is.na(NMD_tx_reason) ~ NMD_tx_status,
                                       TRUE ~ ORFquant_NMD_reason),
         NMD_ORF_50nt_rule = case_when(!is.na(NMD_50nt_rule) ~ NMD_50nt_rule,
                                       TRUE ~ ORFquant_50nt_rule)) %>% 
  # classify based on transcript NMD reason
  group_by(id) %>% 
  mutate(tx_n_NMD = sum(NMD_ORF_50nt_rule == TRUE),
         tx_n_noNMD = sum(NMD_ORF_50nt_rule == FALSE)) %>% 
  mutate(n_unique_NMD_ORF_tx_reason = length(unique(NMD_ORF_tx_reason))) %>% 
  ungroup() %>% 
  mutate(protein_class = case_when(tx_n_NMD > 0 & tx_n_noNMD > 0 ~ "mixed",
                                   tx_n_NMD > 0 & tx_n_noNMD == 0  ~ "NMD",
                                   tx_n_NMD == 0 & tx_n_noNMD > 0 ~ "coding",
                                   source == "UniProt" ~ "UniProt")) %>% 
  mutate(NMD_protein_status = case_when(n_unique_NMD_ORF_tx_reason == 1 ~ "unique",
                                        n_unique_NMD_ORF_tx_reason > 1 ~ "multiple",
                                        source == "UniProt" ~ "UniProt")) %>% 
  mutate(NMD_protein_reason = case_when(NMD_protein_status == "unique" ~ NMD_ORF_tx_reason,
                                        NMD_protein_status == "multiple" ~ "multiple",
                                        source == "UniProt" ~ "UniProt")) %>% 
  # add gene data
  left_join(GENCODE_v42_MainTable %>% 
              dplyr::select(gene_id,
                            gene_name, 
                            DGE_cluster,
                            NMD_bin,
                            L2FC_kdeg,
                            L2FC_ksyn,
                            Mech_conclusion) %>% 
              dplyr::rename("NMD_bin_gene" = "NMD_bin")) %>% 
  # add 24h IAA DGE
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_24h")) %>% 
              dplyr::select(gene_id, log2FoldChange, padj, baseMean) %>% 
              dplyr::rename("gene_log2FC" = "log2FoldChange",
                            "gene_padj" = "padj",
                            "gene_baseMean" = "baseMean")) %>% 
  # aggregate DGE and DTE info per protein ID and 50_nt_rule
  group_by(id) %>% 
  mutate(tx_log2FC_NMD_sigMean = mean(tx_log2FC[NMD_ORF_50nt_rule == TRUE & tx_FDR < 0.0001]),
         tx_log2FC_coding_sigMean = mean(tx_log2FC[NMD_ORF_50nt_rule == FALSE & tx_FDR < 0.0001]),
         .after=tx_logCPM) %>% 
  mutate(tx_dIF_NMD_sigSum = sum(tx_ISAR_dIF[NMD_ORF_50nt_rule == TRUE & tx_ISAR_qval < 0.0001]),
         tx_dIF_coding_sigSum = sum(tx_ISAR_dIF[NMD_ORF_50nt_rule == FALSE & tx_ISAR_qval < 0.0001]),
         .after=tx_ISAR_qval) %>% 
  mutate(gene_log2FC_sigMean = mean(gene_log2FC[gene_padj < 0.0001]),
         .after=gene_baseMean) %>% 
  ungroup()

# How many ORF ids per protein id
WCP_PG_data_20250127_GeneTxORF %>% 
  distinct(id, ORFquant_ORF_id_tr, .keep_all = TRUE) %>% 
  dplyr::count(source)
# Same number as above (10821)

# Sanity Check: How many distinct ids - relation to source
WCP_PG_data_20250127_GeneTxORF %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::count(source) 

##### Long Format -------------------------------------------------------------

# Pivot to long format
WCP_PG_data_20250127_GeneTxORF_long <- WCP_PG_data_20250127_GeneTxORF %>% 
  dplyr::rename("UPF1_06h_MS_log2FC" = "log_fc_depl_06h_over_dmso_06h",
                "UPF1_12h_MS_log2FC" = "log_fc_depl_12h_over_dmso_12h",
                "UPF1_15h_MS_log2FC" = "log_fc_depl_15h_over_dmso_15h",
                "UPF1_18h_MS_log2FC" = "log_fc_depl_18h_over_dmso_18h",
                "UPF1_24h_MS_log2FC" = "log_fc_depl_24h_over_dmso_24h",
                "UPF1_06h_MS_padj" = "adj_p_val_depl_06h_over_dmso_06h",
                "UPF1_12h_MS_padj" = "adj_p_val_depl_12h_over_dmso_12h",
                "UPF1_15h_MS_padj" = "adj_p_val_depl_15h_over_dmso_15h",
                "UPF1_18h_MS_padj" = "adj_p_val_depl_18h_over_dmso_18h",
                "UPF1_24h_MS_padj" = "adj_p_val_depl_24h_over_dmso_24h") %>% 
  pivot_longer(cols = c(starts_with("UPF1_")),
               names_to = c("comparison",".value"),
               names_sep = "_MS_")

save(WCP_PG_data_20250127_GeneTxORF_long ,
     file="Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long.rds")

# Filter for significant hits
WCP_PG_data_20250127_GeneTxORF_long_sig <- WCP_PG_data_20250127_GeneTxORF_long %>% 
  filter(padj < 0.001 & abs(log2FC) > 0.5)

WCP_PG_data_20250127_GeneTxORF_long_sig_simple <- WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
  dplyr::select(gene_name, comparison, log2FC, padj)

##### Overall GO significant proteins----------------------------------------

###### Up ----------------------------------------------------------------------

# Prepare data
WCP_PG_data_20250127_GeneTxORF_long_sig_up_forGO <- WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
  filter(log2FC > 0) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)


# Make lists
WCP_PG_data_20250127_GeneTxORF_long_sig_up_forGO_list <- split(WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
                                                                 filter(log2FC > 0) %>% 
                                                                 distinct(gene_id, comparison) %>% 
                                                                 separate(gene_id, c("gene_id", "Version")) %>% 
                                                                 pull(gene_id),
                                                               WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
                                                                 filter(log2FC > 0) %>% 
                                                                 distinct(gene_id, comparison) %>% 
                                                                 separate(gene_id, c("gene_id", "Version")) %>% 
                                                                 pull(comparison))

# Perform GO analysis
WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres <- gprofiler2::gost(query = WCP_PG_data_20250127_GeneTxORF_long_sig_up_forGO_list,
                                                           custom_bg = Rev_1_F2_GO_DGE_bg_full,
                                                           # multi_query = TRUE,
                                                           sources = "GO:BP",
                                                           domain_scope = "custom",
                                                           organism = "hsapiens",
                                                           correction_method = c("gSCS"),
                                                           evcodes = TRUE,
                                                           # as_short_link = TRUE,
                                                           significant = TRUE)

# As tibble
WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres_result <- WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres$result %>% 
  arrange(p_value)

###### Enrichment visualization ----------------------------------------


##
# sim Matrix
##
simMatrix_WCP_PG_data_up_combined <- rrvgo::calculateSimMatrix(WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres_result %>% 
                                                          distinct(term_id) %>% 
                                                          pull(term_id),
                                                        orgdb="org.Hs.eg.db",
                                                        ont="BP",
                                                        method="Rel")

##
# scores
##

scores_WCP_PG_data_up_combined <- setNames(-log10(WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres_result %>% 
                                                    distinct(term_id, .keep_all = TRUE) %>% 
                                                    pull(p_value)), 
                                           WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres_result %>% 
                                             distinct(term_id) %>% 
                                             pull(term_id))

##
# reduced terms
##

reducedTerms_WCP_PG_data_up_combined <- rrvgo::reduceSimMatrix(simMatrix_WCP_PG_data_up_combined,
                                                        scores_WCP_PG_data_up_combined,
                                                        threshold=0.95,
                                                        orgdb="org.Hs.eg.db")

# Combine G:profiler output with combined&reduced GO IDs
WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO <- WCP_PG_data_20250127_GeneTxORF_long_sig_up_gostres_result %>% 
  left_join(reducedTerms_WCP_PG_data_up_combined %>% 
              dplyr::select(go, parent, parentTerm) %>% 
              dplyr::rename("term_id" = "go",
                            "parent_id" = "parent")) %>% 
  filter(!is.na(parent_id)) %>% 
  group_by(parentTerm,query) %>% 
  mutate(group_score = mean(-log10(p_value))) %>% 
  ungroup()

#save as rds
save(WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO,
     file="Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO.rds")

###### Down ----------------------------------------------------------------------

# Prepare data

WCP_PG_data_20250127_GeneTxORF_long_sig_down_forGO <- WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
  filter(log2FC < 0) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)


# Make lists
WCP_PG_data_20250127_GeneTxORF_long_sig_down_forGO_list <- split(WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
                                                                   filter(log2FC < 0) %>% 
                                                                   distinct(gene_id, comparison) %>% 
                                                                   separate(gene_id, c("gene_id", "Version")) %>% 
                                                                   pull(gene_id),
                                                                 WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
                                                                   filter(log2FC < 0) %>% 
                                                                   distinct(gene_id, comparison) %>% 
                                                                   separate(gene_id, c("gene_id", "Version")) %>% 
                                                                   pull(comparison))

# Perform GO analysis
WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres <- gprofiler2::gost(query = WCP_PG_data_20250127_GeneTxORF_long_sig_down_forGO_list,
                                                             custom_bg = Rev_1_F2_GO_DGE_bg_full,
                                                             # multi_query = TRUE,
                                                             sources = "GO:BP",
                                                             domain_scope = "custom",
                                                             organism = "hsapiens",
                                                             correction_method = c("gSCS"),
                                                             evcodes = TRUE,
                                                             # as_short_link = TRUE,
                                                             significant = TRUE)

# As tibble
WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres_result <- WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres$result %>% 
  arrange(p_value)

###### Enrichment visualization ----------------------------------------

##
# sim Matrix
##
simMatrix_WCP_PG_data_down_combined <- rrvgo::calculateSimMatrix(WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres_result %>% 
                                                            distinct(term_id) %>% 
                                                            pull(term_id),
                                                          orgdb="org.Hs.eg.db",
                                                          ont="BP",
                                                          method="Rel")

##
# scores
##

scores_WCP_PG_data_down_combined <- setNames(-log10(WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres_result %>% 
                                                      distinct(term_id, .keep_all = TRUE) %>% 
                                                      pull(p_value)), 
                                             WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres_result %>% 
                                               distinct(term_id) %>% 
                                               pull(term_id))

##
# reduced terms
##

reducedTerms_WCP_PG_data_down_combined <- rrvgo::reduceSimMatrix(simMatrix_WCP_PG_data_down_combined,
                                                          scores_WCP_PG_data_down_combined,
                                                          threshold=0.9,
                                                          orgdb="org.Hs.eg.db")

# Combine G:profiler output with combined&reduced GO IDs
WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO <- WCP_PG_data_20250127_GeneTxORF_long_sig_down_gostres_result %>% 
  left_join(reducedTerms_WCP_PG_data_down_combined %>% 
              dplyr::select(go, parent, parentTerm) %>% 
              dplyr::rename("term_id" = "go",
                            "parent_id" = "parent")) %>% 
  filter(!is.na(parent_id)) %>% 
  group_by(parentTerm,query) %>% 
  mutate(group_score = sum(p_value)) %>% 
  ungroup()

#save as rds
save(WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO,
     file="Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO.rds")