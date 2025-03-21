#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Annotation
# Objective: Generation of NMD-Regulated Human Transcriptome (NMDRHT) via transcriptome assembly
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)
library(ggh4x)
library(ggridges)
library(GenomicFeatures)

##
# GTF analyses ------------------------------------------------------------
##

# Comment: prepare simplified version of GENCODE version 42 annotation - contains informative tags!

## Prepared only once - GENCODE v42 simplify

gencode.v42.annotation <- rtracklayer::readGFF("Resources/GENCODE/gencode.v42.annotation.gff3.gz",
                                           filter = list(type=c("gene",
                                                                "transcript")))
gencode.v42.annotation <- as_tibble(gencode.v42.annotation)

gencode.v42.annotation_tags <- gencode.v42.annotation %>%
  mutate(ensembl_canonical = case_when(str_detect(tag, "Ensembl_canonical") ~ TRUE,
                                                  TRUE ~ FALSE)) %>%
  mutate(ensembl_basic = case_when(str_detect(tag, "basic") ~ TRUE,
                                       TRUE ~ FALSE)) %>%
  mutate(appris_principal = case_when(str_detect(tag, "appris_principal") ~ TRUE,
                                   TRUE ~ FALSE)) 
  

gtf_gencode_df_short <- gencode.v42.annotation_tags %>%
  dplyr::select("type",
                "gene_id",
                "gene_name",
                "gene_type",
                "transcript_id",
                "transcript_type",
                "transcript_name",
                "transcript_support_level",
                "ensembl_canonical",
                "ensembl_basic",
                "appris_principal") %>%
  distinct()

gtf_gencode_df_short %>% write_csv("Resources/GENCODE/gencode.v42.gtf_df_short.csv")

gtf_gencode_df_short <- read_csv("Resources/GENCODE/gencode.v42.gtf_df_short.csv")

##
# General statistics on GENCODE v42
##

# How many distinct gene_ids?
gtf_gencode_df_short %>% 
  filter(type == "gene") %>% 
  distinct(gene_id, .keep_all = TRUE) %>% 
  dplyr::select(-type) %>% 
  dplyr::count()

# How many distinct transcript_ids?
gtf_gencode_df_short %>% 
  filter(type == "transcript") %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::select(-type) %>% 
  dplyr::count()

#
##
## LR&SR complete approach ---------------------------------------------------------------
##
#

# Comment: LR = long read | SR = short read RNA-Seq data
# Comment: the individual transcriptome assemblies and SQANTI3 analyses were previously generated with individual scripts 
# Comment: the output files (SQANTI3 classification.txt and GTF file) need to be present in the indicated SQANTI_path and GTF_path!
# Comment: alternatively, the checkpoint .rds files could be loaded as intermediate steps

##
### LR complete  ----------------------------------------------------------------
##

# Import LR dataset metadata
NMDRHT_LR_datasets <- read_csv("Resources/Longread/NMDRHT_LR_complete_datasets.csv", 
                                                trim_ws = TRUE) %>% 
  mutate(unique_ID = fct_inorder(unique_ID)) %>% 
  mutate(method = fct_inorder(method)) %>% 
  mutate(sequencing = fct_inorder(sequencing)) %>% 
  mutate(annotation = fct_inorder(annotation)) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(study = "HCT_N_AID_UPF1_LR") %>% 
  mutate(folder = case_when(sequencing == "PacBio" ~ "2024_05_UPF1_PacBio",
                            sequencing == "ONT_dRNA" ~ "2024_06_UPF1_Nanopore")) 

##
### SR combined-only  ----------------------------------------------------------------
##
# Import SR dataset metadata
NMDRHT_SR_datasets <- read_csv("Resources/Shortread/NMDRHT_SR_complete_datasets.csv", 
                                         trim_ws = TRUE) %>% 
  mutate(unique_ID = fct_inorder(unique_ID)) %>% 
  mutate(method = fct_inorder(method)) %>% 
  mutate(sequencing = fct_inorder(sequencing)) %>% 
  mutate(annotation = fct_inorder(annotation)) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(study = fct_inorder(study)) %>% 
  mutate(folder = fct_inorder(folder)) %>% 
  filter(condition_2 == "combined")

##
### LR&SR complete  ----------------------------------------------------------------
##
NMDRHT_datasets <- NMDRHT_LR_datasets %>% 
  bind_rows(NMDRHT_SR_datasets)

##
### Export GTF path ----------------------------------------------------------------
##
# for easier run of gffcompare
NMDRHT_datasets %>% 
  pull(GTF_path) %>% 
  readr::write_lines(file="Resources/NMDRHT/gffcompare/NMDRHT_GTF_list.txt")

##
### Import SQANTI3 output (classification.txt) ----------------------------------------------------------------
##

# Define chromomes to keep
selected_chromosomes = c(paste0("chr",seq_len(22)), "chrX", "chrY", "chrM")

# import SQANTI3 classifications, rename chromosome and isoform columns and filter out all-NA columns
NMDRHT_SQANTI3 <- NMDRHT_datasets %>% 
  left_join(readr::read_tsv(NMDRHT_datasets %>% 
                              pull(SQANTI_path), 
                            id = "SQANTI_path")) %>% 
  dplyr::rename("GTF_transcript_id" = "isoform",
         "seqname" = "chrom") %>%
  dplyr::filter(seqname %in% selected_chromosomes) %>% 
  select_if(~sum(!is.na(.)) > 0)

# How many assembled transcripts overall?
NMDRHT_SQANTI3 %>% 
  dplyr::count()

# How many assembled transcripts per structural_category?
NMDRHT_SQANTI3 %>% 
  dplyr::count(structural_category) %>% 
  arrange(desc(n))

##
### Merge NMDRegHumanTxome_SQANTI3 & GENCODE ----------------------------------------------------------------
##
# Rename GENCODE-based columns
NMDRHT_SQANTI3_GENCODE <- NMDRHT_SQANTI3 %>%
  left_join(gtf_gencode_df_short %>% 
              filter(type == "gene") %>% 
              dplyr::select(gene_id, gene_name, gene_type),
            by=c("associated_gene" = "gene_id")) %>% 
  left_join(gtf_gencode_df_short %>% 
              filter(type == "transcript") %>% 
              dplyr::select(transcript_id, transcript_type, transcript_name, transcript_support_level, ensembl_canonical, ensembl_basic, appris_principal),
            by=c("associated_transcript" = "transcript_id")) %>% 
  dplyr::rename("GENCODE_gene_id" = "associated_gene",
         "GENCODE_gene_name" = "gene_name",
         "GENCODE_gene_type" = "gene_type",
         "GENCODE_transcript_id" = "associated_transcript",
         "GENCODE_transcript_name" = "transcript_name",
         "GENCODE_transcript_type" = "transcript_type",
         "GENCODE_transcript_support_level" = "transcript_support_level",) %>% 
  relocate(GENCODE_gene_id, GENCODE_gene_name, GENCODE_gene_type,
           GENCODE_transcript_id, GENCODE_transcript_name, GENCODE_transcript_type, GENCODE_transcript_support_level,
           ensembl_canonical, ensembl_basic, appris_principal,
           .after=GTF_path)

##
### Import GTF files  ----------------------------------------------------------------
##
# for extracting genomic start/end coordinates of each transcripts
NMDRHT_fromGTF <- NMDRHT_datasets %>% 
  left_join(readr::read_tsv(NMDRHT_datasets %>% 
                              pull(GTF_path),
                            col_names = c("seqname",
                                          "source",
                                          "feature",
                                          "start",
                                          "end",
                                          "score",
                                          "strand",
                                          "frame",
                                          "attribute"),
                            col_select = c("seqname",
                                           "source",
                                           "feature",
                                           "start",
                                           "end",
                                           "strand",
                                           "attribute"),
                            id = "GTF_path") %>% 
              filter(feature == "transcript") %>% 
              separate(attribute, 
                       into = c("GTF_transcript_id",NA, NA), 
                       sep = "\\;",
                       extra = "drop",
                       fill = "right") %>% 
              separate(GTF_transcript_id, 
                       into = c(NA,"GTF_transcript_id", NA), 
                       sep = "\\\"",
                       extra = "drop",
                       fill = "right")
            )

##
### Join GTF & SQANTI3  ----------------------------------------------------------------
##
# Match GTF_transcript_ids from SQANTI3 output with start&end information from GTF files 
NMDRHT_SQANTI3_GENCODE_fromGTF <- NMDRHT_SQANTI3_GENCODE %>% 
  left_join(NMDRHT_fromGTF %>% 
              dplyr::select(unique_ID,start,end,GTF_transcript_id)) %>% 
  relocate(start, end,
           .after=seqname)

##
#### Checkpoint #1 - save data --------------------------------------------------------------
##

save(NMDRHT_datasets,
     NMDRHT_SQANTI3_GENCODE,
     NMDRHT_fromGTF, 
     NMDRHT_SQANTI3_GENCODE_fromGTF,
     file = paste0("Resources/NMDRHT/",Sys.Date(),"_NMDRHT_checkpoint_1_datasources.rds"))

# Load if necessary
# load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_1_datasources.rds")

##
# Run Gffcompare -------------------------------------------------
##
# Obtain unique intron chains (UIC)
# This part uses external scripts and tools like IGV tools - this has to be configured properly!

# Run gffcompare script
command_gffcompare_NMDRHT <- paste(c("bash ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/gffcompare/NMDRHT_gffcompare.sh"), sep="", collapse = "")

run_gffcompare_NMDRHT <- system(command_gffcompare_NMDRHT, intern = TRUE)

# Run sort and index
command_sort_NMDRegHumanTxome <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/gffcompare/NMDRHT.combined.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/gffcompare/NMDRHT.combined.sort.gtf"), sep="", collapse = "")

run_sort_NMDRegHumanTxome <- system(command_sort_NMDRegHumanTxome, intern = TRUE)

command_index_NMDRegHumanTxome <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                         "index ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/gffcompare/NMDRHT.combined.sort.gtf"), sep="", collapse = "")

run_index_NMDRegHumanTxome <- system(command_index_NMDRegHumanTxome, intern = TRUE)

##
## Gffcompare tracking  -------------------------------------------------
##
# The tracking file gives information about which UIC was found in which individual transcriptome
NMDRHT_tracking <- read_delim("/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/gffcompare/NMDRHT.tracking", 
                                      delim = "\t", 
                                      escape_double = FALSE,
                                      trim_ws = TRUE,
                                      na = c("", "NA", "-"),
                                      col_names = c("NMDRHT_transcript_id",
                                                    "NMDRHT_locus_id",
                                                    "reference_id",
                                                    "class_code",
                                                    NMDRHT_datasets %>% 
                                                      mutate(unique_ID = as.character(unique_ID)) %>% 
                                                      dplyr::pull(unique_ID)))

## Clean columns  -------------------------------------------------
NMDRHT_tracking_cleaned <- NMDRHT_tracking %>% 
  separate(NMDRHT_transcript_id, c("NMDRHT_transcript_id",NA, NA), sep = "\\|") %>% 
  mutate(NMDRHT_transcript_id = str_replace(NMDRHT_transcript_id, "_", ""),
         NMDRHT_locus_id = str_replace(NMDRHT_locus_id, "XLOC", "NMDRHL")) %>% 
  mutate(NMDRHT_locus_id = str_replace(NMDRHT_locus_id, "_", "")) %>% 
  pivot_longer(cols=-c("NMDRHT_transcript_id",
                       "NMDRHT_locus_id",
                       "reference_id",
                       "class_code"),
               names_to = "query",
               values_to = "query_transcript_id") %>% 
  separate(query_transcript_id, c(NA, "query_transcript_id", NA), sep = "\\|", extra = "drop") %>% 
  pivot_wider(names_from = "query",
              values_from = "query_transcript_id") %>% 
  mutate(class_code_simple = case_when(class_code %in% c("=",
                                                         "k",
                                                         "j") ~ class_code,
                                       TRUE ~ "other"),
         .after = class_code) %>% 
  separate(reference_id, c("gffcompare_gene_id", "gffcompare_transcript_id"), sep = "\\|", extra = "drop")

## Calculate overlap metrics  -------------------------------------------------
NMDRHT_tracking_cleaned_sum <- NMDRHT_tracking_cleaned %>% 
  mutate(UIC_total_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>% 
                                                     mutate(unique_ID = as.character(unique_ID)) %>% 
                                                     dplyr::pull(unique_ID))))),
         UIC_LR_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>% 
                                                             filter(sequencing != "Illumina") %>%
                                                             mutate(unique_ID = as.character(unique_ID)) %>% 
                                                             dplyr::pull(unique_ID))))),
         UIC_LR_Bambu_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>% 
                                                         filter(method == "Bambu") %>%
                                                         mutate(unique_ID = as.character(unique_ID)) %>% 
                                                         dplyr::pull(unique_ID))))),
         UIC_LR_IsoQuant_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>% 
                                                            filter(method == "IsoQuant") %>%
                                                            mutate(unique_ID = as.character(unique_ID)) %>% 
                                                            dplyr::pull(unique_ID))))),
         UIC_LR_StringTie2_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                              filter(method == "StringTie2" & sequencing != "Illumina") %>%
                                                              mutate(unique_ID = as.character(unique_ID)) %>%
                                                              dplyr::pull(unique_ID))))),
         UIC_SR_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                               filter(sequencing == "Illumina") %>%
                                                               mutate(unique_ID = as.character(unique_ID)) %>%
                                                               dplyr::pull(unique_ID))))),
         UIC_ONT_dRNA_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                       filter(sequencing == "ONT_dRNA") %>%
                                                       mutate(unique_ID = as.character(unique_ID)) %>%
                                                       dplyr::pull(unique_ID))))),
         UIC_PacBio_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                       filter(sequencing == "PacBio") %>%
                                                       mutate(unique_ID = as.character(unique_ID)) %>%
                                                       dplyr::pull(unique_ID))))),
         UIC_HEK293_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                             filter(cell_line == "HEK293") %>%
                                                             mutate(unique_ID = as.character(unique_ID)) %>%
                                                             dplyr::pull(unique_ID))))),
         UIC_HCT116_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                           filter(cell_line == "HCT116") %>%
                                                           mutate(unique_ID = as.character(unique_ID)) %>%
                                                           dplyr::pull(unique_ID))))),
         UIC_HFF_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                                  filter(cell_line == "HFF1") %>%
                                                                  mutate(unique_ID = as.character(unique_ID)) %>%
                                                                  dplyr::pull(unique_ID))))),
         UIC_HUVEC_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                                  filter(cell_line == "HUVEC") %>%
                                                                  mutate(unique_ID = as.character(unique_ID)) %>%
                                                                  dplyr::pull(unique_ID))))),
         UIC_deNovo_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                           filter(annotation == "deNovo") %>%
                                                           mutate(unique_ID = as.character(unique_ID)) %>%
                                                           dplyr::pull(unique_ID))))),
         UIC_GENCODE_support = rowSums(!is.na(dplyr::select(., c(NMDRHT_datasets %>%
                                                           filter(annotation == "GENCODE") %>%
                                                           mutate(unique_ID = as.character(unique_ID)) %>%
                                                           dplyr::pull(unique_ID))))),
         .after = "class_code_simple") %>% 
  left_join(gtf_gencode_df_short %>% 
              dplyr::rename("gffcompare_gene_id" = "gene_id",
                            "gffcompare_gene_name" = "gene_name",
                            "gffcompare_gene_type" = "gene_type",
                            "gffcompare_transcript_id" = "transcript_id",
                            "gffcompare_transcript_name" = "transcript_name",
                            "gffcompare_transcript_type" = "transcript_type")) %>% 
  relocate(gffcompare_gene_id, gffcompare_gene_name, gffcompare_gene_type, gffcompare_transcript_id, gffcompare_transcript_name, gffcompare_transcript_type) 

### Stats -------------------------------------------------------------------

# Total number of UICs
NMDRHT_tracking_cleaned_sum %>% 
  distinct(NMDRHT_transcript_id) %>% 
  dplyr::count(name = "total_NMDRHT_transcripts")

##
# Join Gffcompare & SQANTI3  ---------------------------------------------------
##

##
## Prepare gffcompare output for join  ---------------------------------------------------
##
NMDRHT_tracking_forJoin <- NMDRHT_tracking_cleaned_sum %>% 
  filter(UIC_total_support >= 2) %>% 
  dplyr::select(NMDRHT_transcript_id, 
              NMDRHT_locus_id,
              starts_with("gffcompare"),
              class_code,
              class_code_simple,
              ends_with("_support"),
              NMDRHT_datasets %>% 
                mutate(unique_ID = as.character(unique_ID)) %>% 
                dplyr::pull(unique_ID)) %>% 
  pivot_longer(cols=NMDRHT_datasets %>% 
                 mutate(unique_ID = as.character(unique_ID)) %>% 
                 dplyr::pull(unique_ID),
               names_to = "unique_ID",
               values_to = "GTF_transcript_id") %>% 
  filter(!is.na(GTF_transcript_id))

##
## Join ---------------------------------------------------
##
# Note that gffcompare and SQANTI3 do not always match the same transcript to the same gene -> thus we need to handle conflicts properly
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge <- NMDRHT_SQANTI3_GENCODE_fromGTF %>% 
  right_join(NMDRHT_tracking_forJoin,
            by = c("unique_ID" = "unique_ID",
                   "GTF_transcript_id" = "GTF_transcript_id")) %>% 
  dplyr::filter(!is.na(seqname)) %>% 
  mutate(gene_id_conflict = case_when(gffcompare_gene_id == GENCODE_gene_id ~ "exact_match",
                                      str_detect(GENCODE_gene_id, "novelGene") ~ "novel_gene",
                                      str_detect(GENCODE_gene_id, gffcompare_gene_id) ~ "contained",
                                      structural_category == "fusion" ~ "fusion",
                                      str_count(GENCODE_gene_id, "ENSG") > 1 ~ "SQANTI_not_sure",
                                      TRUE ~ "conflict")) %>% 
  mutate(gene_id = case_when(gene_id_conflict == "exact_match" ~ GENCODE_gene_id,
                             gene_id_conflict == "novel_gene" ~ GENCODE_gene_id,
                             gene_id_conflict == "contained" ~ gffcompare_gene_id,
                             gene_id_conflict == "fusion" ~ gffcompare_gene_id,
                             gene_id_conflict == "SQANTI_not_sure" ~ gffcompare_gene_id,
                             gene_id_conflict == "conflict" ~ GENCODE_gene_id),
         gene_name = case_when(gene_id_conflict == "exact_match" ~ GENCODE_gene_name,
                               gene_id_conflict == "novel_gene" ~ GENCODE_gene_id,
                               gene_id_conflict == "contained" ~ gffcompare_gene_name,
                               gene_id_conflict == "fusion" ~ gffcompare_gene_name,
                               gene_id_conflict == "SQANTI_not_sure" ~ gffcompare_gene_name,
                               gene_id_conflict == "conflict" ~ GENCODE_gene_name)) %>% 
  relocate(NMDRHT_transcript_id,
           NMDRHT_locus_id,
           gene_id,
           gene_name,
           gene_id_conflict,
           starts_with("gffcompare"),
           starts_with("GENCODE"))

# Distribution of conflicts between gffcompare and SQANTI3
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge %>% 
  dplyr::count(gene_id_conflict)

# How many gene_ids are left
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge %>% 
  distinct(gene_id) %>% 
  dplyr::count()

# How many transcript_ids per conflict type
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge %>% 
  distinct(NMDRHT_transcript_id, gene_id_conflict) %>%  
  dplyr::count(gene_id_conflict)

# First select "exact_match" transcripts
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_exactMatch <-  NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge %>% 
  filter(gene_id_conflict == "exact_match")

# Second round: fill with transcripts that do not have exact match
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_conflict <- NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge %>% 
filter(!NMDRHT_transcript_id %in% NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_exactMatch$NMDRHT_transcript_id)

NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total <- NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_exactMatch %>% 
  bind_rows(NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_exactMatch) 

# Check problematic gene_ids
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
  group_by(gene_id) %>% 
  filter(sum(strand == "+") > 1 & sum(strand == "-") > 1)
  
# How many gene_ids are left
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
  distinct(gene_id) %>% 
  dplyr::count()

# How many transcripts are left
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
  distinct(NMDRHT_transcript_id) %>% 
  dplyr::count()

# Check that each transcript is matched with only one gene_id
NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
  distinct(NMDRHT_transcript_id, gene_id) %>% 
  dplyr::count(NMDRHT_transcript_id) %>% 
  dplyr::count(n == 1)

### Different filters ---------------------------------------------------
#### Gene_id-wise *mean* UIC_support filter ---------------------------------------------------
# Require at least UIC_support 4 to consider transcripts for high-support

NMDRHT_tracking_HighSup <- NMDRHT_tracking_cleaned_sum %>% 
  right_join(NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
               dplyr::distinct(NMDRHT_transcript_id,
                        NMDRHT_locus_id,
                        gene_id,
                        gene_name,
                        gene_id_conflict)) %>% 
  filter(UIC_total_support >= 4) %>% 
  group_by(NMDRHT_locus_id) %>% 
  mutate(mean_UIC_support = mean(UIC_total_support)) %>% 
  ungroup() %>%
  filter(UIC_total_support >= mean_UIC_support) %>% 
  mutate(UIC_filter_level = "high")

# How many loci are left?
NMDRHT_tracking_HighSup %>% 
  distinct(gene_id) %>% 
  dplyr::count()

# How many transcripts are left?
NMDRHT_tracking_HighSup %>% 
  distinct(NMDRHT_transcript_id) %>% 
  dplyr::count()

# Add loci with low support
NMDRHT_tracking_filtered <- NMDRHT_tracking_HighSup %>% 
  bind_rows(NMDRHT_tracking_cleaned_sum %>% 
              right_join(NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total %>% 
                           dplyr::distinct(NMDRHT_transcript_id,
                                           NMDRHT_locus_id,
                                           gene_id,
                                           gene_name,
                                           gene_id_conflict)) %>% 
              filter(!gene_id %in% NMDRHT_tracking_HighSup$gene_id) %>% 
              filter(UIC_total_support >= 2) %>% 
              group_by(gene_id) %>% 
              mutate(mean_UIC_support = mean(UIC_total_support)) %>% 
              ungroup() %>%
              filter(UIC_total_support >= mean_UIC_support) %>% 
              mutate(UIC_filter_level = "low"))

# How many loci are left?
NMDRHT_tracking_filtered %>% 
  distinct(gene_id) %>% 
  dplyr::count()

# How many transcripts are left?
NMDRHT_tracking_filtered %>% 
  distinct(NMDRHT_transcript_id) %>% 
  dplyr::count()

##
#### Checkpoint #2 - save data --------------------------------------------------------------
##

save(NMDRHT_tracking_cleaned_sum,
     NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge,
     NMDRHT_SQANTI3_GENCODE_fromGTF_gffcompare_merge_total, 
     NMDRHT_tracking_filtered,
     file = paste0("Resources/NMDRHT/",Sys.Date(),"_NMDRHT_checkpoint_2_datasources.rds"))

# Load if necessary
# load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_2_datasources.rds")

##
# Final-Join Gffcompare & SQANTI3  ---------------------------------------------------
##
# Use filtered UICs and supplement with rich information from SQANTI3

##
## Prepare gffcompare filtered output for join  ---------------------------------------------------
##
NMDRHT_tracking_filtered_ids <- NMDRHT_tracking_filtered %>% 
  dplyr::select(NMDRHT_transcript_id, 
                NMDRHT_locus_id,
                gene_id,
                gene_name,
                gene_id_conflict,
                class_code,
                class_code_simple,
                ends_with("_support"),
                NMDRHT_datasets %>% 
                  mutate(unique_ID = as.character(unique_ID)) %>% 
                  dplyr::pull(unique_ID)) %>% 
  pivot_longer(cols=NMDRHT_datasets %>% 
                 mutate(unique_ID = as.character(unique_ID)) %>% 
                 dplyr::pull(unique_ID),
               names_to = "unique_ID",
               values_to = "GTF_transcript_id") %>% 
  filter(!is.na(GTF_transcript_id))
  
##
## Join ---------------------------------------------------
##
NMDRHT_SQANTI3_GENCODE_fromGTF_Filtered <-NMDRHT_SQANTI3_GENCODE_fromGTF %>% 
  filter(paste0(unique_ID,GTF_transcript_id) %in% (NMDRHT_tracking_filtered_ids %>% 
                                           mutate(unique_ID_GTF_transcript_id = paste0(unique_ID,GTF_transcript_id)) %>% 
                                           pull(unique_ID_GTF_transcript_id)
                                         )) %>% 
  left_join(NMDRHT_tracking_filtered_ids,
            by = c("unique_ID" = "unique_ID",
                   "GTF_transcript_id" = "GTF_transcript_id")) 
 

##
## Median start/end  ---------------------------------------------------
##
# Slice best matching GTF_transcript_id based on distance to median start&end
# First only deNovo & having at least 1 LR_support
NMDRHT_final_selection_part1 <- NMDRHT_SQANTI3_GENCODE_fromGTF_Filtered  %>% 
  filter(annotation == "deNovo" & UIC_LR_support > 0) %>% 
  group_by(NMDRHT_transcript_id) %>% 
  mutate(median_start = median(start),
         median_end = median(end),
         .after = end) %>% 
  ungroup() %>% 
  mutate(diff_to_medianStart = median_start-start,
         diff_to_medianEnd = median_end-end,
         .after = median_end) %>% 
  group_by(NMDRHT_transcript_id) %>%
  dplyr::slice(which.min(abs(diff_to_medianStart)+abs(diff_to_medianEnd))) %>% 
  ungroup()

# Second only deNovo & having no LR_support
NMDRHT_final_selection_part2 <- NMDRHT_SQANTI3_GENCODE_fromGTF_Filtered  %>% 
  filter(!NMDRHT_transcript_id %in% NMDRHT_final_selection_part1$NMDRHT_transcript_id) %>% 
  filter(annotation == "deNovo") %>% 
  group_by(NMDRHT_transcript_id) %>% 
  mutate(median_start = median(start),
         median_end = median(end),
         .after = end) %>% 
  ungroup() %>% 
  mutate(diff_to_medianStart = median_start-start,
         diff_to_medianEnd = median_end-end,
         .after = median_end) %>% 
  group_by(NMDRHT_transcript_id) %>%
  dplyr::slice(which.min(abs(diff_to_medianStart)+abs(diff_to_medianEnd))) %>% 
  ungroup()

# Third also those only GENODE-based
NMDRHT_final_selection_part3 <- NMDRHT_SQANTI3_GENCODE_fromGTF_Filtered  %>% 
  filter(!NMDRHT_transcript_id %in% c(NMDRHT_final_selection_part1$NMDRHT_transcript_id,
                                      NMDRHT_final_selection_part2$NMDRHT_transcript_id)) %>% 
  group_by(NMDRHT_transcript_id) %>% 
  mutate(median_start = median(start),
         median_end = median(end),
         .after = end) %>% 
  ungroup() %>% 
  mutate(diff_to_medianStart = median_start-start,
         diff_to_medianEnd = median_end-end,
         .after = median_end) %>% 
  group_by(NMDRHT_transcript_id) %>%
  dplyr::slice(which.min(abs(diff_to_medianStart)+abs(diff_to_medianEnd))) %>% 
  ungroup()

NMDRHT_final_selection <- NMDRHT_final_selection_part1 %>% 
  bind_rows(NMDRHT_final_selection_part2,
            NMDRHT_final_selection_part3)
  
# How many transcripts per annotation
NMDRHT_final_selection %>% 
  distinct(NMDRHT_transcript_id, annotation, .keep_all = TRUE) %>% 
  dplyr::count(annotation)

# Check if no problematic gene_ids exist (containing transcripts on both strands?) -> should be none
NMDRHT_final_selection %>% 
group_by(gene_id) %>% 
  filter(sum(strand == "+") > 1 & sum(strand == "-") > 1)


##
#### Checkpoint #3 - save data --------------------------------------------------------------
##

save(NMDRHT_SQANTI3_GENCODE_fromGTF_Filtered,
     NMDRHT_final_selection,
     file = paste0("Resources/NMDRHT/",Sys.Date(),"_NMDRHT_checkpoint_3_datasources.rds"))

# Load if necessary
# load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_3_datasources.rds")

##
# Extract from GTF --------------------------------------------------------
##

##
### Import full (transcripts & exons) GTF files  ----------------------------------------------------------------
##
# for generating final annotation file
# WARNING: this might take a long time!
NMDRHT_fromGTF_full <- NMDRHT_datasets %>% 
  right_join(readr::read_tsv(NMDRHT_datasets %>%
                               pull(GTF_path),
                             col_names = c("seqnames",
                                          "source",
                                          "type",
                                          "start",
                                          "end",
                                          "score",
                                          "strand",
                                          "phase",
                                          "attribute"),
                            id = "GTF_path") %>% 
              separate(attribute, 
                       into = c("GTF_transcript_id",NA,NA), 
                       sep = "\\;",
                       extra = "drop",
                       fill = "right") %>% 
              separate(GTF_transcript_id, 
                       into = c(NA,"GTF_transcript_id", NA), 
                       sep = "\\\"",
                       extra = "drop",
                       fill = "right")
  )

# Filter full GTF file for selected transcripts
NMDRHT_fromGTF_full_filtered <- NMDRHT_fromGTF_full %>% 
  mutate(unique_ID_GTF_transcript_id = paste0(unique_ID,GTF_transcript_id)) %>% 
  filter(unique_ID_GTF_transcript_id %in% (NMDRHT_final_selection %>% 
                                             mutate(unique_ID_GTF_transcript_id = paste0(unique_ID,GTF_transcript_id)) %>% 
                                             pull(unique_ID_GTF_transcript_id))) %>% 
  dplyr::select(-unique_ID_GTF_transcript_id)

# Join and extract relevant information - reorder and finish additional tags
NMDRHT_final_selection_forGTFexport <- NMDRHT_fromGTF_full_filtered %>% 
  full_join(NMDRHT_final_selection %>% 
              dplyr::select(NMDRHT_transcript_id,
                            NMDRHT_locus_id,
                            GTF_transcript_id,
                            gene_id,
                            gene_name,
                            gene_id_conflict,
                            length,
                            unique_ID,
                            starts_with("GENCODE"),
                            ensembl_canonical,
                            structural_category,
                            subcategory,
                            UIC_total_support,
                            UIC_LR_support,
                            UIC_SR_support,
                            class_code,
                            class_code_simple,
                            within_CAGE_peak,
                            within_polyA_site,
                            polyA_motif_found
                     )) %>% 
  dplyr::rename("transcript_id" = "NMDRHT_transcript_id",
         "width" = "length") %>% 
  dplyr::select(-c(SQANTI_path,
                   GTF_path,
                   condition_2,
                   unique_ID)) %>% 
  relocate(cell_line,method,sequencing,annotation, .after=class_code) %>% 
  mutate(structural_category_simple = case_when(structural_category == "full-splice_match" ~ "FSM",
                                                structural_category == "incomplete-splice_match" ~ "ISM",
                                                structural_category == "novel_not_in_catalog" ~ "NNIC",
                                                structural_category == "novel_in_catalog" ~ "NIC",
                                                structural_category == "intergenic" ~ "IG",
                                                structural_category == "fusion" ~ "F",
                                                structural_category == "antisense" ~ "AS",
                                                structural_category == "genic" ~ "G"
                                                ),
         .after = structural_category)  %>% 
  mutate(transcript_name = case_when(!is.na(GENCODE_transcript_name) ~ paste0(GENCODE_transcript_name,"_",structural_category_simple),
                                     TRUE ~ paste0(transcript_id,"_",structural_category_simple))) %>% 
  mutate(gene_name = case_when(!is.na(gene_name) ~ gene_name,
                               TRUE ~ gene_id)) %>% 
  relocate(seqnames,
           start,
           end,
           width,
           strand,
           source,
           type,
           score,
           phase,
           gene_id,
           transcript_id,
           gene_name,
           transcript_name,
           GENCODE_gene_id,
           GENCODE_transcript_id,
           GENCODE_transcript_name,
           GENCODE_gene_type,
           GENCODE_transcript_type,
           GENCODE_transcript_support_level,
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
           polyA_motif_found
           ) %>% 
  filter(!is.na(start)) %>% 
  dplyr::select(-c(folder)) %>% 
  mutate(across(.cols = c(GENCODE_gene_id,
                          GENCODE_transcript_id,
                          GENCODE_transcript_name,
                          GENCODE_gene_type,
                          GENCODE_transcript_type,
                          GENCODE_transcript_support_level,
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
                          annotation
                          ),
                .fns = ~case_when(type == "exon" ~ NA,
                                  TRUE ~ (.))))

# Export as gtf file
rtracklayer::export(GenomicRanges::makeGRangesFromDataFrame(NMDRHT_final_selection_forGTFexport,
  keep.extra.columns=TRUE),
  "Resources/NMDRHT/NMDRHT_final_selection.gtf")

# Run sort and index
command_sort <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT_final_selection.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT_final_selection.sort.gtf"), sep="", collapse = "")
output_sort <- system(command_sort, intern = TRUE)

command_index <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "index ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT_final_selection.sort.gtf"), sep="", collapse = "")
output_index <- system(command_index, intern = TRUE)

#
##
# ORF annotation ----------------------------------------------------------------
##
#
# Comment: We want to use the ensembl-canonical transcript for each gene as basis for ORF prediction/rescue
# Comment: Otherwise, many problematic transcripts with strange ORFs are used as basis, which rather complicates things

#
## Prepare GENCODE information ---------------------------------------------
#

# For export to IGV-web - modify to include only canonical transcripts and remove unnecessary information - include only transcripts (drastically reduces file size)
gencode.v42.annotation_forIGV <- rtracklayer::readGFF("Resources/GENCODE/gencode.v42.annotation.gff3.gz") %>% 
  mutate(ensembl_canonical = case_when(str_detect(tag, "Ensembl_canonical") ~ TRUE,
                                       TRUE ~ FALSE)) %>%
  mutate(ensembl_basic = case_when(str_detect(tag, "basic") ~ TRUE,
                                   TRUE ~ FALSE)) %>%
  mutate(appris_principal = case_when(str_detect(tag, "appris_principal") ~ TRUE,
                                      TRUE ~ FALSE)) %>% 
  filter(ensembl_canonical == TRUE) %>% 
  filter(type=="transcript") %>% 
  separate(gene_id, c("gene_id_plain", NA), remove = FALSE) %>% 
  separate(transcript_id, c("transcript_id_plain", NA), remove = FALSE) %>% 
  dplyr::select(-c(artif_dupl,ccdsid,protein_id,ont,havana_gene,hgnc_id,exon_id,exon_number,havana_transcript))

# Export canonical-only for IGV-prepared GENCODE GFF3 file - run once
rtracklayer::export(gencode.v42.annotation_forIGV %>% filter(type=="transcript") , "Resources/GENCODE/gencode.v42.annotation.canonical_forIGV_tx.gtf")

# Load Gencode Annotation - run once
gencode.v42.annotation <- rtracklayer::readGFF("Resources/GENCODE/gencode.v42.annotation.gtf")

# Export **for IGV** canonical-only GENCODE GTF file - only transcripts -> plain IDs allowed as well
rtracklayer::export(gencode.v42.annotation %>%
                      separate(gene_id, c("gene_id_plain", NA), remove = FALSE) %>% 
                      separate(transcript_id, c("transcript_id_plain", NA), remove = FALSE) %>% 
                      filter(transcript_id %in% c(gtf_gencode_df_short %>% filter(ensembl_canonical == TRUE) %>% pull(transcript_id))) , "Resources/GENCODE/gencode.v42.annotation.canonical_forIGV.gtf")


# Export canonical-only GENCODE GTF file - run once
rtracklayer::export(gencode.v42.annotation %>%
                     filter(transcript_id %in% c(gtf_gencode_df_short %>% filter(ensembl_canonical == TRUE) %>% pull(transcript_id))) , "Resources/GENCODE/gencode.v42.annotation.canonical.gtf")

# Export canonical- & NMD-only GENCODE GTF file - run once
rtracklayer::export(gencode.v42.annotation %>%
                     filter(transcript_id %in% c(gtf_gencode_df_short %>% filter(ensembl_canonical == TRUE | transcript_type == "nonsense_mediated_decay") %>% pull(transcript_id))) , "Resources/GENCODE/gencode.v42.annotation.canonical_NMD.gtf")

# Get list of Ensembl canonical transcript ids
GENCODE_canonical <- gtf_gencode_df_short %>% filter(ensembl_canonical == TRUE) %>% pull(transcript_id)

# Import *canonical* annotation as TxDb
gencode.v42.annotation_canonical_txDb <- GenomicFeatures::makeTxDbFromGFF("Resources/GENCODE/gencode.v42.annotation.canonical.gtf")

# Get *canonical* CDS per transcript
gencode.v42.annotation_canonical_txDb_cds <- cdsBy(gencode.v42.annotation_canonical_txDb, by="tx", use.names=TRUE)

# Save CDS information
save(gencode.v42.annotation_canonical_txDb_cds,
     file = paste0("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds.rds"))

# Obtain ranges of CDS
gencode.v42.annotation_canonical_txDb_cds_ranges <- unlist( range(gencode.v42.annotation_canonical_txDb_cds) )

# Restrict to start codon position
gencode.v42.annotation_canonical_txDb_start_codons <- resize(gencode.v42.annotation_canonical_txDb_cds_ranges, 1, fix="start")

# Restrict to start codon position
gencode.v42.annotation_canonical_txDb_stop_codons <- resize(gencode.v42.annotation_canonical_txDb_cds_ranges, 1, fix="end")

# Annotate the start codons and obtain tibble
gencode.v42.annotation_canonical_txDb_start_stop_codons_tbl <- as_tibble(gencode.v42.annotation_canonical_txDb_start_codons) %>% 
  mutate(transcript_id = names(gencode.v42.annotation_canonical_txDb_start_codons)) %>% 
  dplyr::select(-end) %>% 
  dplyr::rename("start_codon_canon" = "start") %>% 
  left_join(as_tibble(gencode.v42.annotation_canonical_txDb_stop_codons) %>% 
              mutate(transcript_id = names(gencode.v42.annotation_canonical_txDb_stop_codons)) %>% 
              dplyr::select(-start) %>% 
              dplyr::rename("stop_codon_canon" = "end")) %>% 
  relocate(stop_codon_canon, .after=start_codon_canon) %>% 
  left_join(gtf_gencode_df_short %>% 
              filter(type=="transcript") %>% 
              dplyr::select(transcript_id,
                            transcript_name,
                            gene_id,
                            gene_name))

# Get genome fasta file - Required to obtain nucleotide sequence
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Obtain GENCODE canonical CDS sequence
gencode.v42.annotation_canonical_txDb_cds_seqs <- extractTranscriptSeqs(Hsapiens, gencode.v42.annotation_canonical_txDb_cds)

## A sanity check:
stopifnot(identical(width(gencode.v42.annotation_canonical_txDb_cds_seqs), unname(sum(width(gencode.v42.annotation_canonical_txDb_cds)))))

gencode.v42.annotation_canonical_txDb_cds_AA <- Biostrings::translate(gencode.v42.annotation_canonical_txDb_cds_seqs, if.fuzzy.codon="solve")

gencode.v42.annotation_canonical_txDb_cds_AA_df <- as.data.frame(gencode.v42.annotation_canonical_txDb_cds_AA) %>% 
  rownames_to_column(var="GENCODE_transcript_id") %>% 
  dplyr::rename("sequence" = "x") %>% 
  # remove stop codon from width calculation
  mutate(width = str_length(sequence)-1) %>% 
  left_join(gencode.v42.annotation_canonical_txDb_start_stop_codons_tbl %>% 
              dplyr::rename("GENCODE_transcript_id" = "transcript_id") %>% 
              dplyr::select(-width))

### Save information to rds -------------------------------------------
save(gencode.v42.annotation_canonical_txDb_cds_AA_df,
     file = paste0("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds_AA_df.rds"))

# prepare cleaned GENCODE v42 canonical protein sequence for export/saving
gencode.v42.annotation_canonical_txDb_cds_AA_df_forExport <- gencode.v42.annotation_canonical_txDb_cds_AA_df %>%
  dplyr::select(GENCODE_transcript_id, sequence) %>% 
  mutate(sequence = stringr::str_replace(sequence, '\\*', '')) %>% 
  dplyr::rename("transcript_id" = "GENCODE_transcript_id",
                "Protein" = "sequence") %>% 
  left_join(gencode.v42.annotation_canonical_txDb_start_stop_codons_tbl %>% dplyr::select(transcript_id, transcript_name, gene_id, gene_name)) %>% 
  relocate(gene_id, gene_name, transcript_id, transcript_name)

save(gencode.v42.annotation_canonical_txDb_cds_AA_df_forExport,
     file = paste0("Resources/GENCODE/gencode.v42.annotation_canonical_txDb_cds_AA_df_forExport.rds"))

#
## factR - domains ---------------------------------------------------------
#

library(factR)

# Load canonical transcript only GENCODE GTF file
gencode.v42.annotation.canonical.gtf_for_factR <- factR::importGTF("Resources/GENCODE/gencode.v42.annotation.canonical.gtf")

# Functions from factR
gencode.v42.annotation_canonical_txDb_cds_cleaned <- .extractCDSchecks(gencode.v42.annotation.canonical.gtf_for_factR, Hsapiens)

gencode.v42.annotation_canonical_txDb_cds_AA_cleaned <- .getSequence(gencode.v42.annotation_canonical_txDb_cds_cleaned, Hsapiens) %>% 
  filter(noATG == FALSE & instop == FALSE) %>% 
  mutate(x = stringr::str_replace(x, '\\*', ''))

gencode.v42.annotation_canonical_txDb_cds_AA_cleaned_fasta <- c(rbind(paste0(">", gencode.v42.annotation_canonical_txDb_cds_AA_cleaned$id), gencode.v42.annotation_canonical_txDb_cds_AA_cleaned$x))

write(x = gencode.v42.annotation_canonical_txDb_cds_AA_cleaned_fasta, file = "Resources/GENCODE/gencode.v42.annotation_canonical.fa")  

#
##
# Run ORFanage  ----------------------------------------------------------------
##
#
# Run ORFanage locally with --mode ALL
# Comment: this allows picking the best matching ORF

# Run gffcompare script
command_gffcompare_ORFanage <- paste(c("bash ",
                                       "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/ORFanage.sh"), sep="", collapse = "")

run_gffcompare_ORFanage <- system(command_gffcompare_ORFanage, intern = TRUE)

# Run sort and index
command_sort_ORFanage <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                                         "sort ",
                                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_all.gtf ",
                                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_all.sort.gtf"), sep="", collapse = "")

run_sort_ORFanage <- system(command_sort_ORFanage, intern = TRUE)

command_index_ORFanage <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                                          "index ",
                                  "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_all.sort.gtf"), sep="", collapse = "")

run_index_ORFanage <- system(command_index_ORFanage, intern = TRUE)

# Import ORFanage stats
NMDRHT_ORFanage_stats_all <- read_delim("Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_stats_all.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE) %>%  
  separate(query_id, 
           into = c("query_id", "ORF_id"), 
           sep = "\\.",
           extra = "drop",
           fill = "right") %>% 
  left_join(NMDRHT_final_selection,
            by = c("query_id" = "NMDRHT_transcript_id")) %>% 
  separate(ORF_id,
           into = c("start", "end"),
           remove = FALSE) %>% 
  mutate(start_codon = case_when(strand == "+" ~ start,
                                 strand == "-" ~ end),
         stop_codon = case_when(strand == "+" ~ end,
                                strand == "-" ~ start),
         .after = end)

# Match ORFanage with canonical start codon positions
NMDRHT_ORFanage_stats_all_canonicalMatch <- NMDRHT_ORFanage_stats_all %>% 
  dplyr::select(-c(start,end)) %>% 
  left_join(gencode.v42.annotation_canonical_txDb_start_stop_codons_tbl %>% 
              dplyr::select(seqnames,start_codon_canon,stop_codon_canon,strand,gene_id) %>% 
              dplyr::rename("GENCODE_gene_id" = "gene_id",
                            "seqname" = "seqnames")) %>% 
  relocate(start_codon_canon, .after = start_codon) %>% 
  relocate(stop_codon_canon, .after = stop_codon)

# Classify ORFanage-annotated ORFs
NMDRHT_ORFanage_filt <- NMDRHT_ORFanage_stats_all_canonicalMatch %>% 
  filter(!notes %in% c("missing_start",
                       "missing_stop",
                       "no-overlap",
                       "dup")) %>% 
  mutate(start_codon_type = case_when(is.na(start_codon) ~ "non_coding",
                                      start_codon == start_codon_canon ~ "start_canonical",
                                      TRUE ~ "start_other"),
         stop_codon_type = case_when(is.na(stop_codon) ~ "non_coding",
                                     stop_codon == stop_codon_canon ~ "stop_canonical",
                                      TRUE ~ "stop_other"),
         .after=stop_codon_canon)  %>% 
  mutate(ORF_type = case_when(start_codon_type == "start_canonical" & stop_codon_type == "stop_canonical" ~ "exact_match",
                              start_codon_type == "start_canonical" ~ "start_match",
                              start_codon_type == "stop_canonical" ~ "stop_match",
                              start_codon_type == "start_other" ~ "start_other",
                              start_codon_type == "stop_other" ~ "stop_other",
                              TRUE ~ "non_coding"),
         .after=stop_codon_type) 


# Collapse annotated ORFs - by classified ORF
NMDRHT_ORFanage_collapsed <- NMDRHT_ORFanage_filt %>% 
  group_by(query_id) %>% 
  mutate(query_ORF_type = case_when(sum(ORF_type == "exact_match") > 0 ~ "exact_match",
                                    sum(ORF_type == "start_match") > 0 ~ "start_match",
                                    sum(ORF_type == "stop_match") > 0 ~ "stop_match",
                                    sum(ORF_type == "start_other") > 0 ~ "start_other",
                                    sum(ORF_type == "stop_other") > 0 ~ "stop_other",
                                    TRUE ~ "non_coding"),
         .after=ORF_type) %>% 
  ungroup()

NMDRHT_ORFanage_collapsed_selected <- NMDRHT_ORFanage_collapsed %>% 
  filter(query_ORF_type == "exact_match" & start_codon_type == "start_canonical" & notes != "-") %>% 
  bind_rows(NMDRHT_ORFanage_collapsed %>% 
              filter(query_ORF_type == "start_match" & start_codon_type == "start_canonical" & notes != "-")) %>% 
  bind_rows(NMDRHT_ORFanage_collapsed %>% 
              filter(query_ORF_type == "stop_match" & stop_codon_type == "stop_canonical" & notes != "-")) %>% 
  bind_rows(NMDRHT_ORFanage_collapsed %>% 
              filter(query_ORF_type == "start_other" & start_codon_type == "start_other" & notes == "best_gtf")) %>% 
  bind_rows(NMDRHT_ORFanage_collapsed %>% 
              filter(query_ORF_type == "stop_other" & stop_codon_type == "start_other" & notes == "best_gtf")) %>% 
  bind_rows(NMDRHT_ORFanage_collapsed %>% 
              filter(query_ORF_type == "non_coding" & start_codon_type == "non_coding")) %>% 
  mutate(transcript_id = case_when(query_ORF_type == "non_coding" ~ paste0(query_id),
                                   TRUE ~ paste0(query_id,".",ORF_id)),
         .after=ORF_id)


# Check if duplicate isoforms exist - should contain no transcripts
NMDRHT_ORFanage_collapsed_selected %>% 
  group_by(query_id) %>% 
  filter(n()>1) %>% 
  dplyr::count()

# Problematic transcripts with strand-problems exist?
NMDRHT_ORFanage_collapsed_selected %>% 
  group_by(gene_id) %>% 
  filter(sum(strand == "+") > 1 & sum(strand == "-") > 1)

NMDRHT_ORFanage_collapsed_selected %>% 
  dplyr::count(ORF_type)

# Load ORFanage Annotation

NMDRHT_ORFanage_all.sort.gtf <- rtracklayer::readGFF("Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_all.sort.gtf") %>% 
  relocate(seqid, start, end, strand,
           source, type, score, phase,
           gene_id,
           transcript_id,
           gene_name,
           transcript_name,
           GENCODE_gene_id,
           GENCODE_transcript_id,
           GENCODE_transcript_name,
           GENCODE_gene_type,
           GENCODE_transcript_type,
           GENCODE_transcript_support_level,
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
           within_CAGE_peak, within_polyA_site, polyA_motif_found,
           orfanage_status, orfanage_duplicity, orfanage_template, orfanage_template_source,
           study, GTF_transcript_id, cell_line, method, sequencing, annotation
           )

# Recover lost gene_id for exons
NMDRHT_ORFanage_all.sort.gtf_fixed <- NMDRHT_ORFanage_all.sort.gtf %>%
  filter(transcript_id %in% NMDRHT_ORFanage_collapsed_selected$transcript_id) %>% 
  left_join(NMDRHT_ORFanage_collapsed_selected %>% dplyr::select(transcript_id,
                                                                  ORF_type,
                                                                  start_codon_type,
                                                                  stop_codon_type)) %>% 
  separate(transcript_id, 
           into = c("transcript_id", "ORF_id"), 
           sep = "\\.",
           extra = "drop",
           fill = "right") %>% 
  arrange(transcript_id) %>% 
  fill(gene_id, gene_name, transcript_name) %>% 
  mutate(across(.cols = c(ORF_type,
                          start_codon_type,
                          stop_codon_type
  ),
  .fns = ~case_when(type %in% c("exon", "CDS") ~ NA,
                    TRUE ~ (.))))

# Export collapsed ORFanage annotation
rtracklayer::export(NMDRHT_ORFanage_all.sort.gtf_fixed,
                    "Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_collapsed.gtf")

# Run sort and index
command_sort <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_collapsed.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_collapsed.sort.gtf"), sep="", collapse = "")
output_sort <- system(command_sort, intern = TRUE)

command_index <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                         "index ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_collapsed.sort.gtf"), sep="", collapse = "")
output_index <- system(command_index, intern = TRUE)

# Check again for entries without gene_id
NMDRHT_ORFanage_all.sort.gtf_fixed  %>%
  filter(is.na(gene_id))

# Check for problematic entries
NMDRHT_ORFanage_all.sort.gtf_fixed %>% 
  group_by(gene_id) %>% 
  filter(sum(strand == "+") > 1 & sum(strand == "-") > 1)

#
##
# FactR -------------------------------------------------------------------
##
#
# Use to predict NMD sensitivity of ORFanage'd transcripts

library(factR) 

NMDRHT_ORFanage_all.sort.gtf_fixed_for_factR <- factR::importGTF("Resources/NMDRHT/ORFanage/NMDRHT_ORFanage_collapsed.sort.gtf") 

NMDRHT_ORFanage_NMDprediction.out <- factR::predictNMD(NMDRHT_ORFanage_all.sort.gtf_fixed_for_factR)  

NMDRHT_ORFanage_NMDprediction.out_annotated <- NMDRHT_ORFanage_NMDprediction.out %>% 
  dplyr::rename("transcript_id" = "transcript") %>% 
  left_join(NMDRHT_ORFanage_all.sort.gtf_fixed %>% 
              dplyr::filter(type == "transcript") %>%
              dplyr::select(gene_id, gene_name, transcript_id, transcript_name, GENCODE_transcript_type, structural_category_simple)) %>% 
  relocate(gene_id, gene_name, transcript_id, transcript_name, GENCODE_transcript_type, structural_category_simple)

# How many NMD-sensitive transcripts
NMDRHT_ORFanage_NMDprediction.out_annotated %>% 
  dplyr::count(is_NMD) %>% 
  arrange(desc(n))

NMDRHT_ORFanage_all.sort.gtf_fixed_NMD <- NMDRHT_ORFanage_all.sort.gtf_fixed %>%
  left_join(NMDRHT_ORFanage_NMDprediction.out_annotated %>% 
              dplyr::select(transcript_id,
                            is_NMD,
                            stop_to_lastEJ,
                            num_of_downEJs,
                            "3'UTR_length"
                            ) %>% 
              dplyr::rename("UTR3_length" = "3'UTR_length",
                            "NMD_status" = "is_NMD")) %>% 
  mutate(gene_biotype = case_when(is.na(GENCODE_gene_type) & ORF_type != "non_coding" ~ "protein_coding",
                                  is.na(GENCODE_gene_type) & ORF_type == "non_coding" ~ "lncRNA",
                                  !is.na(GENCODE_gene_type) ~ GENCODE_gene_type),
         transcript_biotype = case_when(is.na(GENCODE_transcript_type) & ORF_type != "non_coding" & NMD_status == FALSE ~ "protein_coding",
                                        is.na(GENCODE_transcript_type) & ORF_type != "non_coding" & NMD_status == TRUE ~ "nonsense_mediated_decay",
                                        is.na(GENCODE_transcript_type) & ORF_type == "non_coding" ~ "lncRNA",
                                        ORF_type != "non_coding" & NMD_status == FALSE ~ "protein_coding",
                                        ORF_type != "non_coding" & NMD_status == TRUE ~ "nonsense_mediated_decay",
                                        ORF_type == "non_coding" ~ "lncRNA"),
         .after=transcript_name) %>% 
  arrange(transcript_id) %>% 
  fill(gene_biotype, transcript_biotype) %>% 
  relocate(ORF_id,
           .before = orfanage_status) %>% 
  mutate(across(.cols = c(NMD_status,
                          stop_to_lastEJ,
                          num_of_downEJs,
                          UTR3_length,
  ),
  .fns = ~case_when(type %in% c("exon", "CDS") ~ NA,
                    TRUE ~ (.))))

# Export collapsed ORFanage annotation -> NMDRHT version 1.1
rtracklayer::export(NMDRHT_ORFanage_all.sort.gtf_fixed_NMD,
                    "Resources/NMDRHT/NMDRHT.v1.1.gtf")

# Run sort and index
command_sort <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                        "sort ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.1.gtf ",
                        "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.1.sort.gtf"), sep="", collapse = "")
output_sort <- system(command_sort, intern = TRUE)

command_index <- paste(c("/home/volker/Tools/IGV_2.14.1/igvtools ",
                         "index ",
                         "/home/volker/2025_UPF1_NMDRHT/Resources/NMDRHT/NMDRHT.v1.1.sort.gtf"), sep="", collapse = "")
output_index <- system(command_index, intern = TRUE)

NMDRHT_ORFanage_all.sort.gtf_fixed_NMD %>% 
  filter(type=="transcript") %>% 
  dplyr::count(NMD_status, ORF_type)

# FactR for GENCODE canonical ---------------------------------------------

# Load canonical transcript only GENCODE GTF file
gencode.v42.annotation.canonical.gtf_for_factR <- factR::importGTF("Resources/GENCODE/gencode.v42.annotation.canonical.gtf") 

# Predict NMD sensitivity 
gencode.v42.annotation.canonical_NMDprediction.out <- predictNMD(gencode.v42.annotation.canonical.gtf_for_factR)  %>% 
  dplyr::rename("transcript_id" = "transcript") %>% 
  left_join(gtf_gencode_df_short %>% filter(type == "transcript"))

# Stats
gencode.v42.annotation.canonical_NMDprediction.out %>% 
  dplyr::count(is_NMD)

# Export to csv - used by ORFquant later on
gencode.v42.annotation.canonical_NMDprediction.out %>% write_csv("Resources/GENCODE/gencode.v42.annotation.canonical_NMDprediction.csv")

##
#### Checkpoint #4 - save data --------------------------------------------------------------
##

save(NMDRHT_ORFanage_all.sort.gtf_fixed,
     NMDRHT_ORFanage_NMDprediction.out_annotated,
     file = paste0("Resources/NMDRHT/",Sys.Date(),"_NMDRHT_checkpoint_4_datasources.rds"))

# Load if necessary
# load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_4_datasources.rds")

##
# Export simplified NMDRHT tbl ---------------------------------------------
##

NMDRHT.v1.1 <- rtracklayer::readGFF("Resources/NMDRHT/NMDRHT.v1.1.sort.gtf",
                                              filter = list(type=c("transcript")))

NMDRHT.v1.1_tbl <- dplyr::as_tibble(NMDRHT.v1.1) %>% 
  dplyr::select(gene_id, transcript_id,
                gene_name, transcript_name,
                gene_biotype, transcript_biotype,
                GENCODE_gene_id, GENCODE_transcript_id,
                GENCODE_gene_type,
                structural_category, UIC_total_support, ORF_type,
                NMD_status, stop_to_lastEJ, num_of_downEJs, UTR3_length)

NMDRHT.v1.1_tbl %>% write_csv("Resources/NMDRHT/NMDRHT.v1.1_short.csv")

NMDRHT.v1.1_tbl <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_short.csv")


# Stats: How many tx per gene ----------------------------------------------------

NMDRHT.v1.1_tbl %>% 
  distinct(gene_id) %>% 
  dplyr::count()

NMDRHT.v1.1_tbl %>% 
  distinct(transcript_id) %>% 
  dplyr::count()

# 45757 tx per 22014 genes = 2.08 tx/gene

gtf_gencode_df_short %>% 
  filter(type == "gene") %>% 
  distinct(gene_id) %>% 
  dplyr::count()

gtf_gencode_df_short %>% 
  filter(type == "transcript") %>% 
  distinct(transcript_id) %>% 
  dplyr::count()

# 252242 tx per 62649 genes = 4.03 tx/gene