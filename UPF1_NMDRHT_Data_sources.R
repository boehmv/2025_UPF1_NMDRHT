#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Data_sources
# Objective: Import general data sources (DGE, DTE, DTU, AS) and aggregate relevant data for re-use
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)

##
# Load data sources --------------------------------------------
##

# Load experimentSet overview 

UPF1_NMDRHT_datasets_complete <- read_csv("Resources/Metadata/UPF1_NMDRHT_datasets_complete.csv", trim_ws = TRUE) %>% 
  mutate(experiment_path = case_when(seq_type == "Ribo-Seq" ~ paste0(path,"/experiment_forCRSA.txt"),
                                     TRUE ~ paste0(path,"/experiment.txt"))) %>% 
  mutate(DESeq2_DGE_path = paste0(path,"/DESeq2/DGE/DESeq2_DGE_cat.csv")) %>% 
  mutate(Swish_DGE_path = paste0(path,"/Swish/DGE/Swish_DGE_cat.csv")) %>% 
  mutate(edgeR_DTE_path = paste0(path,"/edgeR/DTE/edgeR_DTE_cat.csv")) %>% 
  mutate(edgeR_DTE_NMDRHT_path = paste0(path,"/edgeR/DTE_NMDRegHumanTxome/edgeR_DTE_cat_NMDRegHumanTxome.csv")) %>% 
  mutate(ISAR_path = paste0(path,"/ISAR/SwitchList_filt_Analyzed.csv")) %>%
  mutate(ISAR_NMDRHT_path = paste0(path,"/ISAR_NMD/SwitchList_filt_Analyzed.csv")) %>%
  mutate(LeafCutter_path = paste0(path,"/leafcutter/leafcutter_AS_cat.csv")) %>% 
  mutate(experimentSet = fct_inorder(experimentSet)) %>% 
  mutate(publicationName = fct_inorder(publicationName))

# Export datasets for GitHub
UPF1_NMDRHT_datasets_complete %>% 
  dplyr::select(-c(path, experiment_path, DESeq2_DGE_path, Swish_DGE_path, edgeR_DTE_path, edgeR_DTE_NMDRegHumanTxome_path, ISAR_path, ISAR_NMDRegHumanTxome_path, LeafCutter_path)) %>% 
  write_csv("Resources/Metadata/UPF1_NMDRHT_datasets_forGitHub.csv")

# Load sample and condition mapping
UPF1_NMDRHT_datasets_experiments_complete <- readr::read_tsv(UPF1_NMDRHT_datasets_complete$experiment_path, 
                                                             id = "experiment_path",
                                                             col_names = c("sample", "condition_2")) %>% 
  left_join(UPF1_NMDRHT_datasets_complete %>% dplyr::select(-c(meta_repository, meta_ID, fastq_repository, fastq_ID, DESeq2_DGE_path, Swish_DGE_path, edgeR_DTE_path, edgeR_DTE_NMDRegHumanTxome_path, ISAR_path, ISAR_NMDRegHumanTxome_path, LeafCutter_path))) %>% 
  arrange((experimentSet)) %>% 
  mutate(condition_2 = fct_inorder(condition_2))

# Export datasets for GitHub
UPF1_NMDRHT_datasets_experiments_complete %>% 
  dplyr::select(-c(path, experiment_path)) %>% 
  write_csv("Resources/Metadata/UPF1_NMDRHT_datasets_experiments_forGitHub.csv")

# Filter for SR_RNA-Seq
UPF1_NMDRHT_datasets <- UPF1_NMDRHT_datasets_complete %>% 
  filter(seq_type == "SR_RNA-Seq") %>% 
  dplyr::select(-seq_type)

UPF1_NMDRHT_datasets %>% 
  write_csv("Resources/Metadata/UPF1_NMDRHT_datasets.csv")

UPF1_NMDRHT_datasets_experiments <- UPF1_NMDRHT_datasets_experiments_complete %>% 
  filter(seq_type == "SR_RNA-Seq") %>% 
  dplyr::select(-seq_type)

UPF1_NMDRHT_datasets_experiments %>% 
  write_csv("Resources/Metadata/UPF1_NMDRHT_datasets_experiment.csv")

# Check if all data files exist
UPF1_NMDRHT_datasets_source_check <- UPF1_NMDRHT_datasets %>% 
  dplyr::select(experimentSet, path, publicationName) %>% 
  mutate(DESeq2_DGE_check = file.exists(UPF1_NMDRHT_datasets %>% pull(DESeq2_DGE_path))) %>% 
  mutate(Swish_DGE_check = file.exists(UPF1_NMDRHT_datasets %>% pull(Swish_DGE_path))) %>% 
  mutate(edgeR_DTE_check = file.exists(UPF1_NMDRHT_datasets %>% pull(edgeR_DTE_path))) %>% 
  mutate(edgeR_DTE_NMDRHT_check = file.exists(UPF1_NMDRHT_datasets %>% pull(edgeR_DTE_NMDRHT_path))) %>% 
  mutate(ISAR_DTU_check = file.exists(UPF1_NMDRHT_datasets %>% pull(ISAR_path))) %>% 
  mutate(ISAR_DTU_NMDRHT_check = file.exists(UPF1_NMDRHT_datasets %>% pull(ISAR_NMDRHT_path))) %>% 
  mutate(LeafCutter_AS_check = file.exists(UPF1_NMDRHT_datasets %>% pull(LeafCutter_path)))

##
# Load DESeq2 DGE data -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(DESeq2_DGE_path))))

DESeq2_DGE_combined <- readr::read_csv(UPF1_NMDRHT_datasets %>% 
                                         pull(DESeq2_DGE_path),
                                       id = "DESeq2_DGE_path",
                                       col_select = (-1)) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other",
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(DESeq2_DGE_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-DESeq2_DGE_path)

# Save pre-parsed datasources - as csv
DESeq2_DGE_combined  %>%  write_csv(file.path("Resources", "DESeq2_DGE_combined.csv"))

##
# Load Swish DGE data -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(Swish_DGE_path))))

Swish_DGE_combined <- UPF1_NMDRHT_datasets %>% 
  pull(Swish_DGE_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'Swish_DGE_path') %>% 
  # filter(keep == TRUE) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other",
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(Swish_DGE_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-Swish_DGE_path)

# Save pre-parsed datasources - as csv
Swish_DGE_combined  %>%  write_csv(file.path("Resources", "Swish_DGE_combined.csv"))

##
# Load edgeR DTE data -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(edgeR_DTE_path))))

edgeR_DTE_combined <- UPF1_NMDRHT_datasets %>% 
  pull(edgeR_DTE_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'edgeR_DTE_path') %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          transcript_type == "protein_coding" ~ "coding",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(edgeR_DTE_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-edgeR_DTE_path)

# Save pre-parsed datasources - as csv
edgeR_DTE_combined  %>%  write_csv(file.path("Resources", "edgeR_DTE_combined.csv"))

##
# Load edgeR DTE NMDRHT files   -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(edgeR_DTE_NMDRHT_path))))

edgeR_DTE_NMDRHT_combined <- UPF1_NMDRHT_datasets %>% 
  pull(edgeR_DTE_NMDRHT_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'edgeR_DTE_NMDRHT_path') %>% 
  mutate(type = case_when(!transcript_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          transcript_biotype == "protein_coding" ~ "coding",
                          transcript_biotype == "nonsense_mediated_decay" ~ "NMD",
                          transcript_biotype == "lncRNA" ~ "lncRNA")) %>% 
  # Remove v1.1 (ORFanage-mostly)-derived parameters
  dplyr::select(-c(ORF_type, NMD_status, stop_to_lastEJ, num_of_downEJs, UTR3_length)) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(edgeR_DTE_NMDRHT_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-edgeR_DTE_NMDRHT_path)

# Save pre-parsed datasources - as csv
edgeR_DTE_NMDRHT_combined  %>%  write_csv(file.path("Resources", "edgeR_DTE_NMDRHT_combined.csv"))

##
# Load ISAR DTU data -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(ISAR_path))))

ISAR_DTU_combined <- UPF1_NMDRHT_datasets %>% 
  pull(ISAR_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'ISAR_path')  %>% 
  mutate(type = case_when(!iso_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          iso_biotype == "protein_coding" ~ "coding",
                          iso_biotype == "nonsense_mediated_decay" ~ "NMD",
                          iso_biotype == "lncRNA" ~ "lncRNA")) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(ISAR_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-ISAR_path)

# Save pre-parsed datasources - as csv
ISAR_DTU_combined  %>%  write_csv(file.path("Resources", "ISAR_DTU_combined.csv"))

##
# Load ISAR DTU NMDRHT data   -----------------------------------------------------------
##

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% pull(ISAR_NMDRHT_path))))

ISAR_DTU_NMDRHT_combined <- UPF1_NMDRHT_datasets %>% 
  pull(ISAR_NMDRHT_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'ISAR_NMDRHT_path')  %>% 
  mutate(type = case_when(!iso_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          iso_biotype == "protein_coding" ~ "coding",
                          iso_biotype == "nonsense_mediated_decay" ~ "NMD",
                          iso_biotype == "lncRNA" ~ "lncRNA")) %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(ISAR_NMDRHT_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-ISAR_NMDRHT_path)

# Save pre-parsed datasources - as csv
ISAR_DTU_NMDRHT_combined  %>%  write_csv(file.path("Resources", "ISAR_DTU_NMDRHT_combined.csv"))

##
# Load LeafCutter AS data -----------------------------------------------------------
##

# Comment: Datasets HeLa_UPF1_KD_Longman_2020 (total) and HCT116_UPF1_AID_TestSeq_this_Study were excluded

# Check if all files are present - otherwise stop here
stopifnot("*** Not all files are present! ***" = all(file.exists(UPF1_NMDRHT_datasets %>% 
                                                                   filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                                                                                                "HCT116_UPF1_AID_TestSeq_this_Study"))
                                                                 %>% pull(LeafCutter_path))))

LeafCutter_AS_combined <- UPF1_NMDRHT_datasets %>% 
  filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                               "HCT116_UPF1_AID_TestSeq_this_Study")) %>% 
  pull(LeafCutter_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'LeafCutter_path')  %>% 
  left_join(UPF1_NMDRHT_datasets %>% dplyr::select(LeafCutter_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-LeafCutter_path)

# Save pre-parsed datasources - as csv
LeafCutter_AS_combined  %>%  write_csv(file.path("Resources", "LeafCutter_AS_combined.csv"))

# intron-annotated
LeafCutter_AS_introns_annotated <- UPF1_NMDRHT_datasets %>% 
  filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                               "HCT116_UPF1_AID_TestSeq_this_Study")) %>% 
  mutate(Intron_path = paste0(path,"/leafcutter/leafcutter_introns_cat.csv")) %>% 
  pull(Intron_path) %>% 
  map_dfr(read_csv, id = 'Intron_path')  %>% 
  left_join(UPF1_NMDRHT_datasets %>%
              mutate(Intron_path = paste0(path,"/leafcutter/leafcutter_introns_cat.csv")) %>% 
              dplyr::select(Intron_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-Intron_path)

# Save pre-parsed datasources - as csv
LeafCutter_AS_introns_annotated  %>%  write_csv(file.path("Resources", "LeafCutter_AS_introns_annotated.csv"))

# intron-annotated_summary
LeafCutter_AS_introns_annotated_summary <- UPF1_NMDRHT_datasets %>% 
  filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                               "HCT116_UPF1_AID_TestSeq_this_Study")) %>% 
  mutate(Intron_summary_path = paste0(path,"/leafcutter/leafcutter_intronSummary_cat.csv")) %>% 
  pull(Intron_summary_path) %>% 
  map_dfr(read_csv, id = 'Intron_summary_path')  %>% 
  left_join(UPF1_NMDRHT_datasets %>%
              mutate(Intron_summary_path = paste0(path,"/leafcutter/leafcutter_intronSummary_cat.csv")) %>% 
              dplyr::select(Intron_summary_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-Intron_summary_path)

# Save pre-parsed datasources - as csv
LeafCutter_AS_introns_annotated_summary  %>%  write_csv(file.path("Resources", "LeafCutter_AS_introns_annotated_summary.csv"))

# cluster-annotated
LeafCutter_AS_clusters_annotated <- UPF1_NMDRHT_datasets %>% 
  filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                               "HCT116_UPF1_AID_TestSeq_this_Study")) %>% 
  mutate(Cluster_path = paste0(path,"/leafcutter/leafcutter_clusters_cat.csv")) %>% 
  pull(Cluster_path) %>% 
  map_dfr(read_csv, id = 'Cluster_path')  %>% 
  left_join(UPF1_NMDRHT_datasets %>%
              mutate(Cluster_path = paste0(path,"/leafcutter/leafcutter_clusters_cat.csv")) %>% 
              dplyr::select(Cluster_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-Cluster_path)

# Save pre-parsed datasources - as csv
LeafCutter_AS_clusters_annotated  %>%  write_csv(file.path("Resources", "LeafCutter_AS_clusters_annotated.csv"))

# intron-annotated_summary
LeafCutter_AS_clusters_annotated_summary <- UPF1_NMDRHT_datasets %>% 
  filter(!experimentSet %in% c("HeLa_UPF1_KD_Longman_2020",
                               "HCT116_UPF1_AID_TestSeq_this_Study")) %>% 
  mutate(Cluster_summary_path = paste0(path,"/leafcutter/leafcutter_clusterSummary_cat.csv")) %>% 
  pull(Cluster_summary_path) %>% 
  map_dfr(read_csv, id = 'Cluster_summary_path')  %>% 
  left_join(UPF1_NMDRHT_datasets %>%
              mutate(Cluster_summary_path = paste0(path,"/leafcutter/leafcutter_clusterSummary_cat.csv")) %>% 
              dplyr::select(Cluster_summary_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-Cluster_summary_path)

# Save pre-parsed datasources - as csv
LeafCutter_AS_clusters_annotated_summary  %>%  write_csv(file.path("Resources", "LeafCutter_AS_clusters_annotated_summary.csv"))
