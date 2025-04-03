#!/usr/bin/env Rscript

# Title: UPF1_Revision_Figure7
# Objective: Code for generating Panels of Figure 7 + Figure S7 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_Revision_Analysis"

###
## Rev_1 - Figure 7 - (A) Alternative Splicing -----------------------------------------------------------
###

# Load Main transcript NMDRHT table
NMDRHT.v1.2_MainTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv")

# Load LeafCutter AS cluster summary
LeafCutter_AS_clusters_annotated_summary  <-  read_csv(file.path("Resources", "LeafCutter_AS_clusters_annotated_summary.csv"))

### Cluster summary --------------------------------------------------------
LeafCutter_AS_clusters_annotated_summary %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  filter(!Results == "Number of differentially spliced clusters at FDR = 1e-04") %>% 
  mutate(Results = fct_relevel(Results,
                               "Fully annotated",
                               "Contain unannotated junctions")) %>% 
  ggplot(aes(y=condition_2,
             x=n,
             fill=Results)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_col(color="black",
           position="stack",
           linewidth = 0.2) +
  # geom_text(aes(label = n),
  #           color="black",
  #           position = position_dodge(width = 0.9),
  #           size = 4*0.36,
  #           hjust= -0.25,
  #           # angle = 90
  # ) +
  scale_x_continuous(expand = expansion(c(0.025, 0.25))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("Fully annotated" = "#B80C66",
                             "Contain unannotated junctions" = "#FCA50A")) +
  labs(y="",
       x="Number of cluster",
       fill="Cluster") +
  force_panelsizes(rows = unit(length(LeafCutter_AS_clusters_annotated_summary %>% 
                                        filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                                                    "HCT116_UPF1_FKBP_degradation_this_Study",
                                                                    "HEK293_UPF1_FKBP_degradation_this_Study",
                                                                    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
                                        distinct(condition_2) %>% 
                                        pull(condition_2))*2+0.6, "mm"),
                   cols = unit(20, "mm")) 

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_A_LeafCutter_AS_clusters_annotated_summary_all.pdf"),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###
## Rev_1 - Figure 7 - (B) AS-NMD proteins -----------------------------------------------------------
###

# Load required data
load("Resources/NMDRHT/NMDRHT_cds_seqs_AA_df_GENCODE_match.rds")

### Plot rel. length  -----------------------------------------------------------
NMDRHT_cds_seqs_AA_df_GENCODE_match %>% 
  ggplot(aes(y=rel_proteinWidth_AS_NMD,
             x="AS-NMD")) +
  theme(legend.position="right",
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.major.x = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.y = element_blank(),
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
  geom_violin(alpha=0.25,
              fill="#E5BF85",
              linewidth=0.1,
              position=position_dodge(0.7), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               alpha=0.8,
               fill="#E5BF85",
               position=position_dodge(0.7),
               width=0.2,
               linewidth=0.1,
               color="black") +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(7.5, "mm")) +
  ylim(0,100) +
  guides(fill="none") +
  labs(y="relative protein length (%)\n(NMDRHT/GENCODE-canonical)",
       x="AS-NMD/ntranscripts",
       fill="")

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_B1_AS_NMD_relativeProteinLength.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Plot PFAM_MobiDBLite --------------------------------------------------------------------

# Load data
NMDRHT_interPro_Pfam_MobiDBLite_combined <- read_csv("Resources/NMDRHT/InterProScan/NMDRHT_interPro_Pfam_MobiDBLite_combined.csv")

NMDRHT_interPro_Pfam_MobiDBLite_combined_class <- NMDRHT_interPro_Pfam_MobiDBLite_combined %>% 
  mutate(PFAM_diff_class = case_when(PFAM_diff >= 3 ~ ">=3",
                                     PFAM_diff <= -3 ~ "<=-3",
                                     TRUE ~ as.character(PFAM_diff))) %>% 
  # mutate(PFAM_diff_class = fct_expand(PFAM_diff_class, "3")) %>% 
  mutate(PFAM_diff_class = fct_relevel(PFAM_diff_class,
                                       "<=-3",
                                       "-2",
                                       "-1",
                                       "0",
                                       "1",
                                       "2",
                                       ">=3")) %>% 
  mutate(unstr_diff_class = case_when(unstr_diff >= 3 ~ ">=3",
                                      unstr_diff <= -3 ~ "<=-3",
                                      TRUE ~ as.character(unstr_diff))) %>% 
  mutate(unstr_diff_class = fct_relevel(unstr_diff_class,
                                        "<=-3",
                                        "-2",
                                        "-1",
                                        "0",
                                        "1",
                                        "2",
                                        ">=3")) 

# Plot
NMDRHT_interPro_Pfam_MobiDBLite_combined_class %>% 
  dplyr::count(PFAM_diff_class, unstr_diff_class) %>% 
  ggplot(aes(x=PFAM_diff_class,
             y=unstr_diff_class,
             fill=n,
             size=n)) +
  theme(legend.position="right",
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_line(colour = 'gray', linewidth = 0.1),
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
  geom_point(shape=21) +
  scale_size(range = c(0.5, 2)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_viridis_c(option="mako",
                       begin = 0.1,
                       end=0.9) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm")) +
  labs(x="PFAM domains difference (n)",
       y="Intrinsic disorder difference (n)")

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_B2_AS_NMD_PFAM_ID_region.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 7 - (C) Proteomics - Peptide-level -----------------------------------------------------------
###

### detected-vs-predicted plots ---------------------------------------------
load("Resources/Proteomics/combined_tbl_peptidome_class_unique_detectability.rds")
WCP_peptide_data_20250127_SN_final_log2_combined <- read_csv("Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250127_SN_final_log2_combined.csv")

##### overall numbers ---------------------------------------------
combined_tbl_peptidome_class_unique_detectability %>% 
  dplyr::count(peptide_class) %>% 
  mutate(n_per = n/sum(n)*100) %>% 
  mutate(set = "predicted") %>% 
  bind_rows(WCP_peptide_data_20250127_SN_final_log2_combined %>% 
              dplyr::count(peptide_class) %>% 
              mutate(peptide_class = case_when(is.na(peptide_class) ~ "NA",
                                               TRUE ~ peptide_class)) %>% 
              mutate(n_per = n/sum(n)*100) %>% 
              mutate(set = "detected")) %>% 
  mutate(set = (fct_relevel(set,
                            "detected",
                            "predicted"))) %>% 
  mutate(peptide_class = fct_rev(fct_relevel(peptide_class,
                                             "coding",
                                             "shared",
                                             "NMD",
                                             "NA"))) %>% 
  ggplot(aes(y=set,
             x=n_per,
             fill=peptide_class)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_col(color="black",
           position = position_dodge2(width = 0.9, preserve = "single"),
           linewidth = 0.1) +
  geom_text(aes(label = case_when(peptide_class == "coding" ~ paste0(n,
                                                                     " (",
                                                                     round(n_per,1),
                                                                     "%)"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_dodge2(width = 0.9, 
                                       preserve = "single")
  ) +
  geom_text(aes(label = case_when(peptide_class != "coding" ~ paste0(n,
                                                                     " (",
                                                                     round(n_per,1),
                                                                     "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_dodge2(width = 0.9, 
                                       preserve = "single")
  ) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "shared" = "#85ABE5",
                             "NMD" = "#E5BF86",
                             "NA" = "darkgray")) +
  coord_cartesian(clip = "off") +
  labs(y="",
       x="peptides (%)",
       fill="Peptide class") +
  guides(fill=guide_legend(reverse=T)) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C1_Proteomics_NMD_predicted_SN_detected_peptide_class.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### junction-spanning ---------------------------------------------
combined_tbl_peptidome_class_unique_detectability %>% 
  dplyr::count(peptide_class,peptide_overJunction) %>% 
  group_by(peptide_class) %>% 
  mutate(n_per = n/sum(n)*100) %>% 
  ungroup() %>% 
  mutate(set = "predicted",
         peptide_overJunction = as.character(peptide_overJunction)) %>% 
  bind_rows(WCP_peptide_data_20250127_SN_final_log2_combined %>% 
              filter(!is.na(peptide_overJunction)) %>% 
              dplyr::count(peptide_class,peptide_overJunction) %>% 
              group_by(peptide_class) %>% 
              mutate(n_per = n/sum(n)*100) %>% 
              ungroup() %>% 
              mutate(set = "detected",
                     peptide_overJunction = as.character(peptide_overJunction))) %>% 
  mutate(set = (fct_relevel(set,
                            "detected",
                            "predicted"))) %>% 
  mutate(peptide_overJunction = fct_rev(fct_relevel(peptide_overJunction,
                                                    "FALSE",
                                                    "TRUE"))) %>% 
  mutate(peptide_class = (fct_relevel(peptide_class,
                                      "coding",
                                      "shared",
                                      "NMD"))) %>% 
  ggplot(aes(y=set,
             x=n_per,
             fill=peptide_overJunction)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_col(color="black",
           position = position_stack(),
           linewidth = 0.1) +
  geom_text(aes(label = case_when(peptide_overJunction == "TRUE" & n_per >= 10 ~ paste0(round(n_per,0),
                                                                                        "%"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(peptide_overJunction == "TRUE" & n_per < 10 ~ paste0(round(n_per,0),
                                                                                       "%"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 1)
  ) +
  geom_text(aes(label = case_when(peptide_overJunction != "TRUE" ~ paste0(round(n_per,0),
                                                                          "%"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("TRUE" = "#006666",
                             "FALSE" = "lightgray")) +
  coord_cartesian(clip = "off") +
  labs(y="",
       x="peptides (%)",
       fill="Junction-\nspanning\npeptide") +
  guides(fill=guide_legend(reverse=T)) +
  facet_wrap(~peptide_class,nrow=3) +
  force_panelsizes(rows = unit(7.5, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C2_Proteomics_NMD_predicted_SN_detected_overJunction.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### terminal peptide ---------------------------------------------
combined_tbl_peptidome_class_unique_detectability %>% 
  dplyr::count(peptide_class,terminal_peptide) %>% 
  group_by(peptide_class) %>% 
  mutate(n_per = n/sum(n)*100) %>% 
  ungroup() %>% 
  mutate(set = "predicted",
         terminal_peptide = as.character(terminal_peptide)) %>% 
  bind_rows(WCP_peptide_data_20250127_SN_final_log2_combined %>% 
              filter(!is.na(terminal_peptide)) %>% 
              dplyr::count(peptide_class,terminal_peptide) %>% 
              group_by(peptide_class) %>% 
              mutate(n_per = n/sum(n)*100) %>% 
              ungroup() %>% 
              mutate(set = "detected",
                     terminal_peptide = as.character(terminal_peptide))) %>% 
  mutate(set = (fct_relevel(set,
                            "detected",
                            "predicted"))) %>% 
  mutate(terminal_peptide = fct_rev(fct_relevel(terminal_peptide,
                                                "FALSE",
                                                "TRUE"))) %>% 
  mutate(peptide_class = (fct_relevel(peptide_class,
                                      "coding",
                                      "shared",
                                      "NMD"))) %>% 
  ggplot(aes(y=set,
             x=n_per,
             fill=terminal_peptide)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_col(color="black",
           position = position_stack(),
           linewidth = 0.1) +
  geom_text(aes(label = case_when(terminal_peptide == "TRUE" & n_per >= 10 ~ paste0(round(n_per,0),
                                                                                    "%"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(terminal_peptide == "TRUE" & n_per < 10 ~ paste0(round(n_per,0),
                                                                                   "%"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 1)
  ) +
  geom_text(aes(label = case_when(terminal_peptide != "TRUE" ~ paste0(round(n_per,0),
                                                                      "%"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("TRUE" = "#006666",
                             "FALSE" = "lightgray")) +
  coord_cartesian(clip = "off") +
  labs(y="",
       x="peptides (%)",
       fill="Terminal peptide") +
  guides(fill=guide_legend(reverse=T)) +
  facet_wrap(~peptide_class,nrow=3) +
  force_panelsizes(rows = unit(7.5, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C3_Proteomics_NMD_predicted_SN_detected_terminalPeptide.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### detectability ---------------------------------------------

combined_tbl_peptidome_class_unique_detectability %>% 
  dplyr::select(sequence, peptide_class, RF_score) %>% 
  mutate(set = "predicted") %>% 
  bind_rows(WCP_peptide_data_20250127_SN_final_log2_combined %>% 
              filter(!is.na(peptide_class)) %>% 
              dplyr::select(Sequence, peptide_class, RF_score) %>% 
              dplyr::rename("sequence" = "Sequence") %>% 
              mutate(set = "detected")) %>% 
  mutate(peptide_class = fct_rev(fct_relevel(peptide_class,
                                             "coding",
                                             "shared",
                                             "NMD"))) %>% 
  # group_by(peptide_class, set) %>% 
  # get_summary_stats(RF_score)
  ggplot(aes(x=RF_score,
             y=set,
             fill=peptide_class)) +
  theme(legend.position="right",
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.major.x = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.y = element_blank(),
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
  geom_violin(alpha=0.25,
              linewidth=0.1,
              scale = "width",
              width=0.8,
              position=position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               alpha=0.8,
               position=position_dodge(0.9),
               width=0.2,
               linewidth=0.1,
               color="black") +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "shared" = "#85ABE5",
                             "NMD" = "#E5BF86")) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm")) +
  xlim(0,1) +
  guides(fill=guide_legend(reverse=T)) +
  labs(x="peptide detectability\n(ProteomicsDB-based random forest score)",
       y="peptides",
       fill="peptide class")

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C4_NMD_peptideDetectability_RFscore.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Overlap SN & DIA-NN distinct --------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined <- read_csv("Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250127_SN_final_log2_combined.csv")
WCP_peptide_data_20250206_DiaNN_final_log2_combined <- read_csv("Resources/Proteomics/2025_02_06_Reanalysis/WCP_peptide_data_20250206_DiaNN_final_log2_combined.csv")


#### Overlap all  ----------------------------------------------------------------

##### nVenn -------------------------------------------------------------

# Select ids
WCP_peptide_data_20250127_SN_all_ids <- WCP_peptide_data_20250127_SN_final_log2_combined %>%
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  pull(gene_seq)

WCP_peptide_data_20250206_diaNN_all_ids <- WCP_peptide_data_20250206_DiaNN_final_log2_combined %>%
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  pull(gene_seq)

WCP_peptide_data_20250206_overlap_all_myNV <- nVennR::plotVenn(list(WCP_peptide_data_20250127_SN_all_ids,
                                                                    WCP_peptide_data_20250206_diaNN_all_ids), 
                                                               sNames=c("SN", 
                                                                        "DiaNN"), showPlot = FALSE)

showSVG(WCP_peptide_data_20250206_overlap_all_myNV,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#01665E", "#B2B3B3"))

WCP_peptide_data_20250206_SN_all_only_IDs <- getVennRegion(WCP_peptide_data_20250206_overlap_all_myNV, c(1,0))
WCP_peptide_data_20250206_diaNN_all_only_IDs <-getVennRegion(WCP_peptide_data_20250206_overlap_all_myNV, c(0,1))
WCP_peptide_data_20250206_SN_diaNN_overlap_all_IDs <- getVennRegion(WCP_peptide_data_20250206_overlap_all_myNV, c(1,1))

##### Correlation -------------------------------------------------------------
WCP_peptide_data_20250127_SN_final_log2_combined_expand <- read_csv("Resources/Proteomics/2025_01_27_Reanalysis/WCP_peptide_data_20250127_SN_final_log2_combined_expand.csv")

WCP_peptide_data_20250127_SN_final_log2_combined_NMD <- WCP_peptide_data_20250127_SN_final_log2_combined %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250127_SN_final_log2_combined_expand_NMD <- WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand <- read_csv("Resources/Proteomics/2025_02_06_Reanalysis/WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand.csv")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD <- WCP_peptide_data_20250206_DiaNN_final_log2_combined %>% 
  filter(peptide_class == "NMD")

WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand_NMD <- WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand %>% 
  filter(peptide_class == "NMD")

SN_DiaNN_forCor <- WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  mutate(tool = "Spectronaut") %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_all_IDs) %>% 
  dplyr::select(gene_seq, condition, sample, Intensity) %>% 
  dplyr::rename("SN_intensity" = "Intensity") %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand %>% 
              mutate(tool = "DiaNN") %>% 
              mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
              filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_all_IDs) %>% 
              dplyr::select(gene_seq, condition, sample, Intensity) %>% 
              dplyr::rename("DiaNN_intensity" = "Intensity"))

SN_DiaNN_cor <- SN_DiaNN_forCor %>% 
  group_by(sample) %>% 
  rstatix::cor_test(SN_intensity, DiaNN_intensity) 

#### Overlap NMD  ----------------------------------------------------------------

##### nVenn -------------------------------------------------------------

# Select ids
WCP_peptide_data_20250127_SN_ids <- WCP_peptide_data_20250127_SN_final_log2_combined_NMD %>%
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  pull(gene_seq)

WCP_peptide_data_20250206_diaNN_ids <- WCP_peptide_data_20250206_DiaNN_final_log2_combined_NMD %>%
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  pull(gene_seq)

WCP_peptide_data_20250206_overlap_myNV <- nVennR::plotVenn(list(WCP_peptide_data_20250127_SN_ids,
                                                                WCP_peptide_data_20250206_diaNN_ids), 
                                                           sNames=c("SN", 
                                                                    "DiaNN"), showPlot = FALSE)

showSVG(WCP_peptide_data_20250206_overlap_myNV,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#01665E", "#B2B3B3"))

WCP_peptide_data_20250206_SN_only_IDs <- getVennRegion(WCP_peptide_data_20250206_overlap_myNV, c(1,0))
WCP_peptide_data_20250206_diaNN_only_IDs <-getVennRegion(WCP_peptide_data_20250206_overlap_myNV, c(0,1))
WCP_peptide_data_20250206_SN_diaNN_overlap_IDs <- getVennRegion(WCP_peptide_data_20250206_overlap_myNV, c(1,1))

# Check for which gene multiple HQ NMD peptides are found
WCP_peptide_data_20250127_SN_final_log2_combined_NMD %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  dplyr::count(gene_name) %>% 
  arrange(desc(n))

##### Correlation -------------------------------------------------------------

SN_DiaNN_forCor_NMD <- WCP_peptide_data_20250127_SN_final_log2_combined_expand_NMD %>% 
  mutate(tool = "Spectronaut") %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  dplyr::select(gene_seq, condition, sample, Intensity) %>% 
  dplyr::rename("SN_intensity" = "Intensity") %>% 
  left_join(WCP_peptide_data_20250206_DiaNN_final_log2_combined_expand_NMD  %>% 
              mutate(tool = "DiaNN") %>% 
              mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
              filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
              dplyr::select(gene_seq, condition, sample, Intensity) %>% 
              dplyr::rename("DiaNN_intensity" = "Intensity"))

SN_DiaNN_cor_NMD <- SN_DiaNN_forCor_NMD %>% 
  group_by(sample) %>% 
  rstatix::cor_test(SN_intensity, DiaNN_intensity) 

###### Combined Plot -----------------------------------------------------------

SN_DiaNN_cor %>% 
  mutate(set = "all") %>% 
  bind_rows(SN_DiaNN_cor_NMD %>% 
              mutate(set = "NMD")) %>% 
  mutate(condition = case_when(str_detect(sample, "12h_DMSO") ~ "12h_DMSO",
                               str_detect(sample, "12h_IAA") ~ "12h_IAA",
                               str_detect(sample, "24h_DMSO") ~ "24h_DMSO",
                               str_detect(sample, "24h_IAA") ~ "24h_IAA",
                               str_detect(sample, "06h_DMSO") ~ "06h_DMSO",
                               str_detect(sample, "06h_IAA") ~ "06h_IAA",
                               str_detect(sample, "15h_DMSO") ~ "15h_DMSO",
                               str_detect(sample, "15h_IAA") ~ "15h_IAA",
                               str_detect(sample, "18h_DMSO") ~ "18h_DMSO",
                               str_detect(sample, "18h_IAA") ~ "18h_IAA"))   %>%
  mutate(condition = fct_rev(condition)) %>% 
  group_by(condition) %>% 
  mutate(cor_mean = mean(cor)) %>% 
  ungroup() %>% 
  ggplot(aes(x=cor,
             y=condition,
             color=set)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_segment(aes(x=0, xend=cor_mean),
               linewidth=0.2,
               position = position_dodge(width=0.75),
               alpha=0.75, show.legend=FALSE) +
  geom_point(aes(fill=set),
             color="black",
             stroke=0.1,
             size=0.75,
             shape=21,
             alpha=0.75,
             position = position_dodge(width=0.75)) +
  scale_color_manual(values=c("all" = "#66B2B2",
                              "NMD" = "#B09771")) +
  scale_fill_manual(values=c("all" = "#66B2B2",
                             "NMD" = "#B09771")) +
  xlim(c(0,1)) +
  labs(x="pearson correlation coefficient",
       y="",
       fill="Peptides") +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C7_SN_DiaNN_pearson.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 



### SN - plot raw intensities of HQ peptides --------------------------------

#### With NMD reason ---------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  mutate(condition = (fct_relevel(condition,
                                  "06h_DMSO",
                                  "12h_DMSO",
                                  "15h_DMSO",
                                  "18h_DMSO",
                                  "24h_DMSO",
                                  "06h_IAA",
                                  "12h_IAA",
                                  "15h_IAA",
                                  "18h_IAA",
                                  "24h_IAA"))) %>% 
  arrange(condition) %>% 
  mutate(sample = fct_inorder(sample)) %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  # filter(!is.na(intensity)) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(gene_seq) %>% 
  mutate(mean_log2_intensity = mean(log2_Intensity, na.rm=TRUE),
         valid_values = sum(!is.na(log2_Intensity))) %>% 
  # 35% valid values (allows DMSO and 06h IAA to be empty)
  mutate(selected = case_when(valid_values >= 0.35*43 ~ TRUE,
                              TRUE ~ FALSE)) %>% 
  ungroup() %>% 
  dplyr::arrange(NMD_peptide_reason_simple, selected, (mean_intensity)) %>% 
  mutate(gene_seq = fct_inorder(gene_seq)) %>% 
  # mutate(sequence = fct_rev(fct_inorder(sequence))) %>% 
  # mutate(sample = fct_rev(sample)) %>% 
  ggplot(aes(y=gene_seq,
             x=sample,
             fill=log2_Intensity)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4)) +
  ggside::geom_ysidetile(aes(x = "NMD reason", yfill = NMD_peptide_reason))  +
  ggside::scale_yfill_manual(values=c("none" = "#214D65",
                                      "lncRNA" = "#76716E",
                                      "AS_NMD"="#E5BF85",
                                      "AS_NMD_UTR3" = "#C1AC65",
                                      "novel" = "#9A9A49",
                                      "uORF" = "#708831",
                                      "overl_uORF" = "#43761E",
                                      "multiple" = "gray30",
                                      "other" = "#D9DADA")) +
  # geom_ysidetile(aes(x = "40% valid values", yfill = selected))  +
  # scale_yfill_manual(values=c("TRUE" = "#963B5A",
  #                             "FALSE" = "gray90")) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C8_Proteomics_SN_NMD_HQ_peptides_log2_intensity_NMD_Reason.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### With Valid Values ---------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_expand %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  mutate(condition = (fct_relevel(condition,
                                  "06h_DMSO",
                                  "12h_DMSO",
                                  "15h_DMSO",
                                  "18h_DMSO",
                                  "24h_DMSO",
                                  "06h_IAA",
                                  "12h_IAA",
                                  "15h_IAA",
                                  "18h_IAA",
                                  "24h_IAA"))) %>% 
  arrange(condition) %>% 
  mutate(sample = fct_inorder(sample)) %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  # filter(!is.na(intensity)) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(gene_seq) %>% 
  mutate(mean_log2_intensity = mean(log2_Intensity, na.rm=TRUE),
         valid_values = sum(!is.na(log2_Intensity))) %>% 
  # 35% valid values (allows DMSO and 06h IAA to be empty)
  mutate(selected = case_when(valid_values >= 0.35*43 ~ TRUE,
                              TRUE ~ FALSE)) %>% 
  ungroup() %>% 
  # distinct(gene_seq, selected) %>% 
  # dplyr::count(selected)
  dplyr::arrange(NMD_peptide_reason_simple, selected, (mean_intensity)) %>% 
  mutate(gene_seq = fct_inorder(gene_seq)) %>% 
  # mutate(sequence = fct_rev(fct_inorder(sequence))) %>% 
  # mutate(sample = fct_rev(sample)) %>% 
  ggplot(aes(y=gene_seq,
             x=sample,
             fill=log2_Intensity)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4)) +
  # geom_ysidetile(aes(x = "NMD reason", yfill = NMD_peptide_reason))  +
  # scale_yfill_manual(values=c("none" = "#214D65",
  #                             "lncRNA" = "#76716E",
  #                             "AS_NMD"="#E5BF85",
  #                             "AS_NMD_UTR3" = "#C1AC65",
  #                             "novel" = "#9A9A49",
  #                             "uORF" = "#708831",
  #                             "overl_uORF" = "#43761E",
  #                             "multiple" = "gray30",
  #                             "other" = "#D9DADA")) +
  ggside::geom_ysidetile(aes(x = "35% valid values", yfill = selected))  +
  ggside::scale_yfill_manual(values=c("TRUE" = "#963B5A",
                                      "FALSE" = "gray90")) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C9_Proteomics_SN_NMD_HQ_peptides_log2_intensity_ValidValues.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### SN - Add Gene and tx-level info ------------------------------------------------------

#### Plots -------------------------------------------------------------------
load("Resources/Proteomics/WCP_peptide_data_20250127_SN_all_GeneTxORF.rds")

##### Prepare for Plot -------------------------------------------------------------------
WCP_peptide_data_20250127_SN_all_GeneTxORF_forPlot <- WCP_peptide_data_20250127_SN_all_GeneTxORF %>% 
  filter(peptide_class == "NMD") %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
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
                               str_detect(sample, "18h_IAA") ~ "18h_IAA")) %>%   
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(NMD_peptide_reason_simple, (mean_intensity)) %>% 
  mutate(gene_seq = fct_inorder(gene_seq)) %>% 
  mutate(sample = fct_rev(sample))

####Expand NMD Intensity ----------------------------------------------------

load("Resources/Proteomics/WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA.rds")

WCP_peptide_data_20250127_SN_final_log2_combined_Imputed_expand <- WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  dplyr::select(-starts_with(c("ANOVA", "modF"))) %>% 
  pivot_longer(cols = c(starts_with("Intensity."),
                        starts_with("log2_Intensity."), 
                        starts_with("Imputed_log2_Intensity.")),
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
##### Raw log2 Intensity ------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_Imputed_expand %>% 
  separate(Gene.names, c("gene_names_plain", NA)) %>% 
  mutate(gene_seq = paste0(gene_names_plain," (",Sequence, ")")) %>% 
  # filter(!is.na(intensity)) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  dplyr::arrange(desc(mean_intensity)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(Sequence = fct_rev(fct_inorder(Sequence))) %>% 
  mutate(sample = (sample)) %>% 
  ggplot(aes(x=sample,
             y=gene_seq,
             fill=log2_Intensity)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(100, "mm"),
                   cols = unit(80, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C10_Proteomics_NMD_peptides_intensity_SN_afterModF.pdf"),
       width = 40,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

##### Imputed log2 Intensity ------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_Imputed_expand %>% 
  separate(Gene.names, c("gene_names_plain", NA)) %>% 
  mutate(gene_seq = paste0(gene_names_plain," (",Sequence, ")")) %>% 
  # filter(!is.na(intensity)) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  dplyr::arrange(desc(mean_intensity)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(Sequence = fct_rev(fct_inorder(Sequence))) %>% 
  mutate(sample = (sample)) %>% 
  ggplot(aes(x=sample,
             y=gene_seq,
             fill=Imputed_log2_Intensity)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(100, "mm"),
                   cols = unit(80, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C11_Proteomics_NMD_peptides_Imputed_intensity_SN_afterModF.pdf"),
       width = 40,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

####Expand modF ----------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_modF_expand <- WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  dplyr::select(-starts_with(c("ANOVA", "log2_Intensity.", "Imputed_log2_Intensity."))) %>% 
  pivot_longer(cols = c(starts_with(c("Intensity.",
                                      "modF_Imputed_log2_Intensity."))),
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

##### Raw log2 Intensity ------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_modF_expand %>% 
  separate(Gene.names, c("gene_names_plain", NA)) %>% 
  mutate(gene_seq = paste0(gene_names_plain," (",Sequence, ")")) %>% 
  # filter(!is.na(intensity)) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_intensity = mean(Intensity, na.rm=TRUE)) %>% 
  ungroup() %>% 
  dplyr::arrange(desc(mean_intensity)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(Sequence = fct_rev(fct_inorder(Sequence))) %>% 
  mutate(sample = (sample)) %>% 
  ggplot(aes(x=condition,
             y=gene_seq,
             fill=modF_Imputed_log2_Intensity)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(100, "mm"),
                   cols = unit(80, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C12_Proteomics_NMD_peptides_ModF_intensity_SN_afterModF.pdf"),
       width = 40,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

####Expand ANOVA ----------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand <- WCP_peptide_data_20250127_SN_final_log2_combined_modF_ANOVA %>% 
  mutate(gene_seq = paste0(gene_name," (",Sequence, ")")) %>% 
  filter(gene_seq %in% WCP_peptide_data_20250206_SN_diaNN_overlap_IDs) %>% 
  dplyr::select(-starts_with(c("log2_Intensity.", "Imputed_log2_Intensity."))) %>% 
  dplyr::select(gene_seq, 
                gene_name,
                Sequence,
                terminal_peptide,
                peptide_overJunction,
                NMD_peptide_reason,
                ANOVA_Imputed_log2_Intensity.06h_DMSO.vs.Imputed_log2_Intensity.06h_IAA,
                ANOVA_Imputed_log2_Intensity.12h_DMSO.vs.Imputed_log2_Intensity.12h_IAA,
                ANOVA_Imputed_log2_Intensity.15h_DMSO.vs.Imputed_log2_Intensity.15h_IAA,
                ANOVA_Imputed_log2_Intensity.18h_DMSO.vs.Imputed_log2_Intensity.18h_IAA,
                ANOVA_Imputed_log2_Intensity.24h_DMSO.vs.Imputed_log2_Intensity.24h_IAA,
  ) %>% 
  dplyr::rename("log2FC_ANOVA_06h_IAA" = "ANOVA_Imputed_log2_Intensity.06h_DMSO.vs.Imputed_log2_Intensity.06h_IAA",
                "log2FC_ANOVA_12h_IAA" = "ANOVA_Imputed_log2_Intensity.12h_DMSO.vs.Imputed_log2_Intensity.12h_IAA",
                "log2FC_ANOVA_15h_IAA" = "ANOVA_Imputed_log2_Intensity.15h_DMSO.vs.Imputed_log2_Intensity.15h_IAA",
                "log2FC_ANOVA_18h_IAA" = "ANOVA_Imputed_log2_Intensity.18h_DMSO.vs.Imputed_log2_Intensity.18h_IAA",
                "log2FC_ANOVA_24h_IAA" = "ANOVA_Imputed_log2_Intensity.24h_DMSO.vs.Imputed_log2_Intensity.24h_IAA") %>% 
  pivot_longer(cols = c(starts_with(c("log2FC_ANOVA"))),
               names_to = c(".value", "sample"),
               names_sep = "_ANOVA_") %>% 
  mutate(log2FC = -log2FC) %>% 
  mutate(condition = case_when(str_detect(sample, "12h_DMSO") ~ "12h_DMSO",
                               str_detect(sample, "12h_IAA") ~ "12h_IAA",
                               str_detect(sample, "24h_DMSO") ~ "24h_DMSO",
                               str_detect(sample, "24h_IAA") ~ "24h_IAA",
                               str_detect(sample, "06h_DMSO") ~ "06h_DMSO",
                               str_detect(sample, "06h_IAA") ~ "06h_IAA",
                               str_detect(sample, "15h_DMSO") ~ "15h_DMSO",
                               str_detect(sample, "15h_IAA") ~ "15h_IAA",
                               str_detect(sample, "18h_DMSO") ~ "18h_DMSO",
                               str_detect(sample, "18h_IAA") ~ "18h_IAA")) %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = (fct_relevel(NMD_peptide_reason_simple,
                                                  "AS_NMD",
                                                  "AS_NMD_UTR3",
                                                  "novel",
                                                  "uORF",
                                                  "overl_uORF",
                                                  "multiple",
                                                  "other")))


##### Stats -------------------------------------------------------------------
WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats <- WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  distinct(Sequence, gene_name) %>% 
  dplyr::count(gene_name) %>% 
  filter(n > 1) %>%
  arrange(desc(n))

WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats_single <- WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  distinct(Sequence, gene_name) %>% 
  dplyr::count(gene_name) %>% 
  filter(n == 1) %>%
  arrange(desc(n))


##### ANOVA log2FC ------------------------------------------------------

WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  # filter(!is.na(intensity)) %>% 
  dplyr::arrange(NMD_peptide_reason_simple, desc(log2FC)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(sample = (sample)) %>% 
  ggplot(aes(x=condition,
             y=gene_seq,
             fill=log2FC)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1),
        axis.text.y = element_text(size = 4),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       # limits = c(-5,5),
                       na.value = "grey80") +
  labs(y="peptides",
       x="") +
  ggside::geom_ysidetile(aes(x = "NMD reason", yfill = NMD_peptide_reason_simple))  +
  ggside::scale_yfill_manual(values=c("none" = "#214D65",
                                      "lncRNA" = "#76716E",
                                      "AS_NMD"="#E5BF85",
                                      "AS_NMD_UTR3" = "#C1AC65",
                                      "novel" = "#9A9A49",
                                      "uORF" = "#708831",
                                      "overl_uORF" = "#43761E",
                                      "multiple" = "gray30",
                                      "other" = "#D9DADA")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(100, "mm"),
                   cols = unit(80, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C13_Proteomics_NMD_peptides_ANOVA_log2FC_SN_afterModF.pdf"),
       width = 40,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 


##### Selection Plot ----------------------------------------------------------

###### Selection ---------------------------------------------------------------

WCP_peptide_data_20250127_SN_Selection_ANOVA <-  WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  filter(gene_name %in% WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats$gene_name) %>% 
  # filter(!is.na(intensity)) %>% 
  dplyr::arrange(NMD_peptide_reason_simple, gene_seq, desc(log2FC)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(Sequence = fct_rev(fct_inorder(Sequence))) %>% 
  mutate(sample = (sample)) %>% 
  mutate(condition = (fct_relevel(condition,
                                  "06h_IAA",
                                  "12h_IAA",
                                  "15h_IAA",
                                  "18h_IAA",
                                  "24h_IAA")))

# save 
save(WCP_peptide_data_20250127_SN_Selection_ANOVA,
     file="Resources/Proteomics/WCP_peptide_data_20250127_SN_Selection_ANOVA.rds")

WCP_peptide_data_20250127_SN_Selection_ANOVA_raw <- WCP_peptide_data_20250127_SN_final_log2_combined_Imputed_expand %>% 
  left_join(WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
              dplyr::select(log2FC, gene_seq, condition)) %>% 
  filter(gene_name %in% WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats$gene_name) %>% 
  # filter(!is.na(intensity)) %>% 
  dplyr::arrange(gene_seq, desc(log2FC)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(condition = (fct_relevel(condition,
                                  "06h_DMSO",
                                  "12h_DMSO",
                                  "15h_DMSO",
                                  "18h_DMSO",
                                  "24h_DMSO",
                                  "06h_IAA",
                                  "12h_IAA",
                                  "15h_IAA",
                                  "18h_IAA",
                                  "24h_IAA"))) %>% 
  group_by(condition, gene_seq) %>% 
  mutate(mean_log2_intensity = mean(log2_Intensity, na.rm=TRUE),
         valid_values = sum(!is.na(log2_Intensity))) %>% 
  ungroup()


###### Plot log2FC ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  ggplot(aes(x=condition,
             y=gene_seq,
             fill=log2FC)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  # geom_tile(color="black",
  #           linewidth=0.1) +
  geom_tile() +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       # limits = c(-5,5),
                       na.value = "grey80") +
  labs(y="peptides",
       x="") +
  ggside::geom_ysidetile(aes(x = "NMD reason", yfill = NMD_peptide_reason_simple))  +
  ggside::scale_yfill_manual(values=c("none" = "#214D65",
                                      "lncRNA" = "#76716E",
                                      "AS_NMD"="#E5BF85",
                                      "AS_NMD_UTR3" = "#C1AC65",
                                      "novel" = "#9A9A49",
                                      "uORF" = "#708831",
                                      "overl_uORF" = "#43761E",
                                      "multiple" = "gray30",
                                      "other" = "#D9DADA")) +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(5*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C14_Proteomics_AS_NMD_peptides_HQ_log2FC_imputed_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###### Plot raw log2 intensity ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Selection_ANOVA_raw %>% 
  # distinct(gene_seq, sample, log2_Intensity) %>% 
  ggplot(aes(x=condition,
             y=gene_seq)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile(fill = NA) +
  geom_point(aes(
    fill=mean_log2_intensity,
    size=valid_values),
    shape=21,
    stroke=0.1) +
  scale_size(range = c(0.75, 2)) +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(10*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C15_Proteomics_AS_NMD_peptides_HQ_raw_log2_intesity_imputed_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###### Plot peptide sequence ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  ggplot(aes(x="NMD_peptide_reason",
             y=Sequence,
             fill=NMD_peptide_reason_simple)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile(color="black",
            linewidth=0.1) +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "multiple" = "gray30",
                             "other" = "#D9DADA")) +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(1*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C16_Proteomics_AS_NMD_peptides_HQ_peptide_NMD_reason_imputed_Seq_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###### Plot tx log2FC DTE ---------------------------------------------------------------  

load("Resources/Proteomics/combined_peptidome_complete_class_junction.rds")

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  dplyr::select(-gene_name) %>% 
  left_join(combined_peptidome_complete_class_junction %>% 
              dplyr::select(gene_id, gene_name, transcript_id, sequence, NMD_tx_reason) %>% 
              dplyr::rename("Sequence" = "sequence")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h")) %>%
              dplyr::select(transcript_id, condition_2, logFC, FDR, logCPM) %>% 
              dplyr::rename("tx_log2FC" = "logFC",
                            "tx_FDR" = "FDR",
                            "tx_logCPM" = "logCPM")) %>% 
  group_by(gene_seq, condition_2) %>% 
  mutate(mean_tx_log2FC = mean(tx_log2FC, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x=mean_tx_log2FC,
             y=gene_seq,
             fill=condition_2)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_segment(aes(x=0, color=condition_2, xend=mean_tx_log2FC),
               alpha=0.75, show.legend=FALSE) +
  geom_point(shape=21,
             size=1.5,
             stroke=0.2,
             alpha=0.75) +
  scale_color_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                              "UPF1_Nter_4h" = "#229380",
                              "UPF1_Nter_12h" = "#00696D",
                              "UPF1_Nter_24h"="#005560")) +
  scale_fill_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                             "UPF1_Nter_4h" = "#229380",
                             "UPF1_Nter_12h" = "#00696D",
                             "UPF1_Nter_24h"="#005560")) +
  # geom_tile(color="black",
  #           linewidth=0.1) +
  # # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  # scale_fill_gradient2(low = ("#2166AC"),
  #                      mid = "white",
  #                      high = ("#B2182B"),
  #                      midpoint = 0,
  #                      # limits = c(-5,5),
  #                      na.value = "grey80") +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C17_Proteomics_AS_NMD_peptides_HQ_peptide_tx_log2FC_DTE_Seq_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###### Plot tx log2FC DGE ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  dplyr::select(-gene_name) %>% 
  left_join(combined_peptidome_complete_class_junction %>% 
              dplyr::select(gene_id, gene_name, transcript_id, sequence, NMD_tx_reason) %>% 
              dplyr::rename("Sequence" = "sequence")) %>% 
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h")) %>% 
              dplyr::select(gene_id, log2FoldChange, padj, baseMean, condition_2) %>% 
              dplyr::rename("gene_log2FC" = "log2FoldChange",
                            "gene_padj" = "padj",
                            "gene_baseMean" = "baseMean")) %>% 
  group_by(gene_seq, condition_2) %>% 
  mutate(mean_gene_log2FC = mean(gene_log2FC, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x=mean_gene_log2FC,
             y=gene_seq,
             fill=condition_2)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_segment(aes(x=0, color=condition_2, xend=mean_gene_log2FC),
               alpha=0.75, show.legend=FALSE) +
  geom_point(shape=21,
             size=1.5,
             stroke=0.2,
             alpha=0.75) +
  scale_color_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                              "UPF1_Nter_4h" = "#229380",
                              "UPF1_Nter_12h" = "#00696D",
                              "UPF1_Nter_24h"="#005560")) +
  scale_fill_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                             "UPF1_Nter_4h" = "#229380",
                             "UPF1_Nter_12h" = "#00696D",
                             "UPF1_Nter_24h"="#005560")) +
  # geom_tile(color="black",
  #           linewidth=0.1) +
  # # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  # scale_fill_gradient2(low = ("#2166AC"),
  #                      mid = "white",
  #                      high = ("#B2182B"),
  #                      midpoint = 0,
  #                      # limits = c(-5,5),
  #                      na.value = "grey80") +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C18_Proteomics_AS_NMD_peptides_HQ_peptide_gene_DGE_Seq_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

##
##### Conservation  ---------------------------------------------------------------  
##

load("Resources/Proteomics/WCP_peptide_data_20250127_SN_GRanges_combined_df.rds")

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250127_SN_GRanges_combined_df %>% 
              dplyr::select(gene_seq, protein_start, protein_end, phastScore)) %>% 
  ggplot(aes(x=phastScore,
             y=gene_seq)) +
  geom_point(aes(fill=peptide_overJunction),
             size=1.5,
             stroke=0.2,
             alpha=0.75,
             shape=21,
             color="black") +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "darkgray", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("TRUE" = "#963B5A",
                             "FALSE" = "#3b9779")) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  labs(y="",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(15, "mm")) 

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C19_Proteomics_AS_NMD_peptides_HQ_peptide_PhastCons_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")





##### DepMap KO ---------------------------------------------------------------  

# Filter DepMap data for gene_names in selection
WCP_peptide_data_20250127_SN_depmap <- depmap::depmap_crispr() %>% 
  filter(gene_name %in% WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats$gene_name)

WCP_peptide_data_20250127_SN_Selection_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250127_SN_depmap) %>% 
  ggplot(aes(x=dependency,
             y=gene_seq)) +
  geom_boxplot(linewidth=0.1,
               color="black",
               fill="#3181b4",
               outliers = FALSE) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  labs(y="",
       x="") +
  force_panelsizes(rows = unit(40*1.5, "mm"),
                   cols = unit(15, "mm")) 

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C20_Proteomics_AS_NMD_peptides_HQ_peptide_DepMap_CRISPR_2peptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


##### Single Peptides ----------------------------------------------------------

###### Single Peptide plots ---------------------------------------------------------------

WCP_peptide_data_20250127_SN_Single_ANOVA <-  WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  filter(gene_name %in% WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats_single$gene_name) %>% 
  # filter(!is.na(intensity)) %>% 
  dplyr::arrange(NMD_peptide_reason_simple, gene_seq, desc(log2FC)) %>% 
  mutate(gene_seq = fct_rev(fct_inorder(gene_seq))) %>% 
  mutate(gene_name = fct_rev(fct_inorder(gene_name))) %>% 
  mutate(Sequence = fct_rev(fct_inorder(Sequence))) %>% 
  mutate(sample = (sample)) %>% 
  mutate(condition = (fct_relevel(condition,
                                  "06h_IAA",
                                  "12h_IAA",
                                  "15h_IAA",
                                  "18h_IAA",
                                  "24h_IAA")))

###### Plot log2FC ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  ggplot(aes(x=condition,
             y=gene_seq,
             fill=log2FC)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  # geom_tile(color="black",
  #           linewidth=0.1) +
  geom_tile() +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       # limits = c(-5,5),
                       na.value = "grey80") +
  labs(y="peptides",
       x="") +
  ggside::geom_ysidetile(aes(x = "NMD reason", yfill = NMD_peptide_reason_simple))  +
  ggside::scale_yfill_manual(values=c("none" = "#214D65",
                                      "lncRNA" = "#76716E",
                                      "AS_NMD"="#E5BF85",
                                      "AS_NMD_UTR3" = "#C1AC65",
                                      "novel" = "#9A9A49",
                                      "uORF" = "#708831",
                                      "overl_uORF" = "#43761E",
                                      "multiple" = "gray30",
                                      "other" = "#D9DADA")) +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(5*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C21_Proteomics_AS_NMD_peptides_HQ_log2FC_imputed_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###### Plot peptide sequence ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  ggplot(aes(x="NMD_peptide_reason",
             y=Sequence,
             fill=NMD_peptide_reason_simple)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=4, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile(color="black",
            linewidth=0.1) +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "multiple" = "gray30",
                             "other" = "#D9DADA")) +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(1*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C22_Proteomics_AS_NMD_peptides_HQ_peptide_NMD_reason_imputed_Seq_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###### Plot peptide name ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  mutate(NMD_peptide_reason_simple = case_when(NMD_peptide_reason %in% c("none",
                                                                         "AS_NMD",
                                                                         "AS_NMD_UTR3",
                                                                         "novel",
                                                                         "uORF",
                                                                         "overl_uORF",
                                                                         "lncRNA",
                                                                         "multiple") ~ NMD_peptide_reason,
                                               TRUE ~ "other")) %>% 
  mutate(NMD_peptide_reason_simple = fct_rev(fct_relevel(NMD_peptide_reason_simple,
                                                         "AS_NMD",
                                                         "AS_NMD_UTR3",
                                                         "novel",
                                                         "uORF",
                                                         "overl_uORF",
                                                         "multiple",
                                                         "other"))) %>%
  ggplot(aes(x="NMD_peptide_reason",
             y=gene_name,
             fill=NMD_peptide_reason_simple)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile(color="black",
            linewidth=0.1) +
  # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "multiple" = "gray30",
                             "other" = "#D9DADA")) +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(1*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C22_Proteomics_AS_NMD_peptides_HQ_peptide_NMD_reason_imputed_Name_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###### Plot tx log2FC DTE ---------------------------------------------------------------  

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  dplyr::select(-gene_name) %>% 
  left_join(combined_peptidome_complete_class_junction %>% 
              dplyr::select(gene_id, gene_name, transcript_id, sequence, NMD_tx_reason) %>% 
              dplyr::rename("Sequence" = "sequence")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h")) %>%
              dplyr::select(transcript_id, condition_2, logFC, FDR, logCPM) %>% 
              dplyr::rename("tx_log2FC" = "logFC",
                            "tx_FDR" = "FDR",
                            "tx_logCPM" = "logCPM")) %>% 
  group_by(gene_seq, condition_2) %>% 
  mutate(mean_tx_log2FC = mean(tx_log2FC, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x=mean_tx_log2FC,
             y=gene_seq,
             fill=condition_2)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_segment(aes(x=0, color=condition_2, xend=mean_tx_log2FC),
               alpha=0.75, show.legend=FALSE) +
  geom_point(shape=21,
             size=1.5,
             stroke=0.2,
             alpha=0.75) +
  scale_color_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                              "UPF1_Nter_4h" = "#229380",
                              "UPF1_Nter_12h" = "#00696D",
                              "UPF1_Nter_24h"="#005560")) +
  scale_fill_manual(values=c("UPF1_Nter_0h" = "#65BB8C",
                             "UPF1_Nter_4h" = "#229380",
                             "UPF1_Nter_12h" = "#00696D",
                             "UPF1_Nter_24h"="#005560")) +
  # geom_tile(color="black",
  #           linewidth=0.1) +
  # # geom_text(aes(label=gene_name), color="black", size = 6*0.36, show_legend=FALSE) +
  # scale_fill_gradient2(low = ("#2166AC"),
  #                      mid = "white",
  #                      high = ("#B2182B"),
  #                      midpoint = 0,
  #                      # limits = c(-5,5),
  #                      na.value = "grey80") +
  labs(y="peptides",
       x="") +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C23_Proteomics_AS_NMD_peptides_HQ_peptide_tx_log2FC_DTE_Seq_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

##
###### Conservation  ---------------------------------------------------------------  
##

load("Resources/Proteomics/WCP_peptide_data_20250127_SN_GRanges_combined_df_single.rds")

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250127_SN_GRanges_combined_df_single %>% 
              dplyr::select(gene_seq, protein_start, protein_end, phastScore)) %>% 
  ggplot(aes(x=phastScore,
             y=gene_seq)) +
  geom_point(aes(fill=peptide_overJunction),
             size=1.5,
             stroke=0.2,
             alpha=0.75,
             shape=21,
             color="black") +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "darkgray", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("TRUE" = "#963B5A",
                             "FALSE" = "#3b9779")) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  labs(y="",
       x="") +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(15, "mm")) 

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C24_Proteomics_AS_NMD_peptides_HQ_peptide_PhastCons_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###### DepMap KO ---------------------------------------------------------------  

# Filter DepMap data for gene_names in selection
WCP_peptide_data_20250127_SN_depmap_single <- depmap::depmap_crispr() %>% 
  filter(gene_name %in% WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand_Stats_single$gene_name)

WCP_peptide_data_20250127_SN_Single_ANOVA %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  left_join(WCP_peptide_data_20250127_SN_depmap_single) %>% 
  ggplot(aes(x=dependency,
             y=gene_seq)) +
  geom_boxplot(linewidth=0.1,
               color="black",
               fill="#3181b4",
               outliers = FALSE) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        # axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=5, color="black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  labs(y="",
       x="") +
  force_panelsizes(rows = unit(47*1.5, "mm"),
                   cols = unit(15, "mm")) 

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C25_Proteomics_AS_NMD_peptides_HQ_peptide_DepMap_CRISPR_SinglePeptides.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")




##### Overall GO selected NMD peptides ----------------------------------------

# Prepare data
WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_forGO <- WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_expand %>% 
  distinct(Sequence, .keep_all = TRUE) %>% 
  dplyr::select(-gene_name) %>% 
  left_join(combined_peptidome_complete_class_junction %>% 
              dplyr::select(gene_id, gene_name, transcript_id, sequence, NMD_tx_reason) %>% 
              dplyr::rename("Sequence" = "sequence")) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Obtain background of expressed genes from degradation/recovery HCT116 data
Rev_1_F2_GO_DGE_bg_full <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c(
    "HCT116_UPF1_AID_degradation_this_Study",
    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Perform GO analysis
WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres <- gprofiler2::gost(query = WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_forGO,
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
WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres_result <- WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres$result

WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres_result %>% 
  mutate(query = "NMD peptides") %>% 
  arrange(p_value) %>% 
  mutate(term_name = fct_rev(fct_inorder(term_name))) %>% 
  mutate(term_id = (fct_inorder(term_id))) %>% 
  ggplot(aes(y=-log10(p_value),
             x=term_id,
             size=(intersection_size/query_size)*100,
             fill=query)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "darkgray", linewidth = 0.1),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.1, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.1,0.1,0.1,0.1)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2) +
  scale_size(range = c(0.25,3)) +
  coord_cartesian(clip = "off",
                  ylim = c(0,8)) +
  labs(size="matching\ngene IDs (%)",
       x="GO:BP",
       y="-log10(p.adjust)") +
  # scale_fill_identity() +
  scale_fill_manual(values=c("NMD peptides" = "#963B5A")) +
  guides(fill="none") +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(13*2, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_C26_WCP_peptide_data_gostres_result.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###### Enrichment visualization ----------------------------------------

simMatrix_WCP_peptide_data <- rrvgo::calculateSimMatrix(WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres_result %>% 
                                                          pull(term_id),
                                                        orgdb="org.Hs.eg.db",
                                                        ont="BP",
                                                        method="Rel")

scores_WCP_peptide_data <- setNames(-log10(WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres_result %>% 
                                             pull(p_value)), 
                                    WCP_peptide_data_20250127_SN_final_log2_combined_ANOVA_gostres_result %>% 
                                      pull(term_id))

reducedTerms_WCP_peptide_data <- rrvgo::reduceSimMatrix(simMatrix_WCP_peptide_data,
                                                        scores_WCP_peptide_data,
                                                        threshold=0.7,
                                                        orgdb="org.Hs.eg.db")

rrvgo::treemapPlot(reducedTerms_WCP_peptide_data)


###
## Rev_1 - Figure 7 - (D) Proteomics - Protein Group-level -----------------------------------------------------------
###

load("Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long.rds")

# Filter for significant hits
WCP_PG_data_20250127_GeneTxORF_long_sig <- WCP_PG_data_20250127_GeneTxORF_long %>% 
  filter(padj < 0.001 & abs(log2FC) > 0.5)

WCP_PG_data_20250127_GeneTxORF_long_sig_simple <- WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
  dplyr::select(gene_name, comparison, log2FC, padj)

##### Bar Plot --------------------------------------------------------------------

# Distribution of unique/multiple hits per ID %>% 
WCP_PG_data_20250127_GeneTxORF_long_sig %>% 
  distinct(comparison,id, .keep_all = TRUE) %>% 
  # Filter out UniProt
  filter(source != "UniProt") %>%
  mutate(UpDown = case_when(log2FC >0 ~ "up",
                            log2FC <0 ~ "down")) %>% 
  mutate(UpDown = fct_relevel(UpDown,
                              "up",
                              "down")) %>% 
  dplyr::count( comparison, protein_class, UpDown) %>% 
  # mutate(n = case_when(UpDown == "up" ~ n,
  #                      UpDown == "down" ~ -n)) %>% 
  mutate(protein_class = (fct_relevel(protein_class,
                                      "coding",
                                      "NMD",
                                      "mixed",
                                      "UniProt"))) %>% 
  mutate(comparison = (fct_relevel(comparison,
                                   "UPF1_06h",
                                   "UPF1_12h",
                                   "UPF1_15h",
                                   "UPF1_18h",
                                   "UPF1_24h"))) %>% 
  ggplot(aes(x=comparison,
             y=n,
             fill=protein_class)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_col(color="black",
           position = position_dodge(),
           linewidth = 0.1) +
  # geom_col(color="black",
  #          position = position_dodge2(width = 0.9, preserve = "single"),
  #          linewidth = 0.1) +
  geom_hline(yintercept = 0,
             color="#963B5A",
             linewidth=0.25
  ) +
  # geom_text(aes(label = case_when(protein_class == "coding" & abs(n) > 150 ~ paste0(n),
  #                                 TRUE ~ "")),
  #           color="white",
  #           size = 4*0.36,
  #           position = position_dodgenudge(x = -100, direction = "split", width = 1)
  # ) +
  # geom_text(aes(label = case_when(protein_class != "coding" ~ paste0(n),
  #                                 TRUE ~ "")),
  #           color="black",
  #           size = 4*0.36,
  #           position = position_dodge(width = 0.9)
# ) +
scale_fill_manual(values=c("coding" = "#214D65",
                           "mixed" = "#624B27",
                           "NMD" = "#E5BF86",
                           "UniProt" = "darkgray")) +
  coord_cartesian(clip = "off") +
  labs(x="",
       y="Proteins (n)",
       fill="Protein class") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(reverse=F)) +
  facet_wrap(~UpDown,ncol=1) +
  force_panelsizes(rows = unit(12.5, "mm"),
                   cols = unit(40, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_D1_Proteomics_ProteinGroups_SN_ProteinClass.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### GO -----------------------------------------------------------

load("Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO.rds")
load("Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO.rds")

#### Combined plot -----------------------------------------------------------
WCP_PG_data_GO_top_5_up <- WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO %>% 
  filter(parentTerm != "cellular process") %>% 
  group_by(parentTerm) %>% 
  mutate(sumScore = sum(-log10(p_value))) %>% 
  ungroup() %>%
  arrange(desc(sumScore)) %>% 
  dplyr::select(parentTerm, sumScore) %>%
  distinct() %>% 
  slice_max(sumScore, n=5)

WCP_PG_data_GO_top_5_down <- WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO %>% 
  group_by(parentTerm) %>% 
  mutate(sumScore = sum(-log10(p_value))) %>% 
  ungroup() %>%
  arrange(desc(sumScore)) %>% 
  dplyr::select(parentTerm, sumScore) %>%
  distinct() %>% 
  slice_max(sumScore, n=5)


WCP_PG_data_20250127_GeneTxORF_long_sig_up_GO %>% 
  mutate(UpDown="up") %>% 
  filter(parentTerm %in% WCP_PG_data_GO_top_5_up$parentTerm) %>%
  complete(query, UpDown, parentTerm) %>% 
  bind_rows(WCP_PG_data_20250127_GeneTxORF_long_sig_down_GO %>%
              mutate(UpDown="down") %>% 
              filter(parentTerm %in% WCP_PG_data_GO_top_5_down$parentTerm) %>% 
              complete(query, UpDown, parentTerm)) %>%
  bind_rows(tibble(query="UPF1_06h", UpDown="up", parentTerm="metabolic process")) %>% 
  mutate(score=-log10(p_value)) %>% 
  dplyr::rename("comparison" = "query") %>% 
  group_by(comparison, UpDown, parentTerm) %>% 
  mutate(mean_score=mean(score)) %>% 
  ungroup() %>% 
  # mutate(score = case_when(UpDown == "up" ~ score,
  #                          UpDown == "down" ~ -score)) %>%  
  # distinct(parentTerm, query, group_score) %>% 
  arrange(desc(mean_score)) %>% 
  mutate(parentTerm = fct_inorder(parentTerm)) %>% 
  # filter(parentTerm %in% c("response to stress",
  #                          "metabolic process",
  #                          "regulation of programmed cell death",
  #                          "biological regulation")) %>%
  mutate(comparison = (fct_relevel(comparison,
                                   "UPF1_06h",
                                   "UPF1_12h",
                                   "UPF1_15h",
                                   "UPF1_18h",
                                   "UPF1_24h"))) %>% 
  mutate(UpDown = fct_relevel(UpDown,
                              "up",
                              "down")) %>% 
  ggplot(aes(y=score,
             x=comparison,
             fill=parentTerm)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_bar(position = position_dodge(width=0.9,
                                     preserve = "single"),
           stat = "summary", 
           color="black",
           linewidth= 0.2,
           fun='mean',
           alpha=0.5) +
  geom_point(aes(group = parentTerm),
             shape=21,
             size=0.75,
             alpha=0.75,
             color="gray20",
             stroke=0.2,
             position = position_dodge(width=0.9, preserve = 'total')) +
  scale_fill_manual(values=c("metabolic process" = "#810f7c",
                             "biological regulation" = "#b10026",
                             "response to stress" = "#fc4e2a",
                             "regulation of programmed cell death" = "#feb24c",
                             "cellular localization" = "#ffffb2",
                             "response to stimulus" = "#034e7b",
                             "cell cycle" = "#3690c0",
                             "cellular component organization" = "#a6bddb",
                             "regulation of cellular process" = "#f1eef6")) +
  # geom_boxplot(outlier.colour=NA, 
  #              color="black",
  #              linewidth=0.2,
  #              position=position_dodge(width=0.9),
  #              alpha=0.75) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  guides(fill=guide_legend(reverse=F)) +
  facet_wrap(~UpDown,ncol=1) +
  force_panelsizes(rows = unit(12.5, "mm"),
                   cols = unit(40, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_D2_Proteomics_ProteinGroups_SN_GO_reduced.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 7 - (E) hominid-specific transcripts -----------------------------------------------------------
###

# Load external data
Suzuki_2018_HS_genes <- read_delim("Resources/External/Suzuki2018_PMID_29856955/Suzuki_2018_HS_genes.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

DESeq2_DGE_combined %>% 
  left_join(Suzuki_2018_HS_genes) %>% 
  filter(gene_family %in% c("CROCC",
                            "GUSB",
                            "NBPF",
                            "NOTCH2",
                            "NPIP",
                            "WASH") |
           str_detect(gene_name, 'TBC1D3') |
           gene_name %in% c("ARHGAP11A",
                            "ARHGAP11B")) %>% 
  dplyr::select(gene_name, gene_id) %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  arrange(desc(n))


# Plot
DESeq2_DGE_combined %>% 
  left_join(Suzuki_2018_HS_genes) %>% 
  filter(gene_family %in% c("CROCC",
                            "GUSB",
                            "NBPF",
                            "NOTCH2",
                            "NPIP",
                            "WASH") |
           str_detect(gene_name, 'TBC1D3') |
           gene_name %in% c("ARHGAP11A",
                            "ARHGAP11B")) %>% 
  # remove false gene_ids
  filter(!gene_id %in% c("ENSG00000215908.12",
                         "ENSG00000290842.1",
                         "ENSG00000183666.20",
                         "ENSG00000241549.9",
                         "ENSG00000272150.6")) %>% 
  mutate(gene_family = case_when(str_detect(gene_name, 'TBC1D3') ~ "TBC1D",
                                 str_detect(gene_name, 'ARHGAP11') ~ "ARHGAP11",
                                 TRUE ~ gene_family)) %>% 
  filter(selected_family == TRUE | str_detect(gene_name, 'TBC1D3') | str_detect(gene_name, 'ARHGAP11')) %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "GUSB",
                                 "GUSBP1",
                                 "GUSBP2",
                                 "GUSBP4",
                                 "GUSBP9",
                                 "GUSBP11",
                                 "NBPF1",
                                 "NBPF3",
                                 "NBPF8",
                                 "NBPF9")) %>% 
  ggplot(aes(x=gene_name,
             y=condition_2
  )) +
  geom_tile(aes(fill = case_when(padj < 0.0001 & abs(log2FoldChange) > 1  ~ "white",
                                 TRUE ~ "gray90")),
            color = "white",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    stroke=0.2,
    shape=21) +
  scale_size(range = c(0.5, 2)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey90") +
  facet_wrap(~gene_family, scales="free_x", nrow=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(38.596, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_E1_Hominid_specific_genes_DGE.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### bakR-----------------------------------------------------------
DESeq2_DGE_combined %>% 
  left_join(Suzuki_2018_HS_genes) %>% 
  filter(gene_family %in% c("CROCC",
                            "GUSB",
                            "NBPF",
                            "NOTCH2",
                            "NPIP",
                            "WASH") |
           str_detect(gene_name, 'TBC1D3') |
           gene_name %in% c("ARHGAP11A",
                            "ARHGAP11B")) %>% 
  # remove false gene_ids
  filter(!gene_id %in% c("ENSG00000215908.12",
                         "ENSG00000290842.1",
                         "ENSG00000183666.20",
                         "ENSG00000241549.9",
                         "ENSG00000272150.6")) %>% 
  mutate(gene_family = case_when(str_detect(gene_name, 'TBC1D3') ~ "TBC1D",
                                 str_detect(gene_name, 'ARHGAP11') ~ "ARHGAP11",
                                 TRUE ~ gene_family)) %>% 
  filter(selected_family == TRUE | str_detect(gene_name, 'TBC1D3') | str_detect(gene_name, 'ARHGAP11')) %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "GUSB",
                                 "GUSBP1",
                                 "GUSBP2",
                                 "GUSBP4",
                                 "GUSBP9",
                                 "GUSBP11",
                                 "NBPF1",
                                 "NBPF3",
                                 "NBPF8",
                                 "NBPF9")) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id, L2FC_kdeg, padj_kdeg)) %>% 
  filter(condition_2 == "UPF1_Nter_12h") %>% 
  ggplot(aes(x=gene_name,
             y=condition_2
  )) +
  geom_tile(aes(fill = L2FC_kdeg),
            color = "white",
            lwd = 0.1,
            linetype = 1) +
  scale_size(range = c(0, 2)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey90") +
  facet_wrap(~gene_family, scales="free_x", nrow=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(2, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure7", "Rev_1_F7_E1_Hominid_specific_genes_kdeg.pdf"),
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
