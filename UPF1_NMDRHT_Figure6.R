#!/usr/bin/env Rscript

# Title: UPF1_Revision_Figure6
# Objective: Code for generating Panels of Figure 6 + Figure S6 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_Revision_Analysis"

###
## Rev_1 - Figure 6 - (A) Early Up Transcripts -----------------------------------------------------------
###

##### early Up - logFC vs NMD_reason ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  mutate(NMD_tx_status_simple = case_when(NMD_tx_status %in% c("coding", "predicted_coding") ~ "coding",
                                          NMD_tx_status %in% c("NMD", "predicted_NMD", "mixed") ~ "NMD",
                                          NMD_tx_status %in% c("lncRNA") ~ "lncRNA")) %>% 
  mutate(NMD_tx_status = fct_rev(fct_relevel(NMD_tx_status,
                                             "coding",
                                             "predicted_coding",
                                             "NMD",
                                             "predicted_NMD",
                                             "mixed",
                                             "lncRNA"))) %>%
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  ggplot(aes(x=logFC,
             y=NMD_tx_status,
             fill=NMD_tx_status)) +
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
  ggridges::geom_density_ridges(
    jittered_points = TRUE,
    scale=0.8,
    rel_min_height = 0.005,
    size=0.1,
    position = ggridges::position_points_jitter(width = 0.05, height = 0, yoffset=-0.1),
    point_shape = '|', point_size = 0.25, point_alpha = 0.25, alpha = 0.7,
  ) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD"="#E5BF86",
                             "mixed" = "#624B27",
                             "predicted_NMD" = "#B09771",
                             "predicted_coding" = "#287DAB")) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_A1_12h_UPF1_NMD_status_Ridgeline.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### early Up - all - bakR ------------------------------------------------------------
# Join with gene-level analyses for bakR parameter
NMDRHT.v1.2_MainTable_bakR_gene_all <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, baseMean)) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                      L2FC_kdeg,
                                                      padj_kdeg,
                                                      L2FC_ksyn,
                                                      padj_ksyn,
                                                      Mech_score,
                                                      kdeg_conclusion,
                                                      RNA_conclusion,
                                                      Mech_conclusion))

# How many transcripts with which mechanism conclusion
NMDRHT.v1.2_MainTable_bakR_gene_all %>% 
  mutate(bakR_detected = case_when(is.na(Mech_score) ~ FALSE,
                                   !is.na(Mech_score) ~ TRUE)) %>% 
  dplyr::count(NMD_50nt_rule, Mech_conclusion) %>% 
  group_by(NMD_50nt_rule) %>% 
  mutate(n_per = n/sum(n))

# Plot
NMDRHT.v1.2_MainTable_bakR_gene_all %>% 
  filter(!is.na(L2FC_kdeg)) %>% 
  pivot_longer(cols=c(L2FC_kdeg, L2FC_ksyn),
               values_to = "L2FC",
               names_to = "rate") %>% 
  ggplot(aes(x=Mech_conclusion,
             y=L2FC,
  )) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  geom_violin(aes(fill=rate),
              alpha=0.5,
              linewidth=0.1,
              position=position_dodge(0.7), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               position=position_dodge(0.7),
               width=0.15,
               linewidth=0.1,
               color="black",
               aes(fill=rate)) +
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d")) +
  facet_wrap(~NMD_50nt_rule) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_A2_DTE_up_all_bakR.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### early Up - no NMD - bakR ------------------------------------------------------------
# Join with gene-level analyses for bakR parameter
NMDRHT.v1.2_MainTable_bakR_gene <- NMDRHT.v1.2_MainTable %>%  
  filter(DTE_cluster == "up 1:early" & NMD_tx_reason %in% c("none", "lncRNA")) %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, baseMean)) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  arrange(FDR)

# 1238 transcripts with bakR parameter! 285 without
NMDRHT.v1.2_MainTable_bakR_gene %>% 
  mutate(bakR_detected = case_when(is.na(Mech_score) ~ FALSE,
                                   !is.na(Mech_score) ~ TRUE)) %>% 
  dplyr::count(NMD_tx_reason, Mech_conclusion) %>% 
  group_by(NMD_tx_reason) %>% 
  mutate(n_per = n/sum(n))

# Plot
NMDRHT.v1.2_MainTable_bakR_gene %>% 
  filter(!is.na(L2FC_kdeg)) %>% 
  pivot_longer(cols=c(L2FC_kdeg, L2FC_ksyn),
               values_to = "L2FC",
               names_to = "rate") %>% 
  ggplot(aes(x=Mech_conclusion,
             y=L2FC,
  )) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  geom_violin(aes(fill=rate),
              alpha=0.5,
              linewidth=0.1,
              position=position_dodge(0.7), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               position=position_dodge(0.7),
               width=0.15,
               linewidth=0.1,
               color="black",
               aes(fill=rate)) +
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d")) +
  facet_wrap(~NMD_tx_reason) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_A3_DTE_up_no_NMD_bakR.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")



####### early Up - coding - GO ------------------------------------------------------------
# Define background: all expressed NMD reason = "none" (~ coding) with ENSEMBL IDs
Rev_1_F6_DTE_up_no_NMD_bg <- edgeR_DTE_NMDRHT_combined %>% 
  filter(experimentSet %in% c(
    "HCT116_UPF1_AID_degradation_this_Study",
    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  left_join(NMDRHT.v1.2_MainTable %>% dplyr::select(gene_id, NMD_tx_reason)) %>% 
  filter(NMD_tx_reason == "none") %>% 
  filter(str_detect(gene_id, "^ENSG") == TRUE) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Make lists
Rev_1_F6_DTE_up_no_NMD_geneID <- split(NMDRHT.v1.2_MainTable_bakR_gene %>% 
                                         separate(gene_id, c("gene_id", "Version")) %>% 
                                         arrange(desc(logFC)) %>% 
                                         distinct(gene_id, .keep_all = TRUE) %>% 
                                         pull(gene_id),
                                       NMDRHT.v1.2_MainTable_bakR_gene %>% 
                                         arrange(desc(logFC)) %>% 
                                         distinct(gene_id, .keep_all = TRUE) %>% 
                                         mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "NA",
                                                                            TRUE ~ Mech_conclusion)) %>% 
                                         pull(Mech_conclusion)) 

Rev_1_F6_DTE_up_no_NMD_gostres <- gprofiler2::gost(query = Rev_1_F6_DTE_up_no_NMD_geneID,
                                       custom_bg = Rev_1_F6_DTE_up_no_NMD_bg,
                                       ordered_query = TRUE,
                                       # multi_query = TRUE,
                                       sources = "GO:BP",
                                       domain_scope = "custom",
                                       organism = "hsapiens",
                                       correction_method = c("gSCS"),
                                       # evcodes = TRUE,
                                       # as_short_link = TRUE,
                                       significant = FALSE)

Rev_1_F6_DTE_up_no_NMD_gostres_result <- Rev_1_F6_DTE_up_no_NMD_gostres$result

Rev_1_F6_DTE_up_no_NMD_gostres_result %>% 
  mutate(query = (fct_relevel(query,
                              "Degradation",
                              "Not Sig.",
                              "Synthesis",
                              "NA"
  ))) %>% 
  # filter(p_value != 1) %>% 
  ggplot(aes(y=-log10(p_value),
             x=query,
             color=query)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
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
  geom_boxplot(outlier.size = 0.5,
               outlier.alpha = 0.5) +
  scale_fill_identity() +
  scale_color_manual(values=c("Degradation" = "#018571",
                              "Not Sig." = "#4D4D4D",
                              "Synthesis" = "#A6611A",
                              "NA" = "#5B5B5B")) +
  geom_hline(yintercept = -log10(0.05)) +
  guides(fill="none") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_A4_DTE_up_no_NMD_gostres_result.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

####### Synth GO terms  ---------------------------------------------

simMatrix_F6_DTE_up_no_NMD <- rrvgo::calculateSimMatrix(Rev_1_F6_DTE_up_no_NMD_gostres_result %>% 
                                                   filter(query == "Synthesis") %>% 
                                                   filter(p_value < 0.05) %>% 
                                                   pull(term_id),
                                                 orgdb="org.Hs.eg.db",
                                                 ont="BP",
                                                 method="Rel")

scores_F6_DTE_up_no_NMD <- setNames(-log10(Rev_1_F6_DTE_up_no_NMD_gostres_result %>% 
                                             filter(query == "Synthesis") %>% 
                                             filter(p_value < 0.05) %>% 
                                             pull(p_value)), Rev_1_F6_DTE_up_no_NMD_gostres_result %>% 
                                      filter(query == "Synthesis") %>% 
                                      filter(p_value < 0.05) %>% 
                                      pull(term_id))

reducedTerms_F6_DTE_up_no_NMD <- rrvgo::reduceSimMatrix(simMatrix_F6_DTE_up_no_NMD,
                                                 scores_F6_DTE_up_no_NMD,
                                                 threshold=0.7,
                                                 orgdb="org.Hs.eg.db")

rrvgo::scatterPlot(simMatrix_F6_DTE_up_no_NMD, reducedTerms_F6_DTE_up_no_NMD)
rrvgo::treemapPlot(reducedTerms_F6_DTE_up_no_NMD)
rrvgo::wordcloudPlot(reducedTerms_F6_DTE_up_no_NMD, min.freq=1, colors="black")


##
#### Tx-level EZbakR ---------------------------------------------------------
##

##### early Up - all - EZbakR-tx ------------------------------------------------------------

# Load EZbakR analysis of transcript-level NMDRHT-based
UPF1_NMDRHT_EZbakR_TEC_combined <- read_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))

# Join with gene-level analyses for EZbakR parameter
NMDRHT.v1.2_MainTable_EZbakR_tx_all <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, baseMean)) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>%
              filter(condition == "UPF1_12h") %>%
              dplyr::select(transcript_id,
                            L2FC_kdeg_tx,
                            padj_kdeg_tx,
                            kdeg_tx_conclusion))

# How many transcripts with which mechanism conclusion
NMDRHT.v1.2_MainTable_EZbakR_tx_all %>% 
  mutate(EZbakR_detected = case_when(is.na(kdeg_tx_conclusion) ~ FALSE,
                                   !is.na(kdeg_tx_conclusion) ~ TRUE)) %>% 
  dplyr::count(NMD_50nt_rule, kdeg_tx_conclusion) %>% 
  group_by(NMD_50nt_rule) %>% 
  mutate(n_per = n/sum(n))

# Plot
NMDRHT.v1.2_MainTable_EZbakR_tx_all %>% 
  filter(!is.na(L2FC_kdeg)) %>% 
  pivot_longer(cols=c(L2FC_kdeg, L2FC_ksyn, L2FC_kdeg_tx),
               values_to = "L2FC",
               names_to = "rate") %>% 
  mutate(rate = fct_relevel(rate,
                            "L2FC_kdeg",
                            "L2FC_ksyn",
                            "L2FC_kdeg_tx")) %>% 
  ggplot(aes(x=Mech_conclusion,
             y=L2FC,
  )) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  geom_violin(aes(fill=rate),
              alpha=0.5,
              linewidth=0.1,
              position=position_dodge(0.8), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               position=position_dodge(0.8),
               width=0.15,
               linewidth=0.1,
               color="black",
               aes(fill=rate)) +
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d",
                             "L2FC_kdeg_tx" = "#cd7eb2")) +
  facet_wrap(~NMD_50nt_rule) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_A2_DTE_up_all_bakR_EZbakR.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (AA) Bambu FL -----------------------------------------------------------
###

### Barplot -----------------------------------------------------------------

# Quantification performed on NMDRHT.v1.1 -> load information
qlf_df_NMDRHT_combined <- read_csv("Resources/NMDRHT/Bambu/Bambu_NMDRHT_qlf_df_FLreads.csv")

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>% 
  mutate(significant = case_when(FDR < 0.0001 & abs(logFC) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(logFC>0 ~ "up",
                             logFC<0 ~ "down")) %>% 
  left_join(qlf_df_NMDRHT_combined %>% 
              dplyr::select(transcript_id, logFC, FDR, condition_2) %>% 
              dplyr::rename("logFC_LR" = "logFC",
                            "FDR_LR" = "FDR",
                            "condition_2_LR" = "condition_2") %>% 
              filter(condition_2_LR == "NMDRHT_ONT_dRNA_12h_IAA")) %>% 
  mutate(logFC_LR = case_when(is.na(logFC_LR) ~ 0,
                              TRUE ~ logFC_LR),
         FDR_LR = case_when(is.na(FDR_LR) ~ 1,
                            TRUE ~ FDR_LR),
         condition_2_LR = case_when(is.na(condition_2_LR) ~ "NMDRHT_ONT_dRNA_12h_IAA",
                                    TRUE ~ condition_2_LR)) %>% 
  bind_rows(NMDRHT.v1.2_MainTable %>% 
              left_join(edgeR_DTE_NMDRHT_combined %>% 
                          filter(condition_2 %in% c("UPF1_Nter_12h"))) %>% 
              mutate(significant = case_when(FDR < 0.0001 & abs(logFC) > 1 ~ TRUE,
                                             TRUE ~ FALSE)) %>% 
              mutate(up_down = case_when(logFC>0 ~ "up",
                                         logFC<0 ~ "down")) %>% 
              left_join(qlf_df_NMDRHT_combined %>% 
                          dplyr::select(transcript_id, logFC, FDR, condition_2) %>% 
                          dplyr::rename("logFC_LR" = "logFC",
                                        "FDR_LR" = "FDR",
                                        "condition_2_LR" = "condition_2") %>% 
                          filter(condition_2_LR == "NMDRHT_PacBio_12h_IAA")) %>% 
              mutate(logFC_LR = case_when(is.na(logFC_LR) ~ 0,
                                          TRUE ~ logFC_LR),
                     FDR_LR = case_when(is.na(FDR_LR) ~ 1,
                                        TRUE ~ FDR_LR),
                     condition_2_LR = case_when(is.na(condition_2_LR) ~ "NMDRHT_PacBio_12h_IAA",
                                                TRUE ~ condition_2_LR))) %>% 
  mutate(significant_LR = case_when(FDR_LR < 0.01 & abs(logFC_LR) > 1 ~ TRUE,
                                    TRUE ~ FALSE)) %>% 
  mutate(up_down_LR = case_when(logFC_LR>0 ~ "up",
                                logFC_LR<0 ~ "down")) %>% 
  mutate(NMD_tx_status = fct_rev(fct_relevel(NMD_tx_status,
                                             "coding",
                                             "predicted_coding",
                                             "NMD",
                                             "predicted_NMD",
                                             "mixed",
                                             "lncRNA"))) %>%
  filter(significant==TRUE) %>% 
  mutate(comb_LR = paste0(significant_LR,"_",up_down_LR)) %>%
  dplyr::count(up_down, comb_LR, condition_2_LR, NMD_tx_status) %>% 
  arrange(desc(n)) %>% 
  group_by(up_down, NMD_tx_status, condition_2_LR) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  mutate(up_down = fct_relevel(up_down,
                               "up",
                               "down")) %>% 
  mutate(condition_2_LR = fct_relevel(condition_2_LR,
                                      "NMDRHT_ONT_dRNA_12h_IAA",
                                      "NMDRHT_PacBio_12h_IAA")) %>% 
  mutate(comb_LR = fct_relevel(comb_LR,
                               "TRUE_up",
                               "FALSE_up",
                               "FALSE_NA",
                               "FALSE_down",
                               "TRUE_down")) %>% 
  # mutate(n = case_when(up_down == "up" ~ n,
  #                      up_down == "down"  ~ -n)) %>% 
  ggplot(aes(y=NMD_tx_status,
             x=n_per,
             fill=comb_LR)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  # geom_pointrange(aes(x=median,
  #                     xmin=Q1,
  #                     xmax=Q3),
  #                 position = position_dodge(width = 0.95),
  #                 linewidth = 0.25,
  #                 size = 0.4,
  #                 stroke=0.1,
  #                 alpha=0.75,
  #                 shape=21) +
  geom_col(color="black",
           linewidth=0.1) +
  geom_text(aes(label = case_when(n_per > 40 & comb_LR %in% c("TRUE_up", "TRUE_down") ~ paste0(n,
                                                                                               " (",
                                                                                               round(n_per,0),
                                                                                               "%)"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 40 & !comb_LR %in% c("TRUE_up", "TRUE_down") ~ paste0(n,
                                                                                                " (",
                                                                                                round(n_per,0),
                                                                                                "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("TRUE_up" = "#B2182B",
                             "FALSE_up" = "#FFC08E",
                             "FALSE_NA" = "gray80",
                             "FALSE_down" = "#B7C2DE",
                             "TRUE_down" = "#2166AC")) +
  # theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~up_down+condition_2_LR,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(18, "mm"),
                   cols = unit(18, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_AA_sig_SR_edgeR_longRead_FLreads_NMD_status_barplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Boxplot -----------------------------------------------------------------

qlf_df_NMDRHT_combined %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=NMD_tx_status,
             y=logFC,
             fill=condition_2)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  geom_hline(yintercept=0,
             color="#898989",
             linewidth=0.25,
             linetype="solid") +
  # geom_pointrange(aes(x=median,
  #                     xmin=Q1,
  #                     xmax=Q3),
  #                 position = position_dodge(width = 0.95),
  #                 linewidth = 0.25,
  #                 size = 0.4,
  #                 stroke=0.1,
  #                 alpha=0.75,
  #                 shape=21) +
  geom_boxplot(outliers = FALSE,
               # coef = 0,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("NMDRHT_ONT_dRNA_12h_IAA" = "#A2D9F7",
                             "NMDRHT_PacBio_12h_IAA" = "#A597DB")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(20+0.4, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_AA_log2FC_edgeR_longRead_FLreads_NMD_status_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (B) long 3UTR NMD -----------------------------------------------------------
###

###### early Up - UTR3 length - Mechanism ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>%
  filter(NMD_tx_reason == "none") %>% 
  filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "NA",
                                     TRUE ~ Mech_conclusion)) %>% 
  mutate(Mech_conclusion = fct_rev(fct_relevel(Mech_conclusion,
                                               "Degradation",
                                               "Not Sig.",
                                               "Synthesis",
                                               "NA"
  ))) %>% 
  ggplot(aes(y=DTE_cluster,
             x=utr3_len,
             fill=Mech_conclusion)) +
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
  geom_boxplot(outliers = FALSE,
               linewidth=0.1,
               width = .8, alpha = .5, show.legend = TRUE) +
  coord_cartesian(xlim = c(0,5000)) +
  scale_fill_manual(values=c("Degradation" = "#018571",
                             "Not Sig." = "#4D4D4D",
                             "Synthesis" = "#A6611A",
                             "NA" = "#5B5B5B")) +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_B1_DTE_up_no_NMD_UTR3_length_boxplot.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###### early Up - UTR3 length ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>%
  filter(NMD_tx_reason == "none") %>% 
  filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "NA",
                                     TRUE ~ Mech_conclusion)) %>% 
  mutate(Mech_conclusion = fct_rev(fct_relevel(Mech_conclusion,
                                               "Degradation",
                                               "Not Sig.",
                                               "Synthesis",
                                               "NA"
  ))) %>% 
  ggplot(aes(y=DTE_cluster,
             x=utr3_len,
             color=DTE_cluster)) +
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
  ggdist::stat_pointinterval(alpha = .5, 
                     size=0.2,
                     interval_linewidth = 0.25,
                     stroke=0,
                     point_size = 1.5,
                     show.legend = TRUE) +
  coord_cartesian(xlim = c(0,5000)) +
  scale_color_manual(values=c("up 1:early" = "#67001f",
                              "up 2:delayed" = "#b2182b",
                              "up 3:late" = "#d6604d",
                              "down 3:late" = "#92c5de",
                              "down 2:delayed" = "#4393c3",
                              "down 1:early" = "#053061",
                              "expressed" = "gray30")) +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_B2_DTE_up_no_NMD_UTR3_length_boxplot_woMech.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

#### Stats -------------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>%
  filter(NMD_tx_reason == "none") %>% 
  filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster, 
                                    "up 1:early",
                                    "up 2:delayed",
                                    "up 3:late",
                                    "expressed",
                                    "down 3:late",
                                    "down 2:delayed",
                                    "down 1:early"))) %>% 
  group_by(DTE_cluster) %>% 
  rstatix::get_summary_stats(utr3_len) 

# Check Quantile-Quantile plot
for_plot <- NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>%
  filter(NMD_tx_reason == "none") %>% 
  filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster, 
                                    "up 1:early",
                                    "up 2:delayed",
                                    "up 3:late",
                                    "expressed",
                                    "down 3:late",
                                    "down 2:delayed",
                                    "down 1:early")))

ggpubr::ggqqplot(for_plot, x = "utr3_len",
         color = "DTE_cluster", 
         # palette = c("#0073C2FF", "#FC4E07"),
         ggtheme = ggpubr::theme_pubclean())

# Perform non-parametric analog of the one-way ANOVA -> Dunn's Test
for_plot %>% rstatix::dunn_test(utr3_len ~ DTE_cluster) %>% 
  filter(group1 == "expressed" | group2 == "expressed")

#### Scatter log2FC UTR3 length ------------------------------------------------------------

# Check for overlap at 12h in HCT N-AID-UPF1, HCT FKBP-UPF1, HEK FKBP-UPF1 and long-read RNA-Seq
NMDRHT.v1.2_MainTable_HQ_up_early <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>%
  bind_rows(qlf_df_NMDRHT_combined %>% 
              filter(DTE_cluster == "up 1:early"))

NMDRHT.v1.2_MainTable_HQ_up_early %>% 
  mutate(condition_2 = fct_relevel(condition_2,
                                   "UPF1_Nter_12h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_FKBP_HEK_12h",
                                   "NMDRHT_ONT_dRNA_12h_IAA", "NMDRHT_PacBio_12h_IAA"
  )) %>% 
  filter(NMD_tx_reason == "none") %>% 
  mutate(n_sig_UPF1 = sum(FDR < 0.01)) %>% 
  ggplot(aes(x=log10(utr3_len),
             y=logFC)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 6),
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
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.25),
                     dpi=1200, dev = "cairo") +
  scale_color_viridis_c(option="mako") +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                 alpha=0.75,
                 color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2")),
               size = 5*0.30,
               color="gray20",
               label.y = "bottom", label.x = "right") +
  facet_wrap(~condition_2, nrow=1)  +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_B3_utr3_len_logFC_pointdensity.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

###### early Up - NMD bin _ 3UTR_length ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  filter(NMD_tx_reason == "none") %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  ggplot(aes(x=NMD_bin,
             y=log10(utr3_len))) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  geom_violin(aes(fill=NMD_bin),
              alpha=0.5,
              linewidth=0.1,
              position=position_dodge(0.7), show.legend = FALSE) +
  geom_boxplot(outliers=FALSE,
               position=position_dodge(0.7),
               width=0.15,
               linewidth=0.1,
               color="black",
               aes(fill=NMD_bin)) +
  scale_fill_manual(values=c("(75,100]" = "#E5BF86",
                             "(50,75]" = "#B09771",
                             "(25,50]" = "#647588",
                             "[0,25]" = "#CACFD0")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_B4_DTE_up_no_NMD_bin_UTR3_len.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Stats
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  filter(NMD_tx_reason == "none") %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>%
  # group_by(NMD_bin) %>% 
  rstatix::dunn_test(utr3_len ~ NMD_bin)

###
## Rev_1 - Figure 6 - (C) SRD -----------------------------------------------------------
###

###### SRD - early Up - DTE cluster ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed","down 3:late",
                             "down 2:delayed",
                             "down 1:early")) %>%
  # filter(NMD_tx_reason == "none") %>%
  filter(!is.na(utr3_mfe_nt) & utr3_len != 0) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed"))) %>% 
  ggplot(aes(y=DTE_cluster,
             x=-utr3_mfe_nt,
             color=NMD_50nt_rule)) +
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
  ggdist::stat_halfeye(alpha = .5, 
               size=1,
               stroke=0,
               position = position_dodge(0.9),
               point_size = 2,
               show.legend = TRUE) +
  scale_color_manual(values=c("TRUE" = "#B09771",
                              "FALSE" = "#214D65")) +
  # scale_fill_manual(values=c("up 1:early" = "#67001f",
  #                             "up 2:delayed" = "#b2182b",
  #                             "up 3:late" = "#d6604d",
  #                             "down 3:late" = "#92c5de",
  #                             "down 2:delayed" = "#4393c3",
  #                             "down 1:early" = "#053061",
  #                             "expressed" = "gray30")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(30, "mm")) +
  guides(fill=guide_legend(reverse=T))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C1_DTE_up_no_NMD_structure_DTE_cluster.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Summary stats
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed",
                             "down 3:late",
                             "down 2:delayed",
                             "down 1:early")) %>%
  # filter(NMD_tx_reason == "none") %>%
  mutate(utr3_mfe_nt = -utr3_mfe_nt) %>% 
  filter(!is.na(utr3_mfe_nt) & utr3_len != 0) %>% 
  # filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed"))) %>% 
  group_by(DTE_cluster, NMD_50nt_rule) %>% 
  rstatix::get_summary_stats(utr3_mfe_nt)

# Dunn's Test
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed",
                             "down 3:late",
                             "down 2:delayed",
                             "down 1:early")) %>%
  # filter(NMD_tx_reason == "none") %>%
  mutate(utr3_mfe_nt = -utr3_mfe_nt) %>% 
  filter(!is.na(utr3_mfe_nt) & utr3_len != 0) %>% 
  # filter(utr3_len != 0) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed"))) %>% 
  group_by(NMD_50nt_rule) %>% 
  rstatix::dunn_test(utr3_mfe_nt ~ DTE_cluster) %>% 
  filter(group1 == "expressed" | group2 == "expressed")

####### SRD - Scatter Plot ------------------------------------------------------------
NMDRHT.v1.2_MainTable_HQ_up_early %>% 
  mutate(condition_2 = fct_relevel(condition_2,
                                   "UPF1_Nter_12h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_FKBP_HEK_12h",
                                   "NMDRHT_ONT_dRNA_12h_IAA", "NMDRHT_PacBio_12h_IAA"
  )) %>% 
  filter(NMD_tx_reason == "none") %>% 
  mutate(n_sig_UPF1 = sum(FDR < 0.01)) %>% 
  ggplot(aes(x=-utr3_mfe_nt,
             y=logFC)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 6),
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
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.25),
                     dpi=1200, dev = "cairo") +
  scale_color_viridis_c(option="mako") +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                 alpha=0.75,
                 color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2")),
               size = 5*0.30,
               color="gray20",
               label.y = "bottom", label.x = "right") +
  facet_wrap(~condition_2, nrow=1)  +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C2_utr3_structure_logFC_pointdensity.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### SRD - HSU-vs-PSU ----------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(utr3_mfe_nt_class = case_when(-utr3_mfe_nt>=0.3 ~ "HSU",
                                       -utr3_mfe_nt<=0.2 ~ "PSU",
                                       TRUE ~ "other"), .after=utr3_mfe_nt) %>% 
  mutate(utr3_mfe_nt_class = fct_relevel(utr3_mfe_nt_class,
                                         "HSU",
                                         "other",
                                         "PSU")) %>% 
  group_by(condition_2, utr3_mfe_nt_class, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=utr3_mfe_nt_class,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=paste0("-","\u0394","G/nt (bins)"),
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(3*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C3_utr3_structure_HSU_PSU_median_log2FC_global_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### SRD - HSU-vs-PSU up early ----------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(utr3_mfe_nt_class = case_when(-utr3_mfe_nt>=0.3 ~ "HSU",
                                       -utr3_mfe_nt<=0.2 ~ "PSU",
                                       TRUE ~ "other"), .after=utr3_mfe_nt) %>% 
  mutate(utr3_mfe_nt_class = fct_relevel(utr3_mfe_nt_class,
                                         "HSU",
                                         "other",
                                         "PSU")) %>% 
  group_by(condition_2, utr3_mfe_nt_class, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=utr3_mfe_nt_class,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=paste0("-","\u0394","G/nt (bins)"),
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(3*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C4_utr3_structure_HSU_PSU_median_log2FC_up_early_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf)

#### SRD - UPF1 helicase-dependent targets -----------------------------------------
# Import TableS1 from PMID: 32017897
Fischer_2020_TableS1 <- read_excel("Resources/External/Fischer2020_PMID_32017897/1-s2.0-S1097276520300423-mmc2.xlsx") %>% 
  clean_names()

# Join with NMDRegHumanTxome_tbl -> *NOTE* joining by gene_name is probably not the best option, but sufficient for now
Fischer_2020_TableS1_joined <- NMDRHT.v1.2_MainTable %>% 
  left_join(Fischer_2020_TableS1 %>% 
              filter(upf1_cli_pseq_wt_vs_deaa_lee_et_al_2015 == "+") %>% 
              dplyr::select(gene_name) %>% 
              distinct(gene_name) %>% 
              mutate(UPF1_helicase_dependent = TRUE))

# Try to recover lost gene_names
Fischer_2020_TableS1_missingIDs <- Fischer_2020_TableS1 %>% 
  filter(upf1_cli_pseq_wt_vs_deaa_lee_et_al_2015 == "+") %>% 
  filter(!gene_name %in% Fischer_2020_TableS1_joined$gene_name) %>% 
  distinct(ncbi_refseq) %>% 
  separate(ncbi_refseq, into=c("ncbi_refseq", NA), sep="\\.")

library("biomaRt")

# Initialize bioMart | GENCODE v.42 is equal to ensembl version 108
ensembl_108 = useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 108)

# Perform bioMart request
Fischer_2020_TableS1_ENSG <- getBM(c("ensembl_gene_id_version","hgnc_symbol", "refseq_mrna"), "refseq_mrna", Fischer_2020_TableS1_missingIDs$ncbi_refseq, ensembl_108)

# Supplement recovered gene_names
# Found transcripts for in total 1438 genes (out of 1462 genes in TableS1)
Fischer_2020_TableS1_joined_fixed <- Fischer_2020_TableS1_joined %>% 
  left_join(Fischer_2020_TableS1_ENSG %>% 
              dplyr::rename("gene_id"="ensembl_gene_id_version") %>% 
              dplyr::distinct(gene_id) %>% 
              mutate(UPF1_helicase_dependent_fixed = TRUE)) %>% 
  mutate(UPF1_helicase_dependent = case_when(is.na(UPF1_helicase_dependent) & !is.na(UPF1_helicase_dependent_fixed) ~ UPF1_helicase_dependent_fixed,
                                             is.na(UPF1_helicase_dependent) ~ FALSE,
                                             TRUE ~ UPF1_helicase_dependent)) %>% 
  dplyr::select(-c(UPF1_helicase_dependent_fixed))

# Stats
Fischer_2020_TableS1_joined_fixed %>% 
  mutate(utr3_mfe_nt = -utr3_mfe_nt) %>% 
  group_by(UPF1_helicase_dependent) %>% 
  rstatix::get_summary_stats(utr3_mfe_nt)

# Reproduce higher structure in UPF1_helicase_dependent transcripts
Fischer_2020_TableS1_joined_fixed %>% 
  mutate(UPF1_helicase_dependent = case_when(is.na(UPF1_helicase_dependent) ~ FALSE,
                                             TRUE ~ UPF1_helicase_dependent)) %>% 
  ggplot(aes(x=-utr3_mfe_nt,
             color=UPF1_helicase_dependent)) +
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
  stat_ecdf(linewidth=0.5,
            geom="line") +
  # ggrastr::rasterise(stat_ecdf(),
  #                    dpi=1200, dev = "cairo") +
  scale_color_manual(values=c("TRUE" = "#008080",
                              "FALSE" = "#EDAE49")) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C5_helicase_dependent_ECDF.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


# More coding transcripts in UPF1_helicase_dependent tx
Fischer_2020_TableS1_joined_fixed %>% 
  filter(NMD_tx_status != "lncRNA") %>% 
  dplyr::count(UPF1_helicase_dependent,NMD_tx_status) %>% 
  group_by(UPF1_helicase_dependent) %>% 
  mutate(n_per = n/sum(n)) %>% 
  ggplot(aes(x=UPF1_helicase_dependent,
             y=n_per,
             fill=NMD_tx_status)) +
  geom_col()

# plot - Up early
Fischer_2020_TableS1_joined_fixed %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%  
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  group_by(condition_2, UPF1_helicase_dependent, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=UPF1_helicase_dependent,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="UPF1 helicase-dependent\nbound transcripts\n(Fischer et al. 2020)",
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(2*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C6_utr3_structure_UPF1_helicase_dependent_up_early_median_log2FC_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# plot - Global
Fischer_2020_TableS1_joined_fixed %>% 
  # filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%  
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  group_by(condition_2, UPF1_helicase_dependent, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=UPF1_helicase_dependent,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="UPF1 helicase-dependent\nbound transcripts\n(Fischer et al. 2020)",
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(2*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_C7_utr3_structure_UPF1_helicase_dependent_global_median_log2FC_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (D) RMD -----------------------------------------------------------
###

#### RMD - REG-1 RIP -----------------------------------------
# Supplmental table Table S1 - from PMID: 26000482
RMD_tableS1 <- read_excel("Resources/External/Mino2015_PMID_26000482/mmc2.xlsx", 
                          skip = 1) %>% 
  clean_names() %>% 
  mutate(significant = case_when(for_sorting < 0.001 ~ TRUE,
                                 TRUE ~ FALSE))

# Try to obtain GENCODE IDs from biomart
RMD_tableS1_filt <- RMD_tableS1 %>% 
  dplyr::filter(significant == TRUE) %>% 
  dplyr::select(acc_number, symbol) %>% 
  separate_rows(acc_number, sep = ",")

# Try to recover lost gene_names
RMD_tableS1_filt_missingIDs <- RMD_tableS1_filt %>% 
  dplyr::rename("ncbi_refseq" = "acc_number") %>% 
  distinct(ncbi_refseq)

library("biomaRt")

# Initialize bioMart | GENCODE v.42 is equal to ensembl version 108
ensembl_108 <- biomaRt::useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 108)

# Perform bioMart request
RMD_tableS1_filt_ENSG <- getBM(c("ensembl_gene_id_version","hgnc_symbol", "refseq_mrna"), "refseq_mrna", RMD_tableS1_filt_missingIDs$ncbi_refseq, ensembl_108)

RMD_tableS1_filt_ENSG_supp <- RMD_tableS1_filt_ENSG %>% 
  left_join(RMD_tableS1_filt %>% dplyr::rename("refseq_mrna" = "acc_number")) %>% 
  distinct(ensembl_gene_id_version, hgnc_symbol, .keep_all = TRUE) %>%
  #filter out duplicate gene IDs (for symbols IER3, LSM2 and HHAT)
  filter(!ensembl_gene_id_version %in% c("ENSG00000280680.3",
                                         "ENSG00000225998.6",
                                         "ENSG00000231502.6",
                                         "ENSG00000172850.11",
                                         "ENSG00000224979.6",
                                         "ENSG00000236826.7",
                                         "ENSG00000206478.5",
                                         "ENSG00000227231.3",
                                         "ENSG00000230128.3",
                                         "ENSG00000237155.3",
                                         "ENSG00000235030.3")) %>% 
  dplyr::rename("gene_id" = "ensembl_gene_id_version",
                "gene_name" = "hgnc_symbol")

# How many genes
RMD_tableS1_filt_ENSG_supp %>% 
  group_by(gene_id,gene_name) %>% 
  dplyr::count(gene_name) %>% 
  arrange(desc(n))

# How many NMDRHT transcripts
NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"
              ))) %>% 
  right_join(RMD_tableS1_filt_ENSG_supp) %>% 
  dplyr::count()

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"
              ))) %>% 
  right_join(RMD_tableS1_filt_ENSG_supp) %>% 
  mutate(DTE_cluster = case_when(is.na(DTE_cluster) ~ "not_expressed", 
                                 TRUE ~ DTE_cluster)) %>% 
  dplyr::count(DTE_cluster) %>% 
  filter(!is.na(DTE_cluster)) %>% 
  mutate(DTE_cluster = fct_relevel(DTE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "up 4:inverse",
                                   "expressed",
                                   "not_expressed",
                                   "down 4:inverse",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early")) %>% 
  ggplot(aes(x="REG-1 RIP",
             y=n,
             fill=DTE_cluster)) +
  geom_col(color="black",
           linewidth=0.1) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                             "up 2:delayed" = "#b2182b",
                             "up 3:late" = "#d6604d",
                             "up 4:inverse" = "#cfbeb4",
                             "down 4:inverse" = "#b9c3c8",
                             "down 3:late" = "#92c5de",
                             "down 2:delayed" = "#4393c3",
                             "down 1:early" = "#053061",
                             "expressed" = "gray30",
                             "not_expressed" = "gray5")) +
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
  geom_text(aes(label = case_when(n > 10 & DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n > 10 & !DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  labs(x="",
       y="number of transcripts",
       fill="DTE cluster") +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_D1_Regnase_RIP_DTE_cluster_absolute.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


###### Plot RMD RIP ------------------------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  right_join(RMD_tableS1_filt_ENSG_supp) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  group_by(condition_2, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=as_factor(NMD_50nt_rule),
             y=condition_2,
             fill=median_logFC)) +
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
        # legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  # facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=paste0("50-nt rule"),
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(3*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_D2_Regnase_RIP_median_log2FC_all_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### RMD - REG-1 HITS-CLIP -----------------------------------------
# Supplmental table Table S4 - from PMID: 26000482
RMD_tableS4 <- read_excel("Resources/External/Mino2015_PMID_26000482/mmc4.xlsx", 
                          sheet = "Intersect") %>%  
  clean_names() %>% 
  dplyr::rename("strand" = "x13")

# Convert hg19 to hg38 coordinates
library(liftOver)

# Convert data frame to GRanges
RMD_tableS4_GRanges <- makeGRangesFromDataFrame(RMD_tableS4,
                                                keep.extra.columns=TRUE)

# Import chain from UCSC (https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz)
hg19_to_hg38_ch = import.chain("/home/volker/reference/hg19ToHg38.over.chain")

# Make chr compatible
seqlevelsStyle(RMD_tableS4_GRanges) = "UCSC"

# Perform liftover and unlist
RMD_tableS4_GRanges_hg38 = liftOver(RMD_tableS4_GRanges, hg19_to_hg38_ch)
RMD_tableS4_GRanges_hg38_unlist = unlist(RMD_tableS4_GRanges_hg38)

# Generate tibble
RMD_tableS4_hg38 <- as_tibble(RMD_tableS4_GRanges_hg38_unlist)

# Prepare for intersect
RMD_tableS4_hg38_forMatch <- RMD_tableS4_hg38 %>% 
  dplyr::select(gene_name, seqnames, start, end, strand) %>% 
  dplyr::rename("chr" = "seqnames")

library(tidygenomics)

# read in GTF file
NMDRegHumanTxome.v1.2_GTF <- rtracklayer::readGFF("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf")

# Prepare NMDRHT GTF data for matching
NMDRegHumanTxome.v1.2_GTF_forMatch <- NMDRegHumanTxome.v1.2_GTF %>% 
  filter(type == "exon") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqid") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect
RMD_tableS4_hg38_NMDRHT_intersect <- genome_intersect(RMD_tableS4_hg38_forMatch, NMDRegHumanTxome.v1.2_GTF_forMatch, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
RMD_tableS4_hg38_NMDRHT_intersect_filt <- RMD_tableS4_hg38_NMDRHT_intersect %>% 
  filter(strand.x == strand.y) %>% 
  dplyr::select(-c(strand.y)) %>% 
  dplyr::rename("gene_name_RMD" = "gene_name",
                "strand" = "strand.x") %>% 
  left_join(NMDRHT.v1.2_MainTable)

# How many distinct transcript ids
RMD_tableS4_hg38_NMDRHT_intersect_filt %>% 
  distinct(transcript_id) %>% 
  dplyr::count()

### Check for binding to 5'UTR, CDS or 3'UTR --------------------------------

# Import GTF file
gtf <- rtracklayer::import(file.path("Resources/NMDRHT/NMDRHT.v1.2.sort.gtf"), format="gtf")

# Generate TxDb from annotation
NMDRHT_txdb1 <- makeTxDbFromGRanges(gtf, drop.stop.codons=FALSE)

# Get 5'UTRs from NMDRHT annotation and generate relevant data frame
NMDRHT_utr5 <- fiveUTRsByTranscript(NMDRHT_txdb1, use.names=TRUE)

# Get CDS from NMDRHT annotation and generate relevant data frame
NMDRHT_cds <- cdsBy(NMDRHT_txdb1, by="tx", use.names=TRUE)

# Get 3'UTRs from NMDRHT annotation and generate relevant data frame
NMDRHT_utr3 <- threeUTRsByTranscript(NMDRHT_txdb1, use.names=TRUE)


##### 5'UTR -------------------------------------------------------------------

NMDRHT_utr5_tbl <- as_tibble(NMDRHT_utr5) %>% 
  dplyr::rename("transcript_id" = "group_name") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqnames") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect
RMD_tableS4_hg38_NMDRHT_intersect_UTR5 <- genome_intersect(RMD_tableS4_hg38_forMatch, NMDRHT_utr5_tbl, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
RMD_tableS4_hg38_NMDRHT_intersect_UTR5 <- RMD_tableS4_hg38_NMDRHT_intersect_UTR5 %>% 
  filter(strand.x == strand.y) %>% 
  dplyr::select(-c(strand.y)) %>% 
  dplyr::rename("gene_name_RMD" = "gene_name",
                "strand" = "strand.x") %>% 
  left_join(NMDRHT.v1.2_MainTable)

# How many distinct transcript ids
RMD_tableS4_hg38_NMDRHT_intersect_UTR5 %>% 
  distinct(transcript_id) %>% 
  dplyr::count()

##### CDS -------------------------------------------------------------------

NMDRHT_cds_tbl <- as_tibble(NMDRHT_cds) %>% 
  dplyr::rename("transcript_id" = "group_name") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqnames") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect
RMD_tableS4_hg38_NMDRHT_intersect_cds <- genome_intersect(RMD_tableS4_hg38_forMatch, NMDRHT_cds_tbl, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
RMD_tableS4_hg38_NMDRHT_intersect_cds <- RMD_tableS4_hg38_NMDRHT_intersect_cds %>% 
  filter(strand.x == strand.y) %>% 
  dplyr::select(-c(strand.y)) %>% 
  dplyr::rename("gene_name_RMD" = "gene_name",
                "strand" = "strand.x") %>% 
  left_join(NMDRHT.v1.2_MainTable)

# How many distinct transcript ids
RMD_tableS4_hg38_NMDRHT_intersect_cds %>% 
  distinct(transcript_id) %>% 
  dplyr::count()

##### 3'UTR -------------------------------------------------------------------

NMDRHT_utr3_tbl <- as_tibble(NMDRHT_utr3) %>% 
  dplyr::rename("transcript_id" = "group_name") %>% 
  relocate(transcript_id) %>% 
  dplyr::rename("chr" = "seqnames") %>% 
  dplyr::select(transcript_id, chr, start, end, strand) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# Perform intersect
RMD_tableS4_hg38_NMDRHT_intersect_utr3 <- genome_intersect(RMD_tableS4_hg38_forMatch, NMDRHT_utr3_tbl, by=c("chr", "start", "end"), mode="both")

# Filter for same transcript_id and strand
# Prepare for output and join with ORFquant data
RMD_tableS4_hg38_NMDRHT_intersect_utr3 <- RMD_tableS4_hg38_NMDRHT_intersect_utr3 %>% 
  filter(strand.x == strand.y) %>% 
  dplyr::select(-c(strand.y)) %>% 
  dplyr::rename("gene_name_RMD" = "gene_name",
                "strand" = "strand.x") %>% 
  left_join(NMDRHT.v1.2_MainTable)

# How many distinct transcript ids
RMD_tableS4_hg38_NMDRHT_intersect_utr3 %>% 
  distinct(transcript_id) %>%
  dplyr::count()


### Plot DTE cluster --------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  right_join(RMD_tableS4_hg38_NMDRHT_intersect_UTR5 %>% 
               mutate(region="UTR5") %>% 
               bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_cds %>% 
                           mutate(region="CDS")) %>% 
               bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_utr3 %>% 
                           mutate(region="UTR3")) %>% 
               dplyr::select(transcript_id, region)) %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  mutate(DTE_cluster = case_when(is.na(DTE_cluster) ~ "not_expressed", 
                                 TRUE ~ DTE_cluster)) %>% 
  dplyr::count(DTE_cluster) %>% 
  filter(!is.na(DTE_cluster)) %>% 
  mutate(DTE_cluster = fct_relevel(DTE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "up 4:inverse",
                                   "expressed",
                                   "not_expressed",
                                   "down 4:inverse",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early")) %>% 
  ggplot(aes(x="REG1-CLIP",
             y=n,
             fill=DTE_cluster)) +
  geom_col(color="black",
           linewidth=0.1) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                             "up 2:delayed" = "#b2182b",
                             "up 3:late" = "#d6604d",
                             "up 4:inverse" = "#cfbeb4",
                             "down 4:inverse" = "#b9c3c8",
                             "down 3:late" = "#92c5de",
                             "down 2:delayed" = "#4393c3",
                             "down 1:early" = "#053061",
                             "expressed" = "gray30",
                             "not_expressed" = "gray5")) +
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
  geom_text(aes(label = case_when(n > 10 & DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n > 10 & !DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  labs(x="",
       y="number of transcripts",
       fill="DTE cluster") +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_D3_Regnase_CLIP_DTE_cluster_absolute.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Plot Region--------------------------------------------------------------------

RMD_tableS4_hg38_NMDRHT_intersect_UTR5 %>% 
  mutate(region="UTR5") %>% 
  bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_cds %>% 
              mutate(region="CDS")) %>% 
  bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_utr3 %>% 
              mutate(region="UTR3")) %>% 
  dplyr::count(region) %>% 
  mutate(region = fct_relevel(region,
                              "UTR5",
                              "CDS",
                              "UTR3")) %>% 
  ggplot(aes(x="region",
             y=n,
             fill=region)) +
  geom_col(color="black",
           linewidth=0.1) +
  scale_fill_manual(values=c("UTR5" = "gray90",
                             "CDS" = "#214D65",
                             "UTR3" = "#A9A196")) +
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
  geom_text(aes(label = case_when(n > 10 & region %in% c("CDS")  ~ paste(n),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n > 10 & !region %in% c("CDS")  ~ paste(n),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  labs(x="",
       y="number of binding sites",
       fill="Region") +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_D4_Regnase_CLIP_region.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Plot log2FC -------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  right_join(RMD_tableS4_hg38_NMDRHT_intersect_UTR5 %>% 
               mutate(region="UTR5") %>% 
               bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_cds %>% 
                           mutate(region="CDS")) %>% 
               bind_rows(RMD_tableS4_hg38_NMDRHT_intersect_utr3 %>% 
                           mutate(region="UTR3")) %>% 
               dplyr::select(transcript_id, region)) %>% 
  distinct(transcript_id, condition_2, .keep_all = TRUE) %>% 
  mutate(region = fct_relevel(region,
                              "UTR5",
                              "CDS",
                              "UTR3")) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(region)) %>% 
  group_by(condition_2, region, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=region,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=paste0("region"),
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(3*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_D5_Regnase_CLIP_median_log2FC_all_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (E) SMD -----------------------------------------------------------
###

# STAU1 CLIP targets ------------------------------------------------------
# Import data from Supplementary Table 2 |  PMID: 25799984
STAU1_CLIP_hits <- read_excel("Resources/External/Sugimoto2015_PMID_25799984/41586_2015_BFnature14280_MOESM50_ESM_mod.xlsx") %>% 
  dplyr::rename("gene_id_plain" = "seqnames")

# How many distinct genes are in STAU1 CLIP data? - n = 2913
STAU1_CLIP_hits %>% 
  distinct(gene_id_plain) %>%
  dplyr::count() 

NMDRHT.v1.2_MainTable %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  right_join(STAU1_CLIP_hits %>% 
               distinct(gene_id_plain, .keep_all = TRUE) %>% 
               mutate(STAU1_CLIP = "TRUE") %>% 
               dplyr::select(gene_id_plain, STAU1_CLIP)) %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::count()

### Plot DTE cluster --------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  right_join(STAU1_CLIP_hits %>% 
               distinct(gene_id_plain, .keep_all = TRUE) %>% 
               mutate(STAU1_CLIP = "TRUE") %>% 
               dplyr::select(gene_id_plain, STAU1_CLIP)) %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  mutate(DTE_cluster = case_when(is.na(DTE_cluster) ~ "not_expressed", 
                                 TRUE ~ DTE_cluster)) %>% 
  dplyr::count(DTE_cluster) %>% 
  filter(!is.na(DTE_cluster)) %>% 
  mutate(DTE_cluster = fct_relevel(DTE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "up 4:inverse",
                                   "expressed",
                                   "not_expressed",
                                   "down 4:inverse",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early")) %>% 
  ggplot(aes(x="STAU1-hiCLIP",
             y=n,
             fill=DTE_cluster)) +
  geom_col(color="black",
           linewidth=0.1) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                             "up 2:delayed" = "#b2182b",
                             "up 3:late" = "#d6604d",
                             "up 4:inverse" = "#cfbeb4",
                             "down 4:inverse" = "#b9c3c8",
                             "down 3:late" = "#92c5de",
                             "down 2:delayed" = "#4393c3",
                             "down 1:early" = "#053061",
                             "expressed" = "gray30",
                             "not_expressed" = "gray5")) +
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
  geom_text(aes(label = case_when(n > 10 & DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n > 10 & !DTE_cluster %in% c("up 1:early", "up 2:delayed","expressed", "not_expressed", "down 2:delayed", "down 1:early")  ~ paste(n),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  labs(x="",
       y="number of transcripts",
       fill="DTE cluster") +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_E1_STAU1_hiCLIP_DTE_cluster_absolute.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Plot log2FC -------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  left_join(STAU1_CLIP_hits %>% 
              distinct(gene_id_plain, .keep_all = TRUE) %>% 
              mutate(STAU1_CLIP = "TRUE") %>% 
              dplyr::select(gene_id_plain, STAU1_CLIP)) %>% 
  distinct(transcript_id, condition_2, .keep_all = TRUE) %>% 
  mutate(STAU1_CLIP = case_when(is.na(STAU1_CLIP) ~ "FALSE",
                                STAU1_CLIP == "TRUE" ~"TRUE")) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(STAU1_CLIP)) %>% 
  group_by(condition_2, STAU1_CLIP, NMD_50nt_rule) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=STAU1_CLIP,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x=paste0("STAU1 CLIP"),
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(3*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_E2_STAU1_hiCLIP_median_log2FC_all_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")



###
## Rev_1 - Figure 6 - (F) Individual targets -----------------------------------------------------------
###

# Import list
UPF1_non_NMD_targets <- read_excel("Resources/External/UPF1_non_NMD_targets.xlsx")

# Join with DTE -> log2FC MAXIMUM set to abs(5)!!!
UPF1_non_NMD_targets_forPlot <- UPF1_non_NMD_targets %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h",
                                        "UPF1_FKBP_HCT_0h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_0h",
                                        "UPF1_FKBP_HEK_12h")) %>% 
              bind_rows(qlf_df_NMDRHT_combined) %>% 
              dplyr::select(transcript_id, logFC, FDR, condition_2)) %>%
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "UPF1_Nter_0h",
                                           "UPF1_Nter_2h",
                                           "UPF1_Nter_4h",
                                           "UPF1_Nter_8h",
                                           "UPF1_Nter_12h",
                                           "UPF1_Nter_24h",
                                           "UPF1_Nter_48h",
                                           "UPF1_FKBP_HCT_0h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_FKBP_HEK_0h",
                                           "UPF1_FKBP_HEK_12h",
                                           "NMDRHT_ONT_dRNA_12h_IAA",
                                           "NMDRHT_PacBio_12h_IAA"))) %>% 
  filter(!is.na(condition_2)) %>% 
  mutate(gene_transcript_name = paste0(gene_name," | ",transcript_id)) %>% 
  mutate(logFC = replace(logFC, logFC > 5, 5)) %>% 
  mutate(logFC = replace(logFC, logFC < -5, -5)) 

UPF1_non_NMD_targets_forPlot %>% 
  rstatix::get_summary_stats(logFC)

UPF1_non_NMD_targets_forPlot %>% 
  mutate(FDR = -log10(FDR)) %>% 
  rstatix::get_summary_stats(FDR)

## long 3'UTR --------------------------------------------------------------

UPF1_non_NMD_targets_long_utr3 <- UPF1_non_NMD_targets_forPlot %>% 
  filter(pathway == "long_utr3")

UPF1_non_NMD_targets_long_utr3 %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::count(NMD_50nt_rule)


### 50 nt rule --------------------------------------------------------------

UPF1_non_NMD_targets_long_utr3 %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  ggplot(aes(x=gene_transcript_name,
             y=condition_2,
             fill=NMD_50nt_rule)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("TRUE" = "#E5BF86",
                             "FALSE" = "#214D65")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_long_utr3 %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F1_long_utr3_targets_50nt_rule.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Structural Category --------------------------------------------------------------

UPF1_non_NMD_targets_long_utr3 %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="SQANTI3",
             fill=structural_category_simple_2)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("FSM" = "#3F4889",
                             "NIC" = "#71CF57",
                             "NNIC" = "#BBDF27")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_long_utr3 %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F2_long_utr3_targets_SQANTI3.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### log2FC ------------------------------------------------------------------

UPF1_non_NMD_targets_long_utr3 %>% 
  ggplot(aes(x=gene_transcript_name,
             y=condition_2)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=logFC,
    size=-log10(FDR)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2),
             limits = c(0,25)) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  # theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "right")) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(13*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_long_utr3 %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F3_long_utr3_targets_log2FC.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

## SRD --------------------------------------------------------------

UPF1_non_NMD_targets_SRD <- UPF1_non_NMD_targets_forPlot %>% 
  filter(pathway %in% c("SRD_HSU", "SRD_PSU"))

UPF1_non_NMD_targets_SRD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::count(NMD_50nt_rule)


### 50 nt rule --------------------------------------------------------------

UPF1_non_NMD_targets_SRD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(pathway, gene_transcript_name) %>% 
  mutate(gene_transcript_name = fct_inorder(gene_transcript_name)) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="NMD_50nt_rule",
             fill=NMD_50nt_rule)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("TRUE" = "#E5BF86",
                             "FALSE" = "#214D65")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SRD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F4_SRD_targets_50nt_rule.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Structural Category --------------------------------------------------------------

UPF1_non_NMD_targets_SRD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(pathway, gene_transcript_name) %>% 
  mutate(gene_transcript_name = fct_inorder(gene_transcript_name)) %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="SQANTI3",
             fill=structural_category_simple_2)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("FSM" = "#3F4889",
                             "NIC" = "#71CF57",
                             "NNIC" = "#BBDF27")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SRD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F5_SRD_targets_SQANTI3.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### log2FC ------------------------------------------------------------------

UPF1_non_NMD_targets_SRD %>% 
  # distinct(transcript_id, .keep_all = TRUE) %>% 
  arrange(pathway, gene_transcript_name) %>% 
  mutate(gene_transcript_name = fct_inorder(gene_transcript_name)) %>% 
  ggplot(aes(x=gene_transcript_name,
             y=condition_2)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=logFC,
    size=-log10(FDR)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2),
             limits = c(0,25)) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  # theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "right")) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(13*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SRD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F6_SRD_targets_log2FC.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

## SMD --------------------------------------------------------------

UPF1_non_NMD_targets_SMD <- UPF1_non_NMD_targets_forPlot %>% 
  filter(pathway == "SMD")

UPF1_non_NMD_targets_SMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::count(NMD_50nt_rule)


### 50 nt rule --------------------------------------------------------------

UPF1_non_NMD_targets_SMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="50_nt_rule",
             fill=NMD_50nt_rule)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("TRUE" = "#E5BF86",
                             "FALSE" = "#214D65")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F7_SMD_targets_50nt_rule.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Structural Category --------------------------------------------------------------

UPF1_non_NMD_targets_SMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="SQANTI3",
             fill=structural_category_simple_2)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("FSM" = "#3F4889",
                             "NIC" = "#71CF57",
                             "NNIC" = "#BBDF27")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F8_SMD_targets_SQANTI3.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### log2FC ------------------------------------------------------------------

UPF1_non_NMD_targets_SMD %>% 
  ggplot(aes(x=gene_transcript_name,
             y=condition_2)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=logFC,
    size=-log10(FDR)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2),
             limits = c(0,25)) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  # theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "right")) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(13*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_SMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F9_SMD_targets_log2FC.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

## RMD --------------------------------------------------------------

UPF1_non_NMD_targets_RMD <- UPF1_non_NMD_targets_forPlot %>% 
  filter(pathway == "RMD")

UPF1_non_NMD_targets_RMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::count(NMD_50nt_rule)


### 50 nt rule --------------------------------------------------------------

UPF1_non_NMD_targets_RMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="50_nt_rule",
             fill=NMD_50nt_rule)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("TRUE" = "#E5BF86",
                             "FALSE" = "#214D65")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_RMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F10_RMD_targets_50nt_rule.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Structural Category --------------------------------------------------------------

UPF1_non_NMD_targets_RMD %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  ggplot(aes(x=gene_transcript_name,
             y="SQANTI3",
             fill=structural_category_simple_2)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  scale_fill_manual(values=c("FSM" = "#3F4889",
                             "NIC" = "#71CF57",
                             "NNIC" = "#BBDF27")) +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(1*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_RMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F11_RMD_targets_SQANTI3.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### log2FC ------------------------------------------------------------------

UPF1_non_NMD_targets_RMD %>% 
  ggplot(aes(x=gene_transcript_name,
             y=condition_2)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=logFC,
    size=-log10(FDR)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5,5),
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2),
             limits = c(0,25)) +
  theme_minimal() + 	
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
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  # theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "right")) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(13*2+0.6, "mm"),
                   cols = unit(length(UPF1_non_NMD_targets_RMD %>% distinct(transcript_id) %>% pull(transcript_id))*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_F12_RMD_targets_log2FC.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (G) Boruta Features -----------------------------------------------------------
###

load("Resources/Boruta/Boruta_data.rds")

# Plot Heatmap
combined_boruta_analysis %>% 
  # filter(decision %in% c("Confirmed", "Shadow")) %>%
  distinct(parameter, decision, median_imp, analysis) %>%
  mutate(confirmed = case_when(decision %in% c("Confirmed", "Shadow") ~ "",
                               decision == "Tentative" ~ "?",
                               decision == "Rejected" ~ "X")) %>% 
  complete(parameter, analysis) %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "Max. shadow feature",
                                 "Mean shadow feature",
                                 "Min. shadow feature",
                                 after = Inf)) %>% 
  mutate(analysis = fct_rev(fct_relevel(analysis,
                                        "all_FP",
                                        "all_RP",
                                        "all_MP",
                                        "noNMD_FP",
                                        "noNMD_RP",
                                        "noNMD_MP",
                                        "NMD_FP",
                                        "NMD_RP",
                                        "NMD_MP"))) %>% 
  ggplot(aes(x=parameter,
             y=analysis,
             fill=median_imp)) +
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
  geom_tile() +
  geom_text(aes(label=confirmed), color="white", size = 4*0.36, show_legend=FALSE) +
  scale_fill_viridis_c(option="plasma",na.value = "gray30") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1, size = 5, family="Arial")) +
  labs(x="",
       y="",
       fill="median\nimportance") +
  force_panelsizes(rows = unit(9*2, "mm"),
                   cols = unit(33*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_G1_boruta_combinedAnalysis_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Combined plot top parameters --------------------------------------------------
combined_boruta_analysis %>% 
  filter(analysis %in% c("all_RP",
                         "noNMD_RP",
                         "NMD_RP")) %>% 
  mutate(analysis = case_when(analysis == "all_RP" ~ "all",
                              analysis == "noNMD_RP" ~ "coding",
                              analysis == "NMD_RP" ~ "NMD")) %>% 
  arrange(desc(median_imp)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  mutate(parameter = fct_relevel(parameter,
                                 "(f) Max.",
                                 "(f) Mean",
                                 "(f) Min.",
                                 after = Inf)) %>% 
  mutate(analysis = (fct_relevel(analysis,
                                 "all",
                                 "coding",
                                 "NMD"))) %>% 
  filter(parameter %in% c("length 5'UTR",
                          "length CDS",
                          "length 3'UTR",
                          "50-nucleotide rule",
                          "num. downstream EJs",
                          "distance stop-lastEJ",
                          "NMD relevance (%)",
                          "(f) Max.",
                          "(f) Mean",
                          "(f) Min.")) %>% 
  mutate(parameter = fct_rev(fct_relevel(parameter,
                                         "50-nucleotide rule",
                                         "num. downstream EJs",
                                         "distance stop-lastEJ",
                                         "NMD relevance (%)",
                                         "length 5'UTR",
                                         "length CDS",
                                         "length 3'UTR",
                                         "(f) Max.",
                                         "(f) Mean",
                                         "(f) Min."))) %>% 
  ggplot(aes(y=parameter,
             x=importance,
             fill=analysis)) +
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
  stat_summary(fun.y=median, geom="point", shape=21, size=2, stroke=0.1, alpha=0.75) +
  scale_fill_manual(values=c("NMD" = "#E5BF86",
                             "coding" = "#214D65",
                             "all" = "#67001F")) +
  # labs(title="cluster 1:up | NMD reason=none | stabilized-only\nreduced parameter set | log2FC(DTE) ~ parameter") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(10*3+0.4, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_G2_Boruta_combined_reducedParameter_points.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Correlation -------------------------------------------------------------

NMDRHT.v1.2_MainTable_forCorrelation <- NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  filter(DTE_cluster == "up 1:early" & !is.na(NMD_50nt_rule)) %>% 
  mutate(NMD_50nt_rule = as.factor(NMD_50nt_rule))

NMDRHT.v1.2_MainTable_forCorrelation_bothNMDcoding <- NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  filter(DTE_cluster == "up 1:early" & !is.na(NMD_50nt_rule)) %>% 
  mutate(NMD_50nt_rule = "coding+NMD")

# Plot
NMDRHT.v1.2_MainTable_forCorrelation  %>% 
  # bind_rows(NMDRegHumanTxome.v1.2_tbl_GlobalTable_forCorrelation_bothNMDcoding) %>% 
  dplyr::select(c(logFC, condition_2, NMD_50nt_rule, 
                  num_of_downEJs,
                  stop_to_lastEJ,
                  NMD_n_sig_tx_perc,
                  utr5_len,
                  cds_len,
                  utr3_len)) %>% 
  group_by(condition_2,NMD_50nt_rule) %>% 
  rstatix::cor_test(vars = c("logFC"),
           vars2 = c("num_of_downEJs",
                     "stop_to_lastEJ",
                     "NMD_n_sig_tx_perc",
                     "utr5_len",
                     "cds_len",
                     "utr3_len")) %>% 
  mutate(var2 = case_when(var2 == "num_of_downEJs" ~ "num. downstream EJs",
                          var2 == "stop_to_lastEJ" ~ "distance stop-lastEJ",
                          var2 == "NMD_n_sig_tx_perc" ~ "NMD relevance (%)",
                          var2 == "utr5_len" ~ "length 5'UTR",
                          var2 == "cds_len" ~ "length CDS",
                          var2 == "utr3_len" ~ "length 3'UTR")) %>% 
  mutate(var2 = fct_rev(fct_relevel(var2,
                                    "num. downstream EJs",
                                    "distance stop-lastEJ",
                                    "NMD relevance (%)",
                                    "length 5'UTR",
                                    "length CDS",
                                    "length 3'UTR"))) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "UPF1_Nter_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(sig = case_when(p < 0.00001 ~ 21,
                         p >= 0.00001 ~ 23)) %>% 
  mutate(stroke= case_when(p < 0.00001 ~ 0.5,
                           p >= 0.00001 ~ 0.1)) %>% 
  ggplot(aes(x=cor,
             y=var2,
             fill=condition_2)) +
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
  geom_vline(xintercept = 0, color="gray20", linewidth=0.2) +
  geom_segment(aes(x=0, xend=cor, color=condition_2),
               position=position_dodge(width=0.75),
               alpha=0.75, show.legend=FALSE) +
  geom_point(aes(shape=sig,
                 stroke=stroke),
             alpha=0.75,
             position=position_dodge(width=0.75)) +
  scale_shape_identity(guide="none") +
  scale_color_manual(values=c("UPF1_Nter_12h" = "#00696D",
                              "UPF1_FKBP_HCT_12h" = "#6B56A7",
                              "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  facet_wrap(~NMD_50nt_rule) +
  labs(x="pearson correlation coefficient",
       y="") +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(30+0.4, "mm"),
                   cols = unit(15+0.4, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_G3_pearson_correlation_earlyUP_12h.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (H) Up-Early NMD relevance -----------------------------------------------------------
###

### heatmap ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_0h",
                                        "UPF1_Nter_2h",
                                        "UPF1_Nter_4h",
                                        "UPF1_Nter_8h",
                                        "UPF1_Nter_12h",
                                        "UPF1_Nter_24h",
                                        "UPF1_Nter_48h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(logFC = replace(logFC, logFC > 4, 4)) %>% 
  mutate(logFC = replace(logFC, logFC < -4, -4)) %>% 
  mutate(NMD_bin_tx = fct_rev(fct_relevel(NMD_bin_tx,
                                          "(75,100]",
                                          "(50,75]",
                                          "(25,50]",
                                          "[0,25]"
  ))) %>% 
  group_by(condition_2, NMD_bin_tx, NMD_50nt_rule, Mech_conclusion) %>% 
  dplyr::summarise(median_logFC = median(logFC)) %>% 
  ungroup() %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  # dplyr::count(utr3_mfe_nt_bin,DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=NMD_bin_tx,
             y=condition_2,
             fill=median_logFC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,4),
                       na.value = "grey90") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  # facet_wrap(~NMD_50nt_rule+Mech_conclusion, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="NMD relevance bin",
       y="",
       fill="median\nlog2FC") +
  force_panelsizes(rows = unit(7*2, "mm"),
                   cols = unit(4*2, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_H1_NMD_relevance_median_log2FC_up_early_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")




### Stat halfeye ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  ggplot(aes(y=NMD_bin,
             x=logFC,
             color=NMD_50nt_rule)) +
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
  ggdist::stat_halfeye(alpha = .5, 
               size=1,
               stroke=0,
               position = position_dodge(0.9),
               point_size = 2,
               show.legend = TRUE) +
  scale_color_manual(values=c("TRUE" = "#B09771",
                              "FALSE" = "#214D65")) +
  # scale_fill_manual(values=c("up 1:early" = "#67001f",
  #                             "up 2:delayed" = "#b2182b",
  #                             "up 3:late" = "#d6604d",
  #                             "down 3:late" = "#92c5de",
  #                             "down 2:delayed" = "#4393c3",
  #                             "down 1:early" = "#053061",
  #                             "expressed" = "gray30")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~condition_2, nrow=2) +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(30, "mm")) +
  guides(fill=guide_legend(reverse=T))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_H2_NMD_relevance_log2FC_up_early_distribution.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


#### Stats -------------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  group_by(NMD_bin, NMD_50nt_rule) %>% 
  rstatix::get_summary_stats(logFC)


NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early") %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>%  
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  group_by(NMD_50nt_rule) %>% 
  rstatix::dunn_test(logFC ~ NMD_bin) %>% 
  filter(group1 == "[0,25]" | group2 == "[0,25]")




### Mech - High NMD bin - up early - no 50nt rule ----------------------------------

# Load metabolism data
GSE84722_ST3 <- read_csv("Resources/External/Mukherjee2017_PMID_27870833/GSE84722_ST3.csv")

# Prepare for plot
NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == FALSE) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>% filter(condition == "UPF1_12h") %>% dplyr::select(transcript_id, L2FC_kdeg_tx)) %>% 
  separate(gene_id, c("gene_id_plain", NA), remove = FALSE) %>% 
  # left_join(GSE84722_ST3 %>% dplyr::select(gene_id_plain, Syn, Proc, Deg, CytNuc, PolyCyt, TrP)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  pivot_longer(cols=c(L2FC_kdeg, L2FC_ksyn, L2FC_kdeg_tx),
               names_to = "rate_name",
               values_to = "rate_value") %>% 
  mutate(rate_name = fct_relevel(rate_name,
                                 "L2FC_kdeg",
                                 "L2FC_ksyn", 
                                 "L2FC_kdeg_tx"))

NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  ggplot(aes(x=NMD_bin,
             y=rate_value,
             fill=rate_name)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(colour = 'gray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
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
  # geom_violin(alpha=0.5,
  #             linewidth=0.1,
  #             position=position_dodge(0.7), show.legend = FALSE) +
  # geom_boxplot(outliers=FALSE,
  #              position=position_dodge(0.7),
  #              width=0.15,
  #              linewidth=0.1,
  #              color="black") +
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  coord_cartesian(ylim = c(-3,3)) +
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d",
                             "L2FC_kdeg_tx" = "#cd7eb2")) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~rate_name, nrow=1) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(7*1.8, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_H3_no50ntRule_earlyUp_MetabolicRates_SLAMseq_boxplot.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Stats -------------------------------------------------------------------
NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin_tx,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  filter(rate_name %in% c("L2FC_kdeg_tx", "L2FC_kdeg")) %>% 
  group_by(NMD_bin, rate_name) %>% 
  dplyr::count(is.na(rate_value)) %>% 
  pivot_wider(names_from = `is.na(rate_value)`,
              values_from = n) %>% 
  janitor::clean_names() %>% 
  mutate(perc = 100*round(false/(true+false),2))

NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_forStats <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early" & NMD_50nt_rule == FALSE) %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>% filter(condition == "UPF1_12h") %>% dplyr::select(transcript_id, L2FC_kdeg_tx)) %>% 
  separate(gene_id, c("gene_id_plain", NA), remove = FALSE) %>% 
  # left_join(GSE84722_ST3 %>% dplyr::select(gene_id_plain, Syn, Proc, Deg, CytNuc, PolyCyt, TrP)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(NMD_50nt_rule))

# up 1:early as control - JUST AS TEST
NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_stats <- NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_forStats %>% rstatix::dunn_test(L2FC_kdeg_tx ~ NMD_bin_tx) %>% 
  filter(group1 == "[0,25]" | group2 == "[0,25]") %>% 
  bind_rows(NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_forStats %>% rstatix::dunn_test(L2FC_kdeg ~ NMD_bin_tx) %>% 
              filter(group1 == "[0,25]" | group2 == "[0,25]")) %>% 
  bind_rows(NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_forStats %>% rstatix::dunn_test(L2FC_ksyn ~ NMD_bin_tx) %>% 
              filter(group1 == "[0,25]" | group2 == "[0,25]"))

NMDRHT.v1.2_MainTable_earlyUpNo50nt_Mech_stats %>% 
  mutate(p.adj = replace(p.adj, p.adj < 1e-100, 1e-100)) %>% 
  mutate(NMD_bin = case_when(group1 == "[0,25]" ~ group2,
                             group1 != "[0,25]" ~ group1)) %>% 
  add_row(NMD_bin = "[0,25]", .y. = c("L2FC_kdeg_tx",
                                      "L2FC_kdeg",
                                      "L2FC_ksyn")) %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin,
                                       "(75,100]",
                                       "(50,75]",
                                       "(25,50]",
                                       "[0,25]"
  ))) %>% 
  ggplot(aes(x=NMD_bin,
             y=.y.,
             fill=-log10(p.adj))) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        strip.text = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_viridis_c(option = "magma",
                       begin=0.05,
                       end=0.95,
                       # limits = c(0, 100),
                       na.value = "grey80") +
  # theme(axis.text = element_blank(),
  #       axis.ticks = element_blank()) +
  labs(x="",
       y="",
       fill="-log10(padj)") + 
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(7*1.8, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_H4_MetabolicRates_earlyUp_NMDbin_DunnTest_padj.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 6 - (I) top 10 examples -----------------------------------------------------------
###

## High NMD bin - up early - no 50nt rule ----------------------------------
NMDRHT.v1.2_MainTable_earlyUpNo50nt_highNMDbin_top10_padj <- NMDRHT.v1.2_MainTable %>% 
  filter(DTE_cluster == "up 1:early" & NMD_bin_tx == "(75,100]" & NMD_50nt_rule == FALSE) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"
              ))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                    L2FC_kdeg,
                                                    padj_kdeg,
                                                    L2FC_ksyn,
                                                    padj_ksyn,
                                                    Mech_score,
                                                    kdeg_conclusion,
                                                    RNA_conclusion,
                                                    Mech_conclusion)) %>% 
  mutate(Mech_conclusion = case_when(is.na(Mech_conclusion) ~ "Not Sig.",
                                     TRUE ~ Mech_conclusion)) %>% 
  # filter(!is.na(utr3_mfe_nt)) %>% 
  filter(!is.na(condition_2)) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  filter(str_detect(ORFquant_NMD_tx_status_source, "^ORFquant")) %>% 
  arrange(FDR) %>% 
  distinct(transcript_id, .keep_all = TRUE) %>% 
  slice_min(FDR, n=10, with_ties = FALSE)


### Lollipop FDR&log2FC -----------------------------------------------------

NMDRHT.v1.2_MainTable_earlyUpNo50nt_highNMDbin_top10_padj %>% 
  dplyr::select(transcript_name,logFC,FDR, structural_category_simple, NMD_tx_source, NMD_n_sig_tx, NMD_n_sig_tx_perc, Mech_conclusion) %>% 
  mutate(transcript_name = fct_rev(fct_inorder(transcript_name))) %>% 
  mutate(FDR_log = -log10(FDR)) %>% 
  ggplot(aes(x=FDR_log,
             y=transcript_name,
             fill=logFC)) +
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
  geom_vline(xintercept = 0, color="gray20", linewidth=0.2) +
  geom_segment(aes(x=0, xend=FDR_log),
               color="gray20",
               alpha=0.75, show.legend=FALSE) +
  geom_point(color="black",
             size=2.25,
             shape=21,
             stroke=0.2) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,5.5),
                       na.value = "grey90") +
  labs(x="-log10(FDR)",
       y="") +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(25+0.48, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_I1_12h_top10_FDR.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Heatmap PARE log2FC -----------------------------------------------------

PARE_edgeR_combined <- read_csv("Resources/PARE/PARE_edgeR_combined.csv")

NMDRHT.v1.2_MainTable_earlyUpNo50nt_highNMDbin_top10_padj %>% 
  dplyr::select(transcript_name,logFC,FDR, structural_category_simple, NMD_tx_source, NMD_n_sig_tx, NMD_n_sig_tx_perc, Mech_conclusion) %>% 
  left_join(PARE_edgeR_combined %>% dplyr::select(transcript_name, L2FC_XS6, L2FC_XU1)) %>% 
  pivot_longer(cols=c(L2FC_XS6, L2FC_XU1),
               names_to="condition_2",
               values_to="log2FC") %>% 
  mutate(transcript_name = fct_rev(fct_inorder(transcript_name))) %>% 
  mutate(FDR_log = -log10(FDR)) %>% 
  ggplot(aes(x=condition_2,
             y=transcript_name,
             fill=log2FC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-1,5.5),
                       na.value = "grey90") +
  labs(x="",
       y="") +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(25+0.48, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_I2_12h_top10_PARE_log2FC.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Heatmap SLAM-seq EZbakR log2FC -----------------------------------------------------

NMDRHT.v1.2_MainTable_earlyUpNo50nt_highNMDbin_top10_padj %>% 
  dplyr::select(transcript_name,logFC,FDR, structural_category_simple, NMD_tx_source, NMD_n_sig_tx, NMD_n_sig_tx_perc, Mech_conclusion) %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>% filter(condition == "UPF1_12h") %>% dplyr::select(transcript_name, L2FC_kdeg_tx, L2FC_kdeg)) %>% 
  pivot_longer(cols=c(L2FC_kdeg_tx, L2FC_kdeg),
               names_to="rate",
               values_to="log2FC") %>% 
  mutate(transcript_name = fct_rev(fct_inorder(transcript_name))) %>% 
  mutate(FDR_log = -log10(FDR)) %>% 
  ggplot(aes(x=rate,
             y=transcript_name,
             fill=log2FC)) +
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       limits = c(-5.5,5.5),
                       na.value = "grey90") +
  labs(x="",
       y="") +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(25+0.48, "mm"),
                   cols = unit(5, "mm"))

ggsave(file.path("Plots", "Figure6", "Rev_1_F6_I2_12h_top10_kdeg_log2FC.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
