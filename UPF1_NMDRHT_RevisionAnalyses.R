#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_RevisionAnalyses.R
# Objective: Code for generating Revision Analyses for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_NMDRHT_Analysis"

###
# Correlation differential expression -----------------------------------------------------------
###

###
## Gene-level --------------------------------------------------------------
###

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")


### Correlation UPF1 24h depletion ------------------------------------------

DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  filter(condition_2 %in% c("UPF1_Nter_24h",
                            "UPF1_Nter_24h_R0h")) %>% 
  dplyr::select(gene_id, log2FoldChange, condition_2) %>% 
  pivot_wider(values_from = log2FoldChange,
              names_from = condition_2,
              values_fill = NA) %>% 
  ggplot(aes(x=UPF1_Nter_24h,
             y=UPF1_Nter_24h_R0h)) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  scale_color_viridis() +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                          color="darkred") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "p.value.label")),
                        size = 5*0.30,
                        color="gray20",
                        label.y = "bottom", label.x = "right") +
  theme(strip.background = element_rect(color="black", fill = "white", linewidth=0.1),
        strip.text = element_text(size=6),
        plot.background = element_blank(), 
        panel.background = element_rect(color="black", fill = "white", linewidth=0.1),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-15,15)) +
  scale_y_continuous(limits = c(-15,15)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  labs(title= "Correlation log2FC\nN-AID-UPF1 24h IAA versus\nN-AID-UPF1 24h IAA | 0h Recovery",
       subtitle= "Missing values filled with NA",
       caption = "Data from DESeq2_DGE_combined.csv",
       x="log2FC(N-AID-UPF1 24h IAA)",
       y="log2FC(N-AID-UPF1 24h IAA | 0h Recovery)",
       color="density") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("Plots", "Revision", "R2_DGE_Correlation_24h.pdf"),
       width = 15,
       height = 15,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Correlation UPF1 24h depletion 0-4h Recovery ------------------------------------------

DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_recovery_this_Study")) %>% 
  filter(condition_2 %in% c("UPF1_Nter_24h_R0h",
                            "UPF1_Nter_24h_R2h",
                            "UPF1_Nter_24h_R4h",
                            "UPF1_Nter_24h_R8h",
                            "UPF1_Nter_24h_R12h",
                            "UPF1_Nter_24h_R24h",
                            "UPF1_Nter_24h_R48h")) %>% 
  dplyr::select(gene_id, log2FoldChange, condition_2) %>% 
  pivot_wider(values_from = log2FoldChange,
              names_from = condition_2,
              values_fill = NA) %>% 
  column_to_rownames(var="gene_id") %>% 
  gather(-"UPF1_Nter_24h_R0h", key = "var", value = "value") %>% 
  mutate(var = fct_relevel(var,
                           "UPF1_Nter_24h_R2h",
                           "UPF1_Nter_24h_R4h",
                           "UPF1_Nter_24h_R8h",
                           "UPF1_Nter_24h_R12h",
                           "UPF1_Nter_24h_R24h",
                           "UPF1_Nter_24h_R48h")) %>% 
  ggplot(aes(x=value,
             y=UPF1_Nter_24h_R0h)) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  scale_color_viridis() +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                          color="darkred") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "p.value.label")),
                        size = 5*0.30,
                        color="gray20",
                        label.y = "bottom", label.x = "right") +
  theme(strip.background = element_rect(color="black", fill = "white", linewidth=0.1),
        strip.text = element_text(size=6),
        plot.background = element_blank(), 
        panel.background = element_rect(color="black", fill = "white", linewidth=0.1),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-15,15)) +
  scale_y_continuous(limits = c(-15,15)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  labs(title= "Correlation log2FC\nN-AID-UPF1 24h IAA | 0h Recovery versus\nN-AID-UPF1 24h IAA | longer Recovery",
       subtitle= "Missing values filled with NA",
       caption = "Data from DESeq2_DGE_combined.csv",
       x="log2FC(N-AID-UPF1 24h IAA | 0h Recovery)",
       y="log2FC(Contrast)",
       color="density") +
  facet_wrap(~ var) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("Plots", "Revision", "R2_DGE_Correlation_Recovery.pdf"),
       width = 15,
       height = 15,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Correlation UPF1 12h FKBP depletion ------------------------------------------

DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h")) %>% 
  dplyr::select(gene_id, log2FoldChange, condition_2) %>% 
  pivot_wider(values_from = log2FoldChange,
              names_from = condition_2,
              values_fill = NA) %>% 
  column_to_rownames(var="gene_id") %>% 
  gather(-"UPF1_Nter_12h", key = "var", value = "value") %>% 
  mutate(var = fct_relevel(var,
                           "UPF1_FKBP_HCT_12h",
                           "UPF1_FKBP_HEK_12h")) %>% 
  ggplot(aes(x=value,
             y=UPF1_Nter_12h)) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  scale_color_viridis() +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                          color="darkred") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "p.value.label")),
                        size = 5*0.30,
                        color="gray20",
                        label.y = "bottom", label.x = "right") +
  theme(strip.background = element_rect(color="black", fill = "white", linewidth=0.1),
        strip.text = element_text(size=6),
        plot.background = element_blank(), 
        panel.background = element_rect(color="black", fill = "white", linewidth=0.1),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-15,15)) +
  scale_y_continuous(limits = c(-15,15)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  labs(title= "Correlation log2FC\nN-AID-UPF1 12h IAA versus\nFKBP-UPF1 12h IAA",
       subtitle= "Missing values filled with NA",
       caption = "Data from DESeq2_DGE_combined.csv",
       x="log2FC(N-AID-UPF1 12h IAA)",
       y="log2FC(Contrast)",
       color="density") +
  facet_wrap(~ var) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("Plots", "Revision", "R2_DGE_Correlation_12h_FKBP.pdf"),
       width = 15,
       height = 15,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Correlation DGE overall --------------------------------------------

DESeq2_DGE_combined_ForCor <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  dplyr::select(gene_id, log2FoldChange, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = log2FoldChange) %>% 
  column_to_rownames(var="gene_id")

DESeq2_DGE_combined_Cor = cor(DESeq2_DGE_combined_ForCor, use="complete.obs")

ggcorrplot::ggcorrplot(DESeq2_DGE_combined_Cor, 
                       type = "full",
                       method = "square",
                       lab = FALSE,
                       tl.cex = 6,
                       lab_size = 2.5,
                       pch.cex = 2.5,
                       ggtheme = ggplot2::theme_minimal,
                       colors = c("#6D9EC1", "white", "#E46726"),
                       outline.col = "white")

ggsave(file.path("Plots", "Revision", "R2_DGE_Correlation_complete.pdf"),
       width = 15,
       height = 15,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
# Transcription Factors - protein-level -----------------------------------------------------------
###

# Import CSV file - downloaded on 12.11.2024 from: https://humantfs.ccbr.utoronto.ca/download.php
TF_DatabaseExtract_v_1_01 <- read_csv("Resources/External/TranscriptionFactors/DatabaseExtract_v_1.01.csv")

# Select relevant information
TF_DatabaseExtract_v_1_01_tbl <- TF_DatabaseExtract_v_1_01 %>% 
  janitor::clean_names() %>% 
  dplyr::rename("gene_id" = "ensembl_id")

# load protein group results
load("Resources/Proteomics/WCP_PG_data_20250127_GeneTxORF_long.rds")

# Join Proteomics with TranscriptFactor data - simplify information
WCP_PG_data_20250127_GeneTxORF_long_tf <- WCP_PG_data_20250127_GeneTxORF_long %>% 
  separate(gene_id, c("gene_id", NA)) %>% 
  dplyr::select(gene_id, gene_name, comparison, log2FC, padj) %>% 
  left_join(TF_DatabaseExtract_v_1_01_tbl %>% dplyr::select(gene_id, dbd, is_tf)) %>% 
  filter(is_tf == "Yes") %>% 
  mutate(significant = case_when(padj < 0.001 & abs(log2FC) > 0.5 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  arrange(padj) %>% 
  distinct()

# Determine significantly expressed TFs
WCP_PG_data_20250127_GeneTxORF_long_tf_sig_id <- WCP_PG_data_20250127_GeneTxORF_long_tf %>% 
  filter(significant == TRUE) %>%  
  distinct(gene_id) %>% 
  pull(gene_id)

# Plot Heatmap
WCP_PG_data_20250127_GeneTxORF_long_tf %>% 
  filter(gene_id %in% WCP_PG_data_20250127_GeneTxORF_long_tf_sig_id) %>% 
  arrange(desc(log2FC)) %>% 
  mutate(gene_name = (fct_inorder(gene_name))) %>% 
  ggplot(aes(y=comparison,
             x=gene_name,
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
  labs(x="transcription factors",
       y="") +
  force_panelsizes(rows = unit(5*2, "mm"),
                   cols = unit(104*1.8, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

ggsave(file.path("Plots", "Revision", "R3_TF_Proteomics_Heatmap.pdf"),
       width = 25,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

