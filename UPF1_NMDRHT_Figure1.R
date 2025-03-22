#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Figure1
# Objective: Code for generating Panels of Figure 1 + Figure S1 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_NMDRHT_Analysis"

###
## Rev_1 - Figure 1 - (A1) DGE -----------------------------------------------------------
###

# Comment: Plot A1 = DGE of shallow test-RNA-Seq of N-AID-UPF1 and C-AID-UPF1

### DGE - Prepare for Plot -----------------------------------------------------------

Rev_1_F1_A1_DESeq2_DGE_combined_forPlot <- DESeq2_DGE_combined %>% 
  filter(experimentSet == "HCT116_UPF1_AID_TestSeq_this_Study") %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(significant = case_when(padj < 0.0001 & abs(log2FoldChange) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(log2FoldChange>0 ~ "up",
                             log2FoldChange<0 ~ "down")) %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_gene = length(unique(gene_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, type) %>% 
  mutate(detected_gene_type = length(unique(gene_id))) %>% 
  ungroup()

#### Rev_1_F1_A1_Heatmap -----------------------------------------------------------

Rev_1_F1_A1_DESeq2_DGE_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, type, up_down, detected_gene_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(n_perTotalType = n/detected_gene_type) %>% 
  mutate(type_up_down = paste0(type, "-", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding-down",
                                     "coding-up",
                                     "lncRNA-down",
                                     "lncRNA-up",
                                     "other-down",
                                     "other-up"))) %>% 
  mutate(n = case_when(up_down == "down" ~ -n,
                       up_down == "up" ~ n)) %>% 
  mutate(n_perTotalType = case_when(up_down == "down" ~ -n_perTotalType,
                                    up_down == "up" ~ n_perTotalType)) %>% 
  ggplot(aes(x=type_up_down,
             y=condition_2,
             fill=n_perTotalType)) +
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
  geom_tile(color="black",
            linewidth=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey90") +
  force_panelsizes(rows = unit(length(unique(Rev_1_F1_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(12+0.6, "mm"))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y="",
       x="",
       fill="Fraction of sig.\nregulated genes\nper gene biotype") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.vjust = 0.5))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_A1_DESeq2_DGE_Heatmap_perType.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### Rev_1_F1_A1_Barplot -----------------------------------------------------------

Rev_1_F1_A1_DESeq2_DGE_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, type, up_down, detected_gene_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(type_up_down = paste0(type, "-", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding-up",
                                     "coding-down",
                                     "lncRNA-up",
                                     "lncRNA-down",
                                     "other-up",
                                     "other-down"))) %>% 
  mutate(n = case_when(up_down == "down" ~ -n,
                       up_down == "up" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n,
             fill=type)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6, color="black"),
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
  geom_vline(xintercept = 0,
             color="#963B5A",
             linewidth=0.25) +
  geom_col(color="black",
           linewidth = 0.1) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "other" = "#D9DADA")) +
  # scale_fill_viridis_d(option="G",
  #                      begin=0.25,
  #                      end=0.75) +
  labs(y="",
       x="Number of \nsig. regulated genes",
       fill="Gene biotype",
       size="Sig. DGE events") +
  force_panelsizes(rows = unit(length(unique(Rev_1_F1_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(20, "mm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE,
                           title.position = "top"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_A1_DESeq2_DGE_Barplot_perType.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### Rev_1_F1_A1_Examples -----------------------------------------------------------

# Define representative NMD-targeted genes plus UPF1
standard_heatmap_genes <- c("ZFAS1",
                            "SNHG12",
                            "GAS5",
                            "GADD45B",
                            "UPF1")

Rev_1_F1_A1_DESeq2_DGE_combined_forPlot %>%  
  filter(gene_name %in% standard_heatmap_genes) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, gene_name) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "UPF1",
                                 "ZFAS1",
                                 "SNHG12",
                                 "GAS5",
                                 "GADD45B")) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=gene_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2)) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="DGE log2FC") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                barwidth = unit(8, "mm"),
                                barheight = unit(2, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "bottom")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(unique(Rev_1_F1_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(5*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_A1_DESeq2_DGE_Examples.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 1 - (B) Essentiality -----------------------------------------------------------
###

CRISPR_DepMap_23Q2_Public_Score_Chronos <- read_csv("Resources/DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos).csv")

# All cell lines
CRISPR_DepMap_all <- CRISPR_DepMap_23Q2_Public_Score_Chronos %>% 
  mutate(cellLine = paste0(cell_line_display_name,"-",depmap_id)) %>% 
  dplyr::select(-c(starts_with("lineage_"), cell_line_display_name, depmap_id)) %>% 
  column_to_rownames(var = "cellLine") %>% 
  rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  pivot_wider(names_from = rowname) %>%
  dplyr::rename("gene" = 1) %>%
  mutate(MeanChronosScore = rowMeans(dplyr::select(., contains("ACH")),
                                     na.rm = TRUE)) %>% 
  as.data.frame() 

CRISPR_DepMap_all %>% 
  write_csv("Resources/DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos)_format.csv")

CRISPR_DepMap_all <- read_csv("Resources/DepMap/CRISPR_(DepMap_Public_23Q2+Score,_Chronos)_format.csv")

# Plot
CRISPR_DepMap_all_NMD <- CRISPR_DepMap_all %>% 
  dplyr::select(-MeanChronosScore) %>% 
  filter(gene %in% c("EIF4A3",
                     "MAGOH",
                     "MAGOHB",
                     "RBM8A",
                     "CASC3",
                     "SMG1",
                     "SMG5",
                     "SMG6",
                     "SMG7",
                     "SMG8",
                     "SMG9",
                     "UPF1",
                     "UPF2",
                     "UPF3A",
                     "UPF3B",
                     "STAU1",
                     "STAU2"
  )) %>% 
  pivot_longer(cols = -c("gene"),
               names_to = "cellLine",
               values_to = "score")

CRISPR_DepMap_all_NMD %>% 
  group_by(gene) %>% 
  mutate(meanScore = mean(score,na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(gene=fct_rev(fct_relevel(gene,
                                  "UPF1",
                                  "UPF2",
                                  "UPF3A",
                                  "UPF3B",
                                  "SMG1",
                                  "SMG5",
                                  "SMG6",
                                  "SMG7",
                                  "SMG8",
                                  "SMG9",
                                  "EIF4A3",
                                  "MAGOH",
                                  "MAGOHB",
                                  "RBM8A",
                                  "CASC3",
                                  "STAU1",
                                  "STAU2"))) %>% 
  ggplot(aes(x=score,
             y=gene,
             fill=meanScore)) +
  theme(legend.position="right", 
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
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  geom_vline(xintercept = -1,
             color="darkred",
             linewidth=0.25) +
  geom_boxplot(outlier.shape = NA,
               linewidth=0.2) +
  scale_fill_gradient(high="white",
                      low="#8B0000") +
  labs(x="CRISPR (DepMap Public 23Q2+Score, Chronos)",
       y="Gene",
       fill="mean\nChronos\nscore") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(30, "mm")) 

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_B1_CRISPR_Chronos_NMD.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### RVIS --------------------------------------------------------------------

genic_intolerance_2023_11_17 <- read_csv("Resources/DepMap/genic_intolerance_2023-11-17.csv") 

gnomad_v4_0_constraint_metrics <- read_table("Resources/DepMap/gnomad.v4.0.constraint_metrics.tsv")

gnomad_v4_0_constraint_metrics_NMD <- gnomad_v4_0_constraint_metrics %>% 
  filter(mane_select == TRUE) %>% 
  filter(str_detect(transcript,"ENST")) %>% 
  filter(gene %in% c("EIF4A3",
                     "MAGOH",
                     "MAGOHB",
                     "RBM8A",
                     "CASC3",
                     "SMG1",
                     "SMG5",
                     "SMG6",
                     "SMG7",
                     "SMG8",
                     "SMG9",
                     "UPF1",
                     "UPF2",
                     "UPF3A",
                     "UPF3B",
                     "STAU1",
                     "STAU2"
  )) %>% 
  add_row(gene = "UPF3B") %>% 
  mutate(gene_id=fct_rev(fct_relevel(gene,
                                     "UPF1",
                                     "UPF2",
                                     "UPF3A",
                                     "UPF3B",
                                     "SMG1",
                                     "SMG5",
                                     "SMG6",
                                     "SMG7",
                                     "SMG8",
                                     "SMG9",
                                     "EIF4A3",
                                     "MAGOH",
                                     "MAGOHB",
                                     "RBM8A",
                                     "CASC3",
                                     "STAU1",
                                     "STAU2")))


genic_intolerance_2023_11_17 %>% 
  mutate(RVIS_v4 = as.numeric(str_split_i(`ExAC_v2_0.05%popn`, " ", 1)),
         RVIS_v4_per = as.numeric(parse_number(str_split_i(`ExAC_v2_0.05%popn`, " ", 2))),
         ExAC_LoF_FDR = `LoF-FDR[ExAC]`) %>% 
  dplyr::select(GENE,RVIS_v4,RVIS_v4_per, ExAC_LoF_FDR) %>% 
  dplyr::rename("gene" = "GENE") %>% 
  left_join(gnomad_v4_0_constraint_metrics_NMD) %>% 
  dplyr::select(gene,RVIS_v4,ExAC_LoF_FDR,lof.oe,lof.oe_ci.lower,lof.oe_ci.upper, lof.pLI, lof.z_score) %>% 
  mutate(gene=fct_rev(fct_relevel(gene,
                                  "UPF1",
                                  "UPF2",
                                  "UPF3A",
                                  "UPF3B",
                                  "SMG1",
                                  "SMG5",
                                  "SMG6",
                                  "SMG7",
                                  "SMG8",
                                  "SMG9",
                                  "EIF4A3",
                                  "MAGOH",
                                  "MAGOHB",
                                  "RBM8A",
                                  "CASC3",
                                  "STAU1",
                                  "STAU2"))) %>% 
  ggplot(aes(x=lof.oe,
             y=gene,
             fill=lof.pLI)) +
  theme(legend.position="right", 
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
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_linerange(aes(xmin=lof.oe_ci.lower,
                     xmax=lof.oe_ci.upper)) +
  geom_point(aes(size=lof.pLI),
             shape=21) +
  # scale_fill_viridis_c(begin = 0.2,
  #                      end = 0.8,
  #                      direction = 1) +
  scale_fill_gradient(name = "Probability of\nloss-of-function\nintolerance",
                      labels = c(0,0.25,0.5,0.75,1),
                      guide = "legend",
                      limits = c(0,1),
                      low="white",
                      high="#8B0000") +
  scale_size(name = "Probability of\nloss-of-function\nintolerance",
             labels = c(0,0.25,0.5,0.75,1),
             limits = c(0,1),
             range = c(1,3)) +
  labs(x="loss-of-function observed/expected",
       y="Gene") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(30, "mm")) 

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_B2_LOF_RVIS_NMD.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


###
## Rev_1 - Figure 1 - (C) External Proteomics -----------------------------------------------------------
###

### Zecha 2018 --------------------------------------------------------------

Protein_half_life_HeLa_29414762 <- read_excel("Resources/External/Zecha2018_PMID_29414762/NMD_half-life_HeLa_34915_1_supp_65961_p3f7xq.xlsx", 
                                              sheet = "Properties and functions") %>% 
  janitor::clean_names()

NMD_names <- c("UPF1",
               "UPF2",
               "UPF3A",
               "UPF3B",
               "SMG1",
               "SMG5",
               "SMG6",
               "SMG7",
               "SMG8",
               "SMG9")

EJC_names <- c("EIF4A3",
               "MAGOH;MAGOHB",
               "MAGOH",
               "MAGOHB",
               "RBM8A",
               "CASC3")

STAU_names <- c("STAU1",
                "STAU2")

Names_combined <- c(NMD_names,
                    EJC_names,
                    STAU_names)

# Plot
Protein_half_life_HeLa_29414762 %>% 
  filter(gene_name_s %in% Names_combined) %>% 
  mutate(t1_2_h = as.numeric(t1_2_h)) %>% 
  mutate(copies_per_cell = as.numeric(copies_per_cell)) %>% 
  mutate(class = case_when(gene_name_s %in% NMD_names ~ "NMD",
                           gene_name_s %in% EJC_names ~ "EJC",
                           gene_name_s %in% STAU_names ~ "STAU",
                           TRUE ~ "other")) %>% 
  mutate(gene_name_s=fct_rev(fct_relevel(gene_name_s,
                                         "UPF1",
                                         "UPF2",
                                         "UPF3A",
                                         "UPF3B",
                                         "SMG1",
                                         "SMG5",
                                         "SMG6",
                                         "SMG7",
                                         "SMG8",
                                         "SMG9",
                                         "EIF4A3",
                                         "MAGOH;MAGOHB",
                                         "RBM8A",
                                         "CASC3",
                                         "STAU1",
                                         "STAU2"))) %>% 
  mutate(class = fct_relevel(class,
                             "NMD",
                             "EJC",
                             "STAU")) %>% 
  ggplot(aes(x=t1_2_h,
             y=gene_name_s,
             fill=class,
             size=log10(copies_per_cell))) +
  theme(legend.position="right", 
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
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_point(shape=21,
             color="black") +
  geom_text(aes(label = paste0("T1/2: ", round(t1_2_h,0), " h"),
                hjust = case_when(gene_name_s == "UPF1" ~ 1.25,
                                  TRUE ~ -0.25)),
            size = 4*0.36) +
  coord_cartesian(clip = 'off') +
  # scale_x_continuous(expand = expansion(c(0.05, 0.4))) +
  scale_size(range = c(1,3)) +
  scale_fill_manual(values=c("NMD" = "#CBAC7C",
                             "EJC" = "#50936D",
                             "STAU" = "#CFAE9F",
                             "other" = "darkgray")) +
  labs(x="protein half-life\nT1/2 [h]",
       y="",
       fill="Protein class",
       size="log10\n(Protein copies\nper cell)") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(20, "mm")) 

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_C1_NMD_EJC_STAU_protein_HL_2024.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Wang 2019 --------------------------------------------------------------
# Protein abundance - human tissues - Import from PMID: 30777892 - Table EV5

PMID_30777892_Table_EV5 <- read_excel("Resources/External/Wang2019_PMID_30777892/PMID_30777892_Table_EV5.xlsx", 
                                      sheet = "A. Protein copies") %>% 
  janitor::clean_names()


PMID_30777892_Table_EV5 %>%   
  filter(gene_name %in% Names_combined) %>% 
  mutate(class = case_when(gene_name %in% NMD_names ~ "NMD",
                           gene_name %in% EJC_names ~ "EJC",
                           gene_name %in% STAU_names ~ "STAU",
                           TRUE ~ "other")) %>% 
  mutate(gene_name=fct_rev(fct_relevel(gene_name,
                                       "UPF1",
                                       "UPF2",
                                       "UPF3A",
                                       "UPF3B",
                                       "SMG1",
                                       "SMG5",
                                       "SMG6",
                                       "SMG7",
                                       "SMG8",
                                       "SMG9",
                                       "EIF4A3",
                                       "MAGOH",
                                       "MAGOHB",
                                       "RBM8A",
                                       "CASC3",
                                       "STAU1",
                                       "STAU2"))) %>% 
  mutate(class = fct_relevel(class,
                             "NMD",
                             "EJC",
                             "STAU")) %>% 
  pivot_longer(cols=-c(protein_id,
                       gene_id,
                       gene_name,
                       class),
               values_to = "copies",
               names_to = "tissue") %>% 
  group_by(gene_name, class) %>% 
  summarize(mean_copies = mean(copies, na.rm = TRUE),
            SD_copies = sd(copies, na.rm = TRUE)) %>% 
  ggplot(aes(x=log10(mean_copies),
             y=gene_name,
             fill=class)) +
  theme(legend.position="right", 
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
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_point(shape=21,
             color="black") +
  geom_text(aes(label = paste0(round(mean_copies,0))),
            hjust = -0.3,
            size = 4*0.36) +
  xlim(2,8.75) +
  scale_fill_manual(values=c("NMD" = "#CBAC7C",
                             "EJC" = "#50936D",
                             "STAU" = "#CFAE9F",
                             "other" = "darkgray")) +
  labs(x=" mean\nlog10(Protein copies per cell)\nfrom 29 human tissues",
       y="",
       fill="Protein class") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(20, "mm")) 

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_C2_NMD_EJC_STAU_protein_copies_tissues_2024.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 1 - (D) UniProt-based Proteomics -----------------------------------------------------------
###

# Load UniProt-only proteomics data (Spectronaut)
UPF1_degron_proteomics_Reanalysis <- read_excel("Resources/Proteomics/2025_03_11_Reanalysis/UPF1_UniProt_WCP_protein_group_SN.xlsx") %>% 
  janitor::clean_names()

# Define names for individual heatmaps
NMD_names <- c("UPF1",
               "UPF2",
               "UPF3A",
               "UPF3B",
               "SMG1",
               "SMG5",
               "SMG6",
               "SMG7",
               "SMG8",
               "SMG9")

EJC_names <- c("EIF4A3",
               "MAGOH;MAGOHB",
               "MAGOH",
               "MAGOHB",
               "RBM8A",
               "CASC3")

STAU_names <- c("STAU1",
                "STAU2")

Other_names <- c("DHX34",
                 "MOV10",
                 "G3BP1",
                 "G3BP2",
                 "GSPT1",
                 "GSPT2",
                 "ETF1",
                 "NBAS",
                 "ATM",
                 "SLBP",
                 "PNRC2",
                 "ATR",
                 "ZC3H12A",
                 "TRIM71",
                 "GW182",
                 "AGO2",
                 "DCP1A",
                 "DCP1B",
                 "XRN1"
)

Names_combined <- c(NMD_names,
                    EJC_names,
                    STAU_names)

Names_combined_other <- c("UPF1",
                          Other_names)

### Other-UPF1 interactors -----------------------------------------------------------------

UPF1_degron_proteomics_Reanalysis %>% 
  filter(genes %in% Names_combined_other) %>% 
  filter(!id %in% c("STAU1_O95793-3",
                    "SMG1_H3BPS6")) %>% 
  dplyr::select(genes,
                starts_with("log_fc"),
                starts_with("adj_p")) %>% 
  dplyr::rename("UPF1_06h_MS_log2FC" = "log_fc_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_log2FC" = "log_fc_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_log2FC" = "log_fc_24h_nmd_inhibition_over_24h_dmso",
                "UPF1_06h_MS_padj" = "adj_p_val_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_padj" = "adj_p_val_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_padj" = "adj_p_val_24h_nmd_inhibition_over_24h_dmso") %>% 
  pivot_longer(cols = c(starts_with("UPF1_")),
               names_to = c("comparison",".value"),
               names_sep = "_MS_") %>% 
  mutate(genes = fct_rev(fct_relevel(genes,
                                     "UPF1",
                                     "GSPT1",
                                     "GSPT2",
                                     "ETF1",
                                     "DHX34",
                                     "NBAS",
                                     "MOV10",
                                     "G3BP1",
                                     "G3BP2",
                                     "DCP1A",
                                     "DCP1B",
                                     "XRN1",
                                     "ATM",
                                     "SLBP",
                                     "PNRC2",
                                     "ATR",
                                     "ZC3H12A",
                                     "TRIM71",
                                     "GW182",
                                     "AGO2"
  ))) %>% 
  mutate(comparison = (fct_relevel(comparison,
                                   "UPF1_06h",
                                   "UPF1_12h",
                                   "UPF1_24h"))) %>% 
  ggplot(aes(x=comparison,
             y=genes,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    size=-log10(padj)),
    stroke=0.2,
    color="black",
    shape=21) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2.5)) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(3*3+0.6, "mm"),
                   rows = unit(15*3+0.6, "mm"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_D1_Other_interactors_log2FC_flipped.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### NMD_EJC_SMD -----------------------------------------------------------------

UPF1_degron_proteomics_Reanalysis %>% 
  filter(genes %in% Names_combined) %>% 
  filter(!id %in% c("STAU1_O95793-3",
                    "SMG1_H3BPS6")) %>% 
  dplyr::select(genes,
                starts_with("log_fc"),
                starts_with("adj_p")) %>% 
  dplyr::rename("UPF1_06h_MS_log2FC" = "log_fc_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_log2FC" = "log_fc_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_log2FC" = "log_fc_24h_nmd_inhibition_over_24h_dmso",
                "UPF1_06h_MS_padj" = "adj_p_val_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_padj" = "adj_p_val_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_padj" = "adj_p_val_24h_nmd_inhibition_over_24h_dmso") %>% 
  pivot_longer(cols = c(starts_with("UPF1_")),
               names_to = c("comparison",".value"),
               names_sep = "_MS_") %>% 
  mutate(genes = fct_rev(fct_relevel(genes,
                                     "UPF1",
                                     "UPF2",
                                     "UPF3B",
                                     "SMG1",
                                     "SMG5",
                                     "SMG6",
                                     "SMG7",
                                     "SMG8",
                                     "SMG9",
                                     "EIF4A3",
                                     "MAGOH",
                                     "MAGOHB",
                                     "RBM8A",
                                     "CASC3",
                                     "STAU1",
                                     "STAU2"))) %>% 
  mutate(comparison = (fct_relevel(comparison,
                                   "UPF1_06h",
                                   "UPF1_12h",
                                   "UPF1_24h"))) %>% 
  ggplot(aes(x=comparison,
             y=genes,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    size=-log10(padj)),
    stroke=0.2,
    color="black",
    shape=21) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(0.5, 2.5)) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(3*3+0.6, "mm"),
                   rows = unit(16*3+0.6, "mm"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_D2_NMD_factors_log2FC_flipped.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### UPF1 peptides -----------------------------------------------------------

# Load peptide data
UPF1_degron_20250312_112231_peptides <- read_delim("Resources/Proteomics/2025_03_11_Reanalysis/UPF1_UniProt_WCP_peptides_SN.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE) %>% 
  clean_names()

# Prepare for plotting
distinct_peptides_SN_UPF1 <- UPF1_degron_20250312_112231_peptides %>% 
  filter(gene_names == "UPF1") %>% 
  dplyr::rename("UPF1_12h_DMSO_A_MS_intensity" = "intensity_a01_12h_dmso_a",
                "UPF1_12h_DMSO_B_MS_intensity" = "intensity_a02_12h_dmso_b",
                "UPF1_12h_DMSO_C_MS_intensity" = "intensity_a03_12h_dmso_c",
                "UPF1_12h_DMSO_D_MS_intensity" = "intensity_a04_12h_dmso_d",
                "UPF1_12h_IAA_A_MS_intensity" = "intensity_a05_12h_depl_a",
                "UPF1_12h_IAA_B_MS_intensity" = "intensity_a06_12h_depl_b",
                "UPF1_12h_IAA_C_MS_intensity" = "intensity_a07_12h_depl_c",
                "UPF1_12h_IAA_D_MS_intensity" = "intensity_a08_12h_depl_d",
                "UPF1_24h_DMSO_A_MS_intensity" = "intensity_a09_24h_dmso_a",
                "UPF1_24h_DMSO_B_MS_intensity" = "intensity_a10_24h_dmso_b",
                "UPF1_24h_DMSO_C_MS_intensity" = "intensity_a11_24h_dmso_c",
                "UPF1_24h_DMSO_D_MS_intensity" = "intensity_a12_24h_dmso_d",
                "UPF1_24h_IAA_A_MS_intensity" = "intensity_b01_24h_depl_a",
                "UPF1_24h_IAA_B_MS_intensity" = "intensity_b02_24h_depl_b",
                "UPF1_24h_IAA_C_MS_intensity" = "intensity_b03_24h_depl_c",
                "UPF1_24h_IAA_D_MS_intensity" = "intensity_b04_24h_depl_d",
                "UPF1_06h_DMSO_A_MS_intensity" = "intensity_b05_6h_dmso_a",
                "UPF1_06h_DMSO_B_MS_intensity" = "intensity_b06_6h_dmso_b",
                "UPF1_06h_DMSO_C_MS_intensity" = "intensity_b07_6h_dmso_c",
                "UPF1_06h_DMSO_D_MS_intensity" = "intensity_b08_6h_dmso_d",
                "UPF1_06h_IAA_A_MS_intensity" = "intensity_b09_6h_depl_a",
                "UPF1_06h_IAA_B_MS_intensity" = "intensity_b10_6h_depl_b",
                "UPF1_06h_IAA_C_MS_intensity" = "intensity_b11_6h_depl_c",
                "UPF1_06h_IAA_D_MS_intensity" = "intensity_b12_6h_depl_d") %>% 
  pivot_longer(cols = c(starts_with("UPF1_")),
               names_to = c("sample",".value"),
               names_sep = "_MS_") %>% 
  # remove outlier sample
  filter(sample != "UPF1_12h_IAA_A") %>% 
  mutate(condition = case_when(str_detect(sample, "12h_DMSO") ~ "12h_DMSO",
                               str_detect(sample, "12h_IAA") ~ "12h_IAA",
                               str_detect(sample, "24h_DMSO") ~ "24h_DMSO",
                               str_detect(sample, "24h_IAA") ~ "24h_IAA",
                               str_detect(sample, "06h_DMSO") ~ "06h_DMSO",
                               str_detect(sample, "06h_IAA") ~ "06h_IAA"))

#### Plot --------------------------------------------------------------------

distinct_peptides_SN_UPF1 %>% 
  dplyr::arrange(desc(intensity)) %>% 
  mutate(sequence = fct_rev(fct_inorder(sequence))) %>% 
  mutate(sample = fct_rev(sample)) %>% 
  ggplot(aes(y=sample,
             x=sequence,
             fill=log2(intensity))) +
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
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_tile() +
  scale_fill_viridis_c(option="mako",
                       na.value = "white") +
  labs(x="peptides",
       y="") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(30, "mm"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_D3_Proteomics_UPF1_peptides.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Check fraction of replicates with NA peptide intensities
distinct_peptides_SN_UPF1 %>% 
  group_by(condition) %>% 
  dplyr::count(is.na(intensity)) %>% 
  mutate(n_per = n/sum(n)) %>% 
  ungroup() %>% 
  filter(`is.na(intensity)` == TRUE)

### Violin plot of significant proteins -------------------------------------------------------

# Obtain tidy data for log2FC and p.adj per condition and ID
UPF1_degron_proteomics_Reanalysis_all <- UPF1_degron_proteomics_Reanalysis %>% 
  dplyr::select(id,
                starts_with("log_fc"),
                starts_with("adj_p")) %>% 
  dplyr::rename("UPF1_06h_MS_log2FC" = "log_fc_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_log2FC" = "log_fc_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_log2FC" = "log_fc_24h_nmd_inhibition_over_24h_dmso",
                "UPF1_06h_MS_padj" = "adj_p_val_06h_nmd_inhibition_over_06h_dmso",
                "UPF1_12h_MS_padj" = "adj_p_val_12h_nmd_inhibition_over_12h_dmso",
                "UPF1_24h_MS_padj" = "adj_p_val_24h_nmd_inhibition_over_24h_dmso") %>% 
  pivot_longer(cols = c(starts_with("UPF1_")),
               names_to = c("comparison",".value"),
               names_sep = "_MS_") %>% 
  separate(id, c("id_name", "id_ID"),sep = "_",remove = FALSE, extra = "merge") %>% 
  filter(comparison %in% c("UPF1_06h",
                           "UPF1_12h",
                           "UPF1_24h"))

# How many significant proteins
UPF1_degron_proteomics_Reanalysis_all_numbers <- UPF1_degron_proteomics_Reanalysis_all %>% 
  distinct(comparison,id, .keep_all = TRUE) %>% 
  filter(padj < 0.001 & abs(log2FC) > 0.5) %>% 
  mutate(up_down = case_when(log2FC > 0 & padj < 0.001 ~ "up",
                             log2FC < 0 & padj < 0.001 ~ "down",
                             TRUE ~ "ns")) %>% 
  count(comparison, up_down)

# Violin plot
UPF1_degron_proteomics_Reanalysis_all %>% 
  distinct(comparison,id_ID, .keep_all = TRUE) %>% 
  filter(padj < 0.001 & abs(log2FC) > 0.5) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "UPF1_06h",
                                          "UPF1_12h",
                                          "UPF1_24h"))) %>% 
  ggplot(aes(x=log2FC,
             y=comparison)) +
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
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "darkgray") +
  # geom_vline(xintercept = 0.5, linewidth = 0.2, color = "#963B5A") +
  # geom_vline(xintercept = -0.5, linewidth = 0.2, color = "#963B5A") +
  geom_violin(aes(fill=comparison),
              linewidth=0.2) +
  geom_point(data = (UPF1_degron_proteomics_Reanalysis_all %>% filter(id_name == "UPF1")), aes(color = case_when(id_name == "UPF1" ~ "#963B5A",
                                                                                                                TRUE ~ "darkgray")),
             alpha=0.75) +
  scale_color_identity() +
  scale_fill_manual(values = c("UPF1_06h" = "#11897C",
                               "UPF1_12h" = "#00696D",
                               "UPF1_15h" = "#006268",
                               "UPF1_18h" = "#005c64",
                               "UPF1_24h" = "#005560")) +
  guides(fill="none") +
  geom_label(data=UPF1_degron_proteomics_Reanalysis_all_numbers,
             aes(x=4,
                 label=case_when(up_down == "up" ~ paste0("up: n = ",n),
                                 up_down == "down" ~ "")),
             size = 4*0.36,
             color="black",
             nudge_y = 0.4,
             # label.padding = unit(0.1, "lines"),
             label.size = 0,
             parse=F
  ) +
  geom_text(data=UPF1_degron_proteomics_Reanalysis_all_numbers,
            aes(x=-4,
                label=case_when(up_down == "down" ~ paste0("down: n = ",n),
                                up_down == "up" ~ "")),
            size = 4*0.36,
            color="black",
            nudge_y = 0.4,
            # label.padding = unit(0.1, "lines"),
            # label.size = 0,
            parse=F
  ) +
  scale_x_continuous(limits=c(-5,5)) +
  labs(subtitle="padj < 0.001 & |log2FC| > 0.5",
       y="",
       x="log2FC\n(IAA/DMSO)") +
  force_panelsizes(cols = unit(45, "mm"),
                   rows = unit(15, "mm"))

ggsave(file.path("Plots", "Figure1", "Rev_1_F1_D4_Proteomics_violin_significant.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
