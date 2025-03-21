#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_RevisionAnalyses.R
# Objective: Code for generating Revision Analyses for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_NMDRHT_Analysis"

###
# Correlation -----------------------------------------------------------
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
## Transcript-level --------------------------------------------------------
###

# Read EZbakR data
UPF1_NMDRHT_EZbakR_TEC_combined <- read_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))

# Read data - if necessary
NMDRHT.v1.2_MainTable_forTable <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_MainTable.csv") %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>% 
              filter(condition == "UPF1_12h") %>% 
              dplyr::select(transcript_id, L2FC_kdeg_tx, padj_kdeg_tx, kdeg_tx_conclusion)) %>% 
  dplyr::select(-c(DTE_cluster_up, DTE_cluster_down))

### Description -------------------------------------------------------------

NMDRHT.v1.2_MainTable_description <- tibble(Field=colnames(NMDRHT.v1.2_MainTable_forTable),
                                            Description=c("NMDRHT.v1.2 gene id - should match GENCODE gene id",
                                                          "NMDRHT.v1.2 transcript id - unique identified of NMDRHT annotation",
                                                          "NMDRHT.v1.2 gene name - should match GENCODE gene name",
                                                          "NMDRHT.v1.2 transcript name - if full-splice match: GENCODE gene name & structural_category_simple; otherwise NMDRHT.v1.2 transcript id & structural_category_simple",
                                                          "NMDRHT.v1.2 gene biotype - should match GENCODE biotype",
                                                          "NMDRHT.v1.2 transcript biotype - protein_coding (50-nt rule = FALSE); nonsense_mediated_decay = (50-nt rule = TRUE); lncRNA = no ORF detected/predicted",
                                                          "Best matching GENCODE.v42 gene id with version number (see: https://www.gencodegenes.org/pages/data_format.html)",
                                                          "Best matching GENCODE.v42 transcript id with version number (see: https://www.gencodegenes.org/pages/data_format.html)",
                                                          "Best matching GENCODE.v42 gene name",
                                                          "Best matching GENCODE.v42 transcript name",
                                                          "Best matching GENCODE.v42 gene biotype (see: https://www.gencodegenes.org/pages/biotypes.html)",
                                                          "Best matching GENCODE.v42 transcript biotype (see: https://www.gencodegenes.org/pages/biotypes.html)",
                                                          "Best matching GENCODE.v42 transcript support level tag (see: https://www.gencodegenes.org/pages/data_format.html)",
                                                          "Best matching GENCODE.v42 ensembl canonical tag - is MANE_Select in most cases (see: https://www.gencodegenes.org/pages/tags.html)",
                                                          "Locus ID given by gffcompare when comparing different transcriptome assemblies - used as proxy for gene id",
                                                          "Transcript genomic coordinates - GRCh38/hg38",
                                                          "SQANTI3-determined structural category (see: https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories)",
                                                          "Simplified SQANTI3 structural category",
                                                          "SQANTI3-determined structural sub-category",
                                                          "gffcompare class code - similar to SQANTI3 structural category but different (see: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#transfrag-class-codes)",
                                                          "Simplified gffcompare class code",
                                                          "Absolute number of unique intron chain (UIC) support - i.e. in how many transcriptomes was this transcript detected?",
                                                          "UIC support based on only long-read transcriptome assemblies",
                                                          "UIC support based on only short-read transcriptome assemblies",
                                                          "UIC filter level - high: gene with well-supported transcripts; low: gene with poorly supported transcripts",
                                                          "(from SQANTI3 wiki): TRUE if the transcript start site is within a CAGE Peak",
                                                          "(from SQANTI3 wiki): TRUE if the transcript start site is within a polyA site",
                                                          "(from SQANTI3 wiki): distance to closest TSS based on CAGE Peak data",
                                                          "(from SQANTI3 wiki): distance to the closest polyA site",
                                                          "(from SQANTI3 wiki):TRUE if a polyA motif is detected in the 3'end of the transcript.",
                                                          "NMD-related ORF status",
                                                          "Data source for NMD-related ORF status determination",
                                                          "Most likely reason for NMD-sensitivity. Based on ORF status, Start/Stop codon matching reference and 50-nucleotide rule",
                                                          "50-nucleotide rule: if the stop codon is more than 50 nucleotides upstream of the last exon-exon junction = TRUE",
                                                          "Distance between stop codon and last exon junction. Counted from the first nucleotide after the stop codon to the first nucleotide of the intron",
                                                          "Number of downstream exon-exon junctions",
                                                          "Length of 3' UTR",
                                                          "ORFquant ORF id",
                                                          "ORF genomic coordinates - GRCh38/hg38",
                                                          "Absolute number of detected Ribo-Seq P-sites for this ORF",
                                                          "Data source for ORFquant ORF status determination",
                                                          "Percentage of selected ORF length-normalized P-sites, determined per transcript",
                                                          "Percentage of NMD-activating ORF length-normalized P-sites, determined per transcript. If 100 = all ORFs are NMD-activating; if 0 = no ORFs are NMD-activating",
                                                          "log2FC of length-normalized P-sites; 12h-vs-0h IAA; missing values imputed with length-normalized P-sites = 0.001",
                                                          "GENCODE canonical transcript used for ORF prediction via ORFanage",
                                                          "ORF classification based on matching predicted start and stop codons to canonical ORF",
                                                          "Start codon matching to canonical ORF",
                                                          "Stop codon matching to canonical ORF",
                                                          "ORFanage ORF id",
                                                          "Number of exons in this transcripts",
                                                          "Length of transcript",
                                                          "Length of CDS/ORF",
                                                          "Length of 5'UTR",
                                                          "Length of 3'UTR",
                                                          "GC content (%) of transcript",
                                                          "GC content (%) of 5'UTR",
                                                          "GC content (%) of CDS/ORF",
                                                          "GC content (%) of 3'UTR",
                                                          "predicted minimum thermodynamic free energy (–ΔG/nt) of the 3’UTR",
                                                          "predicted length-normalized, minimum thermodynamic free energy (–ΔG/nt) of the 3’UTR",
                                                          "Differential transcript expression (DTE) cluster; determined by hierachical clustering",
                                                          "For NMD relevance determination; is the transcript up- or downregulated?",
                                                          "Absolute Number of RNA-Seq conditions (out of 18) in which this transcript is concordantly and significantly up- or downregulated",
                                                          "Percentage (scaled between 0-100) of RNA-Seq conditions in which this transcript is concordantly and significantly up- or downregulated",
                                                          "Binned NMD relevance (n=4)",
                                                          "Median DTE log2FC of concordant and significant NMD relevance conditions",
                                                          "Median DTE adjusted p-value of concordant and significant NMD relevance conditions",
                                                          "Adjusted p-value of ImpulseDE2 modeling",
                                                          "Significant ImpulseDE2 modeling? If TRUE = sig. otherwise = n.s.",
                                                          "Best ImpulseDE2 fit, either Impulse, Sigmoid or noFit",
                                                          "ImpulseDE2 slope parameter of both first and second transition",
                                                          "ImpulseDE2 initial amplitude",
                                                          "ImpulseDE2 peak amplitude",
                                                          "ImpulseDE2 steady state amplitude",
                                                          "ImpulseDE2 onset time",
                                                          "ImpulseDE2 offset time",
                                                          "EZbakR log2FC of degradation rate constant for transcript-level analysis",
                                                          "EZbakR adjusted p-value of degradation rate constant for transcript-level analysis",
                                                          "EZbakR conclusion about degradation rate for transcript-level analysis"
                                            ),
                                            Source=c(rep("NMDRHT.v1.2", 6),
                                                     rep("GENCODE.v42", 8),
                                                     rep("NMDRHT.v1.2", 2),
                                                     rep("SQANTI3", 3),
                                                     rep("gffcompare", 2),
                                                     rep("NMDRHT.v1.2", 4),
                                                     rep("SQANTI3", 5),
                                                     rep("NMD sensitivity", 7),
                                                     rep("ORFquant", 7),
                                                     rep("ORFanage", 5),
                                                     rep("NMDRHT.v1.2", 11),
                                                     "Hierachical clustering",
                                                     rep("NMD relevance", 6),
                                                     rep("ImpulseDE2", 9),
                                                     rep("SLAM-Seq", 3)
                                            ))

# Create a gene-level workbook and add description sheet
UPF1_NMDRHT_NMDRHT.v1.2_MainTable <- createWorkbook()

addWorksheet(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, "Description")

#### Title -------------------------------------------------------------------

writeData(UPF1_NMDRHT_NMDRHT.v1.2_MainTable,
          "Description", 
          "Boehm et al. 2025 - Supplemental Table 3 - Transcript-level analyses", 
          startCol = 1, 
          startRow = 1)

mergeCells(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
           "Description", 
           cols = 1:3, rows = 1)

# Define title style
title_style <- createStyle(
  fontName = "Calibri",
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#32475B",
  textDecoration = "bold",
  border = "TopBottom",
  borderColour = "#FFFFFF"
)

# add style to Title
addStyle(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
         "Description", 
         title_style, rows = 1, cols = 1:3)

setColWidths(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
             "Description", 
             cols = 1:1, widths = 35)

setColWidths(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
             "Description", 
             cols = 2:2, widths = 140)

setColWidths(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
             "Description", 
             cols = 3:3, widths = 30)

#### Subtitle -------------------------------------------------------------------
writeData(UPF1_NMDRHT_NMDRHT.v1.2_MainTable,
          "Description", 
          "Description of all field of transcript-level analyses", 
          startCol = 1, 
          startRow = 2)

mergeCells(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
           "Description", 
           cols = 1:3, rows = 2)

# Define title style
subtitle_style <- createStyle(
  fontName = "Calibri",
  fontSize = 11,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#287DAB",
  border = "TopBottom",
  borderColour = "#FFFFFF"
)

# add style to Title
addStyle(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
         "Description", 
         subtitle_style, rows = 2, cols = 1:3)

# Write data to the sheet
writeData(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
          "Description",
          startCol = 1, startRow = 3,
          NMDRHT.v1.2_MainTable_description, 
          headerStyle = title_style)

### Data --------------------------------------------------------------------

addWorksheet(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, "NMDRHT_MainTable")

# Define alt title style
alt_title_style <- createStyle(
  fontName = "Calibri",
  fontSize = 11,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#32475B",
  textDecoration = "bold",
  border = "TopBottom",
  borderColour = "#FFFFFF"
)

# Write data to the sheet
writeData(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
          "NMDRHT_MainTable",
          startCol = 1, startRow = 1,
          NMDRHT.v1.2_MainTable_forTable, 
          headerStyle = alt_title_style)

# Set automatic column width
setColWidths(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, 
             "NMDRHT_MainTable", 
             cols = 1:ncol(NMDRHT.v1.2_MainTable_forTable), 
             widths = "auto")

# Save the workbook
saveWorkbook(UPF1_NMDRHT_NMDRHT.v1.2_MainTable, "Tables/Boehm_et_al_2025_Revision1_Table_S3.xlsx", overwrite = TRUE)