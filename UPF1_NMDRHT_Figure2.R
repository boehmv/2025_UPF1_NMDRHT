#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Figure2
# Objective: Code for generating Panels of Figure 2 + Figure S2 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_NMDRHT_Analysis"

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")

###
## Rev_1 - Figure 2 - (A1) DGE -----------------------------------------------------------
###

# Comment: Plot F2_A1 = Overview DGE figure

### Prepare for Plot -----------------------------------------------------------

Rev_1_F2_A1_DESeq2_DGE_combined_forPlot <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c("HEK293_UPF1_KD_Kishor_2019",
                              "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
                              "HeLa_UPF1_KD_Longman_2020",
                              "HEK293_UPF1_KD_Fritz_2022",
                              "K562_UPF1_KD_Hug_2022",
                              "HUH7_UPF1_KD_Lee_2022",
                              "HepG2_UPF1_KD_He_2023",
                              "HEK293_SMG567_KO_KD_Boehm_2021",
                              "HEK293_SMG67_KD_Britto_Borges_2024",
                              "HeLa_SMG67_KD_Britto_Borges_2024",
                              "MCF7_SMG67_KD_Britto_Borges_2024",
                              "U2OS_SMG67_KD_Britto_Borges_2024",
                              # "HEK293_UPF3_dKO_Wallmeroth_2022",
                              # "HEK293_CASC3_KO_Gerbracht_2020",
                              "HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study",
                              "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                              "HFF_SMG1i_this_Study",
                              "HUVEC_SMG1i_this_Study")) %>% 
  filter(!condition_2 %in% c("SMG8_KO_0uM",
                             "SMG9_KO_0uM",
                             "SMG8_delKID_0uM",
                             "SMG8_KO_01uM",
                             "SMG9_KO_01uM",
                             "SMG8_delKID_01uM",
                             "SMG8_KO_1uM",
                             "SMG9_KO_1uM",
                             "SMG8_delKID_1uM")) %>% 
  mutate(condition_2 = case_when(condition_2 == "UPF1" ~ publicationName,
                                 condition_2 == "SMG1i" ~ experimentSet,
                                 TRUE ~ condition_2)) %>% 
  mutate(experiment_set = fct_relevel(experimentSet,
                                      "HEK293_SMG567_KO_KD_Boehm_2021",
                                      "HEK293_SMG67_KD_Britto_Borges_2024",
                                      "HeLa_SMG67_KD_Britto_Borges_2024",
                                      "MCF7_SMG67_KD_Britto_Borges_2024",
                                      "U2OS_SMG67_KD_Britto_Borges_2024",
                                      "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                                      "HFF_SMG1i_this_Study",
                                      "HUVEC_SMG1i_this_Study")) %>% 
  arrange(experiment_set) %>% 
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

#### Rev_1_F2_A1_Heatmap -----------------------------------------------------------

Rev_1_F2_A1_DESeq2_DGE_combined_forPlot %>% 
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
  force_panelsizes(rows = unit(length(unique(Rev_1_F2_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.4, "mm"),
                   cols = unit(12+0.4, "mm"))  +
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

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_A1_DESeq2_DGE_Heatmap_perType.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### Rev_1_F2_A1_Barplot -----------------------------------------------------------

Rev_1_F2_A1_DESeq2_DGE_combined_forPlot %>% 
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
        panel.grid.minor = element_blank(),
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
  geom_col(color="black",
           linewidth = 0.2) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "other" = "#D9DADA")) +
  geom_vline(xintercept = 0,
             color="#963B5A",
             linewidth=0.25) +
  # scale_fill_viridis_d(option="G",
  #                      begin=0.25,
  #                      end=0.75) +
  labs(y="",
       x="Number of \nsig. regulated genes",
       fill="Gene biotype",
       size="Sig. DGE events") +
  force_panelsizes(rows = unit(length(unique(Rev_1_F2_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.4, "mm"),
                   cols = unit(20, "mm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE,
                           title.position = "top"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_A1_DESeq2_DGE_Barplot_perType.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### Rev_1_F2_A1_Examples -----------------------------------------------------------

# Define representative NMD-targeted genes plus UPF1
standard_heatmap_genes <- c("ZFAS1",
                            "SNHG12",
                            "GAS5",
                            "GADD45B",
                            "UPF1")

Rev_1_F2_A1_DESeq2_DGE_combined_forPlot %>%  
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
  force_panelsizes(rows = unit(length(unique(Rev_1_F2_A1_DESeq2_DGE_combined_forPlot$condition_2))*2+0.4, "mm"),
                   cols = unit(5*2+0.4, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_A1_DESeq2_DGE_Examples.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 2 - (A2) QC -----------------------------------------------------------
###

# Load rRNA mapping QC multiQC output
rRNA_QC  <-  read_csv(file.path("Resources/QC", "rRNA_QC.csv"))

# Calculate mean overall alignment rate of rRNA
Rev_1_F2_A2_UPF1_rRNA_QC_summary <- rRNA_QC %>% 
  group_by(publicationName) %>% 
  summarize(meanAlignment = mean(general.overall_alignment_rate),
            sdAlignment = sd(general.overall_alignment_rate)) %>% 
  filter(!is.na(publicationName))

# Plot
Rev_1_F2_A2_UPF1_rRNA_QC <- Rev_1_F2_A2_UPF1_rRNA_QC_summary %>% 
  mutate(sample="rRNA") %>% 
  ggplot(aes(x=sample,
             y=fct_rev(publicationName),
             fill=meanAlignment)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(option = "magma")(9),
                       limits=c(0, 100),
                       na.value = "#2B2A29") +
  geom_text(aes(label=round(meanAlignment,0),
                color=ifelse(meanAlignment<50, 'white', 'black')),
            size = 5*0.36) +
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
       y="") +
  coord_fixed(ratio=1) +
  scale_color_identity() +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = "none") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(Rev_1_F2_A2_UPF1_rRNA_QC_summary$publicationName)*2+0.4, "mm"),
                   cols = unit(2+0.4, "mm"))


### Load Salmon QC data ---------------------------------------------------------------

# Load rRNA mapping QC multiQC output
UPF1_Salmon_QC  <-  read_csv(file.path("Resources/QC", "UPF1_Salmon_QC.csv"))

# Calculate mean total counts and percent mapped per experimentSet
Rev_1_F2_A2_UPF1_Salmon_QC_summary <- UPF1_Salmon_QC %>% 
  group_by(publicationName) %>% 
  summarize(reads = mean(general.num_processed),
            percent_mapped = mean(general.percent_mapped),
            libType = paste(unique(general.library_types), collapse = ', ')) %>% 
  mutate(stranded = case_when(libType %in% c("ISR", "ISF") ~ "S",
                              libType %in% c("IU") ~ "U",
                              TRUE ~ "other")) %>% 
  filter(!is.na(publicationName))

# Plot
Rev_1_F2_A2_UPF1_Salmon_QC <- UPF1_Salmon_QC %>% 
  dplyr::select(metadata.sample_id,experimentSet,condition_2, publicationName, general.num_processed, general.num_mapped, general.percent_mapped) %>% 
  mutate(general.num_unmapped = general.num_processed - general.num_mapped) %>% 
  dplyr::rename("total" = "general.num_processed",
                "mapped" = "general.num_mapped",
                "unmapped" = "general.num_unmapped",
                "percent_mapped" = "general.percent_mapped") %>% 
  pivot_longer(cols=c(total, mapped, unmapped),
               names_to = "class",
               values_to = "reads") %>% 
  filter(!is.na(publicationName)) %>% 
  filter(!class %in% c("mapped", "unmapped")) %>% 
  ggplot(aes(x=reads,
             y=fct_rev(publicationName))) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "top",
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  ggdist::stat_pointinterval(aes(interval_color = after_stat(level)),
                     .width = c(0.9),
                     point_interval = "mean_qi",
                     linewidth = 0.5,
                     point_color = "black",
                     point_size = 0.5,
  ) +
  geom_point(aes(fill=percent_mapped),
             color="black",
             position=position_jitter(height = 0.25),
             stroke=0.25,
             size=0.5,
             alpha=0.75,
             shape = 21) +
  geom_label(data=Rev_1_F2_A2_UPF1_Salmon_QC_summary,
             aes(x=10^8.5,
                 label=paste(round(reads/10^6,0))),
             size = 4*0.36,
             fill="gray",
             label.padding = unit(0.1, "lines"),
             label.size = 0,
             parse=F
  ) +
  geom_label(data=Rev_1_F2_A2_UPF1_Salmon_QC_summary,
             aes(x=10^8.7,
                 fill=percent_mapped,
                 color=ifelse(percent_mapped<40, 'white', 'black'),
                 label=paste(round(percent_mapped,0),"%", sep = "")),
             size = 4*0.36,
             label.padding = unit(0.1, "lines"),
             label.size = 0,
             parse=F
  ) +
  geom_label(data=Rev_1_F2_A2_UPF1_Salmon_QC_summary %>% filter(stranded == "U"),
             aes(x=10^8.9,
                 label=paste(stranded)),
             size = 4*0.36,
             fill="lightgray",
             color="black",
             label.padding = unit(0.1, "lines"),
             label.size = 0,
             parse=F
  ) +
  geom_label(data=Rev_1_F2_A2_UPF1_Salmon_QC_summary %>% filter(stranded == "S"),
             aes(x=10^8.9,
                 label=paste(stranded)),
             size = 4*0.36,
             fill="black",
             color="white",
             label.padding = unit(0.1, "lines"),
             label.size = 0,
             parse=F
  ) +
  scale_color_identity() +
  scale_color_manual(
    values = c("#636363"),
    aesthetics = "interval_color",
    guide = guide_legend(override.aes = list(point_fill = NA))
  ) +
  scale_fill_gradientn(colors = viridis_pal()(9), limits=c(0, 100)) +
  scale_x_log10(limits = c(10^7,10^9),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x="Number of reads",
       y="",
       interval_color="Mean +\nProbability interval",
       fill="Percent mapped (%)\n(salmon)") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.vjust = 0.5)) +
  force_panelsizes(rows = unit(length(Rev_1_F2_A2_UPF1_Salmon_QC_summary$stranded)*2+0.4, "mm"))

### Check for non-polyA genes -----------------------------------------------
# Publication: PMID: 21324177
# Select 5 genes ("H3C2", "RPPH1", "RMRP", "TERC", "SNORD3A")

Rev_1_F2_A2_UPF1_noncodingRNA <- Swish_DGE_combined %>% 
  filter(!experimentSet %in% c("HEK293_UPF1_KD_Baird_2018", "HEK293_NMD_Boehm_2018")) %>% 
  mutate(experimentSet = fct_rev(factor(experimentSet, levels=levels(UPF1_NMDRHT_datasets$experimentSet)))) %>% 
  filter(gene_name %in% c("H3C2", "RPPH1", "RMRP", "TERC", "SNORD3A")) %>% 
  group_by(publicationName, gene_name) %>% 
  summarize(mean_log10mean = mean(log10mean)) %>% 
  ggplot(aes(x=gene_name,
             y=fct_rev(publicationName),
             fill=mean_log10mean)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_gradientn(colors = viridis_pal(option = "mako")(9),
                       limits=c(0, 6),
                       na.value = "#2B2A29") +
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
       y="") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title="mean log10\n expression",
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(Rev_1_F2_A2_UPF1_Salmon_QC_summary$stranded)*2+0.4, "mm"),
                   cols = unit(5*2+0.4, "mm"))


### Rev_1_F2_A2 final plot -----------------------------------------------------------

Rev_1_F2_A2_UPF1_Salmon_QC | Rev_1_F2_A2_UPF1_noncodingRNA | Rev_1_F2_A2_UPF1_rRNA_QC

ggsave(file.path("Plots", "Figure2", paste0("Rev_1_F2_A2_UPF1_Salmon_ncRNA_QC.pdf")),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 2 - (A3) PCA -----------------------------------------------------------
###

### PCA Degradation -----------------------------------------------------------

load("Resources/QC/PCA_UPF1_degradation.rds")

plt_plot_deg <- ggplot(plt_deg %>% mutate(condition =(fct_inorder(as_factor(condition)))), aes(PC1, PC2, fill=condition)) +
  theme_classic() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=6, colour = 'black'), 
        axis.title=element_text(size=6), 
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2,
             size=1) +
  ggrepel::geom_text_repel(data=plt_deg %>% distinct(condition, .keep_all = TRUE),
                           aes(label=condition,
                               color=condition),
                           box.padding = 0.0,
                           min.segment.length = Inf,
                           size=4*0.36) +
  scale_fill_manual(values=c("control\n0h" = "#A0A0A0",
                             "control\n48h" = "#2D2D2D",
                             "N-AID-\nUPF1\n0h" = "#65BB8C",
                             "N-AID-\nUPF1\n2h" = "#43A787",
                             "N-AID-\nUPF1\n4h" = "#229380",
                             "N-AID-\nUPF1\n8h" = "#007E78",
                             "N-AID-\nUPF1\n12h" = "#00696D",
                             "N-AID-\nUPF1\n24h" = "#005560",
                             "N-AID-\nUPF1\n48h" = "#0B4151")) +
  scale_color_manual(values=c("control\n0h" = "#A0A0A0",
                              "control\n48h" = "#2D2D2D",
                              "N-AID-\nUPF1\n0h" = "#65BB8C",
                              "N-AID-\nUPF1\n2h" = "#43A787",
                              "N-AID-\nUPF1\n4h" = "#229380",
                              "N-AID-\nUPF1\n8h" = "#007E78",
                              "N-AID-\nUPF1\n12h" = "#00696D",
                              "N-AID-\nUPF1\n24h" = "#005560",
                              "N-AID-\nUPF1\n48h" = "#0B4151")) +
  labs(x = paste0("PC1: ",percentVar_deg[1],"% variance"), 
       y = paste0("PC2: ",percentVar_deg[2],"% variance")) + 
  coord_fixed() +
  guides(fill = "none",
         color = "none") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols =  unit(30, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_A3_DESeq2_N-AID_Degradation_DGE_PCA.pdf"),
       plt_plot_deg,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### PCA Recovery -----------------------------------------------------------

load("Resources/QC/PCA_UPF1_recovery.rds")

plt_plot_rec <- ggplot(plt_rec %>% mutate(condition =(fct_inorder(as_factor(condition)))), aes(PC1, PC2, fill=condition)) +
  theme_classic() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=6, colour = 'black'), 
        axis.title=element_text(size=6), 
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2,
             size=1) +
  ggrepel::geom_text_repel(data=plt_rec %>% distinct(condition, .keep_all = TRUE),
                           aes(label=condition,
                               color=condition),
                           box.padding = 0.0,
                           min.segment.length = Inf,
                           size=4*0.36) +
  scale_fill_manual(values=c("control\n0h" = "#A0A0A0",
                             "N-AID-\nUPF1\n24h-R0h" = "#682714",
                             "N-AID-\nUPF1\n24h-R2h" = "#953C00",
                             "N-AID-\nUPF1\n24h-R4h" = "#BD5500",
                             "N-AID-\nUPF1\n24h-R8h" = "#E17100",
                             "N-AID-\nUPF1\n24h-R12h" = "#EF9300",
                             "N-AID-\nUPF1\n24h-R24h" = "#F5B300",
                             "N-AID-\nUPF1\n24h-R48h" = "#F8CF69")) +
  scale_color_manual(values=c("control\n0h" = "#A0A0A0",
                              "N-AID-\nUPF1\n24h-R0h" = "#682714",
                              "N-AID-\nUPF1\n24h-R2h" = "#953C00",
                              "N-AID-\nUPF1\n24h-R4h" = "#BD5500",
                              "N-AID-\nUPF1\n24h-R8h" = "#E17100",
                              "N-AID-\nUPF1\n24h-R12h" = "#EF9300",
                              "N-AID-\nUPF1\n24h-R24h" = "#F5B300",
                              "N-AID-\nUPF1\n24h-R48h" = "#F8CF69")) +
  labs(x = paste0("PC1: ",percentVar_rec[1],"% variance"), 
       y = paste0("PC2: ",percentVar_rec[2],"% variance")) + 
  coord_fixed() +
  guides(fill = "none",
         color = "none") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols =  unit(30, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_A3_DESeq2_N-AID_Recovery_DGE_PCA.pdf"),
       plt_plot_rec,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### PCA FKBP-HCT -----------------------------------------------------------

load("Resources/QC/PCA_UPF1_FKBP_HCT.rds")

plt_plot_FKBP_HCT <- ggplot(plt_FKBP_HCT %>% mutate(condition =(fct_inorder(as_factor(condition)))), aes(PC1, PC2, fill=condition)) +
  theme_classic() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=6, colour = 'black'), 
        axis.title=element_text(size=6), 
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2,
             size=1) +
  ggrepel::geom_text_repel(data=plt_FKBP_HCT %>% distinct(condition, .keep_all = TRUE),
                           aes(label=condition,
                               color=condition),
                           box.padding = 0.0,
                           min.segment.length = Inf,
                           size=4*0.36) +
  scale_fill_manual(values=c("control\n0h" = "#A0A0A0",
                             "HCT116\nFKBP-UPF1\n0h" = "#C8C3E2",
                             "HCT116\nFKBP-UPF1\n12h" = "#6B56A7")) +
  scale_color_manual(values=c("control\n0h" = "#A0A0A0",
                              "HCT116\nFKBP-UPF1\n0h" = "#C8C3E2",
                              "HCT116\nFKBP-UPF1\n12h" = "#6B56A7")) +
  labs(x = paste0("PC1: ",percentVar_FKBP_HCT[1],"% variance"), 
       y = paste0("PC2: ",percentVar_FKBP_HCT[2],"% variance")) + 
  coord_fixed() +
  guides(fill = "none",
         color = "none") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols =  unit(30, "mm"))

ggsave(file.path("Plots", "Figure2",  "Rev_1_F2_A3_DESeq2_FKBP_HCT_DGE_PCA.pdf"),
       plt_plot_FKBP_HCT,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### PCA FKBP-HEK -----------------------------------------------------------

load("Resources/QC/PCA_UPF1_FKBP_HEK.rds")

plt_plot_FKBP_HEK <- ggplot(plt_FKBP_HEK %>% mutate(condition =(fct_inorder(as_factor(condition)))), aes(PC1, PC2, fill=condition)) +
  theme_classic() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=6, colour = 'black'), 
        axis.title=element_text(size=6), 
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2,
             size=1) +
  ggrepel::geom_text_repel(data=plt_FKBP_HEK %>% distinct(condition, .keep_all = TRUE),
                           aes(label=condition,
                               color=condition),
                           box.padding = 0.0,
                           min.segment.length = Inf,
                           size=4*0.36) +
  scale_fill_manual(values=c("control\n0h" = "#A0A0A0",
                             "HEK293\nFKBP-UPF1\n0h" = "#DD8492",
                             "HEK293\nFKBP-UPF1\n12h" = "#963B5A")) +
  scale_color_manual(values=c("control\n0h" = "#A0A0A0",
                              "HEK293\nFKBP-UPF1\n0h" = "#DD8492",
                              "HEK293\nFKBP-UPF1\n12h" = "#963B5A")) +
  labs(x = paste0("PC1: ",percentVar_FKBP_HEK[1],"% variance"), 
       y = paste0("PC2: ",percentVar_FKBP_HEK[2],"% variance")) + 
  coord_fixed() +
  guides(fill = "none",
         color = "none") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols =  unit(30, "mm"))

ggsave(file.path("Plots", "Figure2",  "Rev_1_F2_A3_DESeq2_FKBP_HEK_DGE_PCA.pdf"),
       plt_plot_FKBP_HEK,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 2 - (B) Cluster Gene -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/Cluster/Rev_1_F2_B_gene_cluster_information.rds")

### Plot Heatmaps -----------------------------------------------------------

#### Up ----------------------------------------------------------------------

Rev_1_F2_B_gene_hclust_split_up <- factor(Rev_1_F2_B_gene_hclust_up_gene_cluster_fix, levels=c("up 1:early",
                                                                                               "up 2:delayed",
                                                                                               "up 3:late",
                                                                                               "up 4:inverse"))

Rev_1_F2_B_gene_split_col = rep(1:17, each = 3)

ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 0,
                             col = 0),
                   labels = unique(Rev_1_F2_B_coldata$condition),
                   labels_rot = 90,
                   labels_just = "right",
                   labels_offset = unit(1, "npc"),
                   labels_gp = gpar(fontsize = 6),
                   height = unit(1, "cm"),
                   which = c("column"))
)

# Color setting for Heatmap
Rev_1_F2_B_gene_col_up = colorRamp2(c(-6,-4,-2, -1, 0, 1, 2, 4, 6), rev(c("#88002D","#D4332A", "#F9834C", "#FFC08E", "white", "#B7C2DE", "#7797C9", "#0C6DA8", "#0E3F5C")))

# log2-transformed *Z-SCALED* heatmap of upregulated DGE genes
Rev_1_F2_B_gene_hclust_up_heat <- Heatmap(matrix = Rev_1_F2_B_gene_test_up_fix,
                                          col =  Rev_1_F2_B_gene_col_up,
                                          cluster_rows = FALSE,
                                          cluster_columns = FALSE,
                                          name="Z-score\n(log2-normalized\ncounts)",
                                          row_split=Rev_1_F2_B_gene_hclust_split_up, 
                                          column_split = Rev_1_F2_B_gene_split_col,
                                          cluster_row_slices = FALSE,
                                          row_dend_reorder = FALSE,
                                          bottom_annotation = ha,
                                          border = TRUE,
                                          border_gp = gpar(col = "black", lwd = 0.2),
                                          row_title_rot = 0,
                                          row_title = "%s", 
                                          use_raster = TRUE,
                                          raster_device = "CairoPNG",
                                          raster_by_magick = FALSE,
                                          show_row_names = FALSE, 
                                          show_column_names = FALSE,
                                          row_title_gp = gpar(fontsize=6),
                                          column_title_gp = gpar(fontsize = 6, fontface = "bold"),
                                          heatmap_legend_param = list(at = c(-6, -4, -2, 0, 2, 4, 6),
                                                                      title_gp = gpar(fontsize = 6), 
                                                                      labels_gp = gpar(fontsize = 6),
                                                                      legend_height = unit(2, "cm"),
                                                                      grid_width = unit(0.2, "cm")),
                                          width = unit((length(unique(Rev_1_F2_B_coldata$condition))*2)+(-1+length(unique(Rev_1_F2_B_coldata$condition))*1), "mm"),
                                          height = unit(6, "cm"),
                                          column_title = "GENCODE - DESeq2 DGE sig. upregulated")

Rev_1_F2_B_gene_hclust_up_heat = draw(Rev_1_F2_B_gene_hclust_up_heat)

pdf(file = file.path("Plots", "Figure2", "Rev_1_F2_B_only_gencode.v42_DESeq2_DGE_sig_upregulated_hclust_heatmap.pdf"),
    useDingbats = FALSE,
    family="Arial",
    fonts="Arial")

draw(Rev_1_F2_B_gene_hclust_up_heat)

dev.off()

#### Down ----------------------------------------------------------------------

Rev_1_F2_B_gene_hclust_split_down <- factor(Rev_1_F2_B_gene_hclust_down_gene_cluster_fix, levels=c("down 1:early",
                                                                                                   "down 2:delayed",
                                                                                                   "down 3:late",
                                                                                                   "down 4:inverse"))

Rev_1_F2_B_gene_split_col = rep(1:17, each = 3)

ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 0,
                             col = 0),
                   labels = unique(Rev_1_F2_B_coldata$condition),
                   labels_rot = 90,
                   labels_just = "right",
                   labels_offset = unit(1, "npc"),
                   labels_gp = gpar(fontsize = 6),
                   height = unit(1, "cm"),
                   which = c("column"))
)

# Color setting for Heatmap
Rev_1_F2_B_gene_col_down = colorRamp2(c(-6,-4,-2, -1, 0, 1, 2, 4, 6), rev(c("#88002D","#D4332A", "#F9834C", "#FFC08E", "white", "#B7C2DE", "#7797C9", "#0C6DA8", "#0E3F5C")))

# log2-transformed *Z-SCALED* heatmap of downregulated DGE genes
Rev_1_F2_B_gene_hclust_down_heat <- Heatmap(matrix = Rev_1_F2_B_gene_test_down_fix,
                                            col =  Rev_1_F2_B_gene_col_down,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            name="Z-score\n(log2-normalized\ncounts)",
                                            row_split=Rev_1_F2_B_gene_hclust_split_down, 
                                            column_split = Rev_1_F2_B_gene_split_col,
                                            cluster_row_slices = FALSE,
                                            row_dend_reorder = FALSE,
                                            bottom_annotation = ha,
                                            border = TRUE,
                                            border_gp = gpar(col = "black", lwd = 0.2),
                                            row_title_rot = 0,
                                            row_title = "%s", 
                                            use_raster = TRUE,
                                            raster_device = "CairoPNG",
                                            raster_by_magick = FALSE,
                                            # raster_quality = 2,
                                            # raster_resize_mat = TRUE,
                                            show_row_names = FALSE, 
                                            show_column_names = FALSE,
                                            row_title_gp = gpar(fontsize=6),
                                            column_title_gp = gpar(fontsize = 6, fontface = "bold"),
                                            heatmap_legend_param = list(at = c(-6, -4, -2, 0, 2, 4, 6),
                                                                        title_gp = gpar(fontsize = 6), 
                                                                        labels_gp = gpar(fontsize = 6),
                                                                        legend_height = unit(2, "cm"),
                                                                        grid_width = unit(0.2, "cm")),
                                            width = unit((length(unique(Rev_1_F2_B_coldata$condition))*2)+(-1+length(unique(Rev_1_F2_B_coldata$condition))*1), "mm"),
                                            height = unit(6, "cm"),
                                            column_title = "GENCODE - DESeq2 DGE sig. downregulated")

Rev_1_F2_B_gene_hclust_down_heat = draw(Rev_1_F2_B_gene_hclust_down_heat)

pdf(file = file.path("Plots", "Figure2", "Rev_1_F2_B_only_gencode.v42_DESeq2_DGE_sig_downregulated_hclust_heatmap.pdf"),
    useDingbats = FALSE,
    family="Arial",
    fonts="Arial")

draw(Rev_1_F2_B_gene_hclust_down_heat)

dev.off()

### Fit the log2-norm counts ------------------------------------------------

#### Up ------------------------------------------------

# Prepare for plotting
Rev_1_F2_B_gene_test_up_klus_cond_forPlot <- Rev_1_F2_B_gene_test_up_klus_cond %>% 
  pivot_longer(cols=-c(gene_id, DGE_cluster_up),
               names_to = "condition",
               values_to = "score") %>% 
  mutate(timepoint = case_when(condition == "control_0h" ~ -20,
                               condition == "control_48h" ~ -12,
                               condition == "UPF1_Nter_0h" ~ 0,
                               condition == "UPF1_Nter_2h" ~ 2,
                               condition == "UPF1_Nter_4h" ~ 4,
                               condition == "UPF1_Nter_8h" ~ 8,
                               condition == "UPF1_Nter_12h" ~ 12,
                               condition == "UPF1_Nter_24h" ~ 24,
                               condition == "UPF1_Nter_48h" ~ 48,
                               condition == "control_R0h" ~ 84,
                               condition == "UPF1_Nter_24h_R0h" ~ 24,
                               condition == "UPF1_Nter_24h_R2h" ~ 26,
                               condition == "UPF1_Nter_24h_R4h" ~ 28,
                               condition == "UPF1_Nter_24h_R8h" ~ 32,
                               condition == "UPF1_Nter_24h_R12h" ~ 36,
                               condition == "UPF1_Nter_24h_R24h" ~ 48,
                               condition == "UPF1_Nter_24h_R48h" ~ 72)) %>% 
  mutate(exp_cond = case_when(condition == "control_0h" ~ "control",
                              condition == "control_48h" ~ "control",
                              condition == "UPF1_Nter_0h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_2h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_4h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_8h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_12h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_24h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_48h" ~ "UPF1_AID",
                              condition == "control_R0h" ~ "control",
                              condition == "UPF1_Nter_24h_R0h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R2h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R4h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R8h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R12h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R24h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R48h" ~ "UPF1_AID_Rec")) %>% 
  mutate(condition = fct_relevel(condition,
                                 "control_0h",
                                 "control_48h",
                                 "UPF1_Nter_0h",
                                 "UPF1_Nter_2h",
                                 "UPF1_Nter_4h",
                                 "UPF1_Nter_8h",
                                 "UPF1_Nter_12h",
                                 "UPF1_Nter_24h",
                                 "UPF1_Nter_48h",
                                 "control_R0h",
                                 "UPF1_Nter_24h_R0h",
                                 "UPF1_Nter_24h_R2h",
                                 "UPF1_Nter_24h_R4h",
                                 "UPF1_Nter_24h_R8h",
                                 "UPF1_Nter_24h_R12h",
                                 "UPF1_Nter_24h_R24h",
                                 "UPF1_Nter_24h_R48h"))

Rev_1_F2_B_gene_test_up_klus_cond_n <- Rev_1_F2_B_gene_test_up_klus_cond_forPlot %>% 
  group_by(DGE_cluster_up) %>% 
  distinct(gene_id,.keep_all = TRUE) %>% 
  summarize(n=n()) %>% 
  dplyr::rename("DGE_cluster" = "DGE_cluster_up")

#### Down ------------------------------------------------

# Prepare for plotting
Rev_1_F2_B_gene_test_down_klus_cond_forPlot <- Rev_1_F2_B_gene_test_down_klus_cond %>% 
  pivot_longer(cols=-c(gene_id, DGE_cluster_down),
               names_to = "condition",
               values_to = "score") %>% 
  mutate(timepoint = case_when(condition == "control_0h" ~ -20,
                               condition == "control_48h" ~ -12,
                               condition == "UPF1_Nter_0h" ~ 0,
                               condition == "UPF1_Nter_2h" ~ 2,
                               condition == "UPF1_Nter_4h" ~ 4,
                               condition == "UPF1_Nter_8h" ~ 8,
                               condition == "UPF1_Nter_12h" ~ 12,
                               condition == "UPF1_Nter_24h" ~ 24,
                               condition == "UPF1_Nter_48h" ~ 48,
                               condition == "control_R0h" ~ 84,
                               condition == "UPF1_Nter_24h_R0h" ~ 24,
                               condition == "UPF1_Nter_24h_R2h" ~ 26,
                               condition == "UPF1_Nter_24h_R4h" ~ 28,
                               condition == "UPF1_Nter_24h_R8h" ~ 32,
                               condition == "UPF1_Nter_24h_R12h" ~ 36,
                               condition == "UPF1_Nter_24h_R24h" ~ 48,
                               condition == "UPF1_Nter_24h_R48h" ~ 72)) %>% 
  mutate(exp_cond = case_when(condition == "control_0h" ~ "control",
                              condition == "control_48h" ~ "control",
                              condition == "UPF1_Nter_0h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_2h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_4h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_8h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_12h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_24h" ~ "UPF1_AID",
                              condition == "UPF1_Nter_48h" ~ "UPF1_AID",
                              condition == "control_R0h" ~ "control",
                              condition == "UPF1_Nter_24h_R0h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R2h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R4h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R8h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R12h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R24h" ~ "UPF1_AID_Rec",
                              condition == "UPF1_Nter_24h_R48h" ~ "UPF1_AID_Rec")) %>% 
  mutate(condition = fct_relevel(condition,
                                 "control_0h",
                                 "control_48h",
                                 "UPF1_Nter_0h",
                                 "UPF1_Nter_2h",
                                 "UPF1_Nter_4h",
                                 "UPF1_Nter_8h",
                                 "UPF1_Nter_12h",
                                 "UPF1_Nter_24h",
                                 "UPF1_Nter_48h",
                                 "control_R0h",
                                 "UPF1_Nter_24h_R0h",
                                 "UPF1_Nter_24h_R2h",
                                 "UPF1_Nter_24h_R4h",
                                 "UPF1_Nter_24h_R8h",
                                 "UPF1_Nter_24h_R12h",
                                 "UPF1_Nter_24h_R24h",
                                 "UPF1_Nter_24h_R48h"))

Rev_1_F2_B_gene_test_down_klus_cond_n <- Rev_1_F2_B_gene_test_down_klus_cond_forPlot %>% 
  group_by(DGE_cluster_down) %>% 
  distinct(gene_id,.keep_all = TRUE) %>% 
  summarize(n=n()) %>% 
  dplyr::rename("DGE_cluster" = "DGE_cluster_down")

#### Combined Up/Down ------------------------------------------------

Rev_1_F2_B_gene_test_combined_klus <- Rev_1_F2_B_gene_test_down_klus_cond_forPlot %>% 
  # filter(exp_cond %in% c("UPF1_AID",
  #                        "UPF1_AID_Rec")) %>% 
  filter(!condition %in% c("UPF1_Nter_48h",
                           "UPF1_Nter_24h_R0h")) %>% 
  dplyr::rename("DGE_cluster" = "DGE_cluster_down") %>% 
  bind_rows(Rev_1_F2_B_gene_test_up_klus_cond_forPlot %>% 
              # filter(exp_cond %in% c("UPF1_AID",
              #                        "UPF1_AID_Rec")) %>% 
              filter(!condition %in% c("UPF1_Nter_48h",
                                       "UPF1_Nter_24h_R0h")) %>% 
              dplyr::rename("DGE_cluster" = "DGE_cluster_up")) %>% 
  mutate(DGE_cluster = (fct_relevel(DGE_cluster,
                                    "up 1:early",
                                    "down 1:early",
                                    "up 2:delayed",
                                    "down 2:delayed",
                                    "up 3:late",
                                    "down 3:late",
                                    "up 4:inverse",
                                    "down 4:inverse"))) %>% 
  # group_by(timepoint,DGE_cluster_down) %>% 
  # summarize(median_exp = median(score)) %>% 
  ggplot(aes(x=timepoint,
             y=score)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_vline(xintercept = 24,
             color="gray15",
             linewidth=0.25,
             linetype = "dashed") +
  # geom_point() +
  # geom_boxplot(outlier.shape=NA) +
  ggrastr::rasterise(geom_line(
    # data=Rev_1_F2_B_gene_test_down_klus_cond_forPlot %>% 
    #                              filter(exp_cond %in% c("UPF1_AID",
    #                                                     "UPF1_AID_Rec")) %>%
    #                              filter(!condition %in% c("UPF1_Nter_48h")) %>% 
    #                              dplyr::rename("DGE_cluster" = "DGE_cluster_down") %>% 
    #                              bind_rows(Rev_1_F2_B_gene_test_up_klus_cond_forPlot %>% 
    #                                          filter(exp_cond %in% c("UPF1_AID",
    #                                                                 "UPF1_AID_Rec")) %>%
    #                                          filter(condition != "UPF1_Nter_48h") %>% 
    #                                          dplyr::rename("DGE_cluster" = "DGE_cluster_up")),
    aes(group=gene_id,
        color=DGE_cluster),
    alpha=0.5,
    linewidth=0.1),
    dpi=600, dev = "cairo") +
  stat_smooth(
    # data=Rev_1_F2_B_gene_test_down_klus_cond_forPlot %>% 
    #             filter(exp_cond %in% c("UPF1_AID",
    #                                    "UPF1_AID_Rec")) %>%
    #             filter(!condition %in% c("UPF1_Nter_48h")) %>% 
    #             dplyr::rename("DGE_cluster" = "DGE_cluster_down") %>% 
    #             bind_rows(Rev_1_F2_B_gene_test_up_klus_cond_forPlot %>% 
    #                         filter(exp_cond %in% c("UPF1_AID",
    #                                                "UPF1_AID_Rec")) %>%
    #                         filter(condition != "UPF1_Nter_48h") %>% 
    #                         dplyr::rename("DGE_cluster" = "DGE_cluster_up")),
    color="black",method = "gam", formula = y ~ s(x,k=7, bs="cs"), size = 1) +
  geom_boxplot(aes(group=interaction(timepoint,exp_cond),
                   y=score,
                   fill=exp_cond),
               outlier.shape = NA,
               fatten=2,
               linewidth=0.25) +
  geom_text(data=bind_rows(Rev_1_F2_B_gene_test_up_klus_cond_n,
                           Rev_1_F2_B_gene_test_down_klus_cond_n),
            aes(label=paste0("n=",n)),
            size = 5*0.36,
            x=36,
            y=3.5) +
  scale_x_continuous(breaks = c(-20,-12,0,12,24,36,48,60,72,84),
                     labels = c("0",
                                "48",
                                "0",
                                "12",
                                "24",
                                "36",
                                "48",
                                "60",
                                "72",
                                "0")) +
  scale_fill_manual(values=c("control" = "white",
                             "UPF1_AID_Rec" = "#F8CF69",
                             "UPF1_AID" = "#65BB8C")) +
  scale_color_manual(values=c("up 1:early" = "#67001f",
                              "up 2:delayed" = "#b2182b",
                              "up 3:late" = "#d6604d",
                              "up 4:inverse" = "#cfbeb4",
                              "down 4:inverse" = "#b9c3c8",
                              "down 3:late" = "#92c5de",
                              "down 2:delayed" = "#4393c3",
                              "down 1:early" = "#053061")) +
  guides(color="none",
         # fill="none"
  ) +
  labs(x = "Time (h)",
       y = "Z-score\n(log2 normalized counts)") +
  facet_wrap(~DGE_cluster,ncol=2) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(total_width = unit(65, "mm"),
                   total_height = unit(100, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_B_gene_test_combined_klus.pdf"),
       Rev_1_F2_B_gene_test_combined_klus,
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### GO cluster -----------------------------------------------------------------

# Obtain background of expressed genes from degradation/recovery HCT116 data
Rev_1_F2_GO_DGE_bg_full <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c(
    "HCT116_UPF1_AID_degradation_this_Study",
    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Make lists of DGE cluster - excluding the non-informative sets
Rev_1_F2_GO_DGE_NMD_cluster <- split(gtf_gencode_df_short_DGE_cluster %>% 
                                       filter(!DGE_cluster %in% c("expressed", "not_expressed")) %>% 
                                       separate(gene_id, c("gene_id", "Version")) %>% 
                                       pull(gene_id),
                                     gtf_gencode_df_short_DGE_cluster %>% 
                                       filter(!DGE_cluster %in% c("expressed", "not_expressed")) %>% 
                                       pull(DGE_cluster))

# Perform GO analysis
Rev_1_F2_GO_DGE_NMD_cluster_gostres <- gprofiler2::gost(query = Rev_1_F2_GO_DGE_NMD_cluster,
                                            custom_bg = Rev_1_F2_GO_DGE_bg_full,
                                            # multi_query = TRUE,
                                            sources = "GO:BP",
                                            domain_scope = "custom",
                                            organism = "hsapiens",
                                            correction_method = c("gSCS"),
                                            # evcodes = TRUE,
                                            # as_short_link = TRUE,
                                            significant = FALSE)

# Results as dataframe
Rev_1_F2_GO_DGE_NMD_cluster_gostres_result <- Rev_1_F2_GO_DGE_NMD_cluster_gostres$result

 # Plot significance of GO terms per Cluster
Rev_1_F2_GO_DGE_NMD_cluster_gostres_result %>% 
  mutate(query = fct_rev(fct_relevel(query,
                                     "up 1:early",
                                     "up 2:delayed",
                                     "up 3:late",
                                     "up 4:inverse",
                                     "down 4:inverse",
                                     "down 3:late",
                                     "down 2:delayed",
                                     "down 1:early"
  ))) %>% 
  # filter(p_value != 1) %>% 
  ggplot(aes(x=-log10(p_value),
             y=query)) +
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
  geom_vline(xintercept = -log10(0.05)) +
  guides(fill="none") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_B_DGE_cluster_GO_absolute.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###
## Rev_1 - Figure 2 - (C) ImpulseDE2 -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/ImpulseDE2/Rev_1_F2_C_ImpulseDE2.rds")

### Overarching analyses ----------------------------------------------------

# How many significant genes per biotype
Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  group_by(type, bestFit) %>% 
  dplyr::count() %>% 
  mutate(timecourse="combined") %>% 
  mutate(timecourse=fct_relevel(timecourse,
                                "combined")) %>% 
  ggplot(aes(y=type,
             x=n,
             fill=bestFit)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        strip.text=element_text(size=5), 
        axis.text=element_text(size=5, colour = 'black'), 
        axis.title=element_text(size=5), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5, colour = 'black'), 
        legend.key.size = unit(0.15, "cm"),
        legend.margin=margin(c(0.15,0.15,0.15,0.15)),
        plot.title = element_text(size = 5), 
        plot.subtitle = element_text(size = 5), 
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial")) + 
  geom_col(position="stack",
           color="black",
           linewidth=0.2) +
  scale_fill_manual(values=c("Impulse" = "#4D5F8E",
                             "Sigmoid" = "#C582B2",
                             "noFit" = "darkgray")) +
  force_panelsizes(rows = unit(5, "mm"),
                   cols = unit(30, "mm")) +
  labs(x="Number of significant genes",
       y="GENCODE\ngene biotype",
       fill="Best fit") +
  guides(fill=guide_legend(reverse=T))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_C_ImpulseDE2_biotypes_sig.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Frequency of onset/offset time

Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  mutate(timecourse=fct_rev(fct_relevel(timecourse,
                                        "combined"))) %>% 
  filter(bestFit != "noFit") %>% 
  dplyr::count(timecourse,bestFit)

# How many sig. genes per DGE_cluster
Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  # filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  mutate(timecourse=(fct_relevel(timecourse,
                                 "combined"))) %>% 
  mutate(bestFit = case_when(sigImpulseDE2 == "n.s." ~ "n.s.",
                             TRUE ~ bestFit)) %>% 
  filter(!DGE_cluster %in% c("not_expressed", "complex", "expressed")) %>% 
  dplyr::count(DGE_cluster, timecourse, bestFit) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  group_by(DGE_cluster, timecourse) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  mutate(bestFit = fct_relevel(bestFit,
                               "Sigmoid",
                               "Impulse",
                               "noFit",
                               "n.s.")) %>%
  ggplot(aes(y=DGE_cluster,
             x=n_per,
             fill=bestFit)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5), 
        legend.key.size = unit(0.15, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_col(position="stack",
           color="black",
           linewidth=0.1) +
  geom_text(aes(label = case_when(n_per > 50 ~ paste0(n,
                                                      " (",
                                                      round(n_per,0),
                                                      "%)"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(
    name = "Best fit", 
    values = c("Impulse" = "#4D5F8E",
               "Sigmoid" = "#C582B2",
               "noFit" = "darkgray",
               "n.s." = "white")
  ) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) +
  labs(x="Fraction of significant genes",
       y="Sig. regulated DGE cluster",
       fill="Best fit") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, reverse=T))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_C_ImpulseDE2_DGE_cluster_percent.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Boxplot padj per DGE cluster
Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  # filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(!DGE_cluster %in% c("not_expressed", "complex", "expressed")) %>% 
  filter(bestFit %in% c("Impulse", "Sigmoid")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-300)) %>% 
  ggplot(aes(x=-log10(padj),
             y=DGE_cluster,
             fill=bestFit)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5), 
        legend.key.size = unit(0.15, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="dashed") +
  geom_boxplot(outliers = FALSE,
               color="black",
               linewidth=0.2) +
  scale_fill_manual(values=c("Impulse" = "#4D5F8E",
                             "Sigmoid" = "#C582B2")) +
  facet_wrap(~timecourse) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_C_ImpulseDE2_DGE_cluster_boxplot_padj.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  


# Boxplot per DGE cluster

Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late",
                            "up 4:inverse",
                            "down 4:inverse"
  )) %>% 
  pivot_longer(cols=c(t1,t2),
               names_to = "Time",
               values_to = "Values") %>% 
  mutate(Time = fct_rev(fct_relevel(Time,
                                    "t1",
                                    "t2"))) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  filter(bestFit == "Impulse") %>% 
  ggplot(aes(x=Values,
             y=DGE_cluster,
             fill=Time)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
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
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="dashed") +
  geom_vline(xintercept=24,
             linewidth=0.25,
             linetype="dashed") +
  geom_boxplot(outliers = FALSE,
               color="black",
               linewidth=0.2) +
  scale_x_continuous(name = "time (h)",
                     breaks = c(-24,0,24,48,72),
                     labels = c(-24,0,24,48,72), 
                     limits = c(-26,74)) +
  scale_fill_manual(
    name = "Parameter", 
    values = c("t1" = "#00696D", 
               "t2" = "#FDF5CA"),
    labels = c("offset time t2", "onset time t1")
  ) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, reverse=T)) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_C_ImpulseDE2_DGE_cluster_boxplot_times.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

Rev_1_F2_C_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late",
                            "up 4:inverse",
                            "down 4:inverse")) %>% 
  pivot_longer(cols=c(t1,t2),
               names_to = "Time",
               values_to = "Values") %>% 
  filter(bestFit == "Impulse") %>% 
  group_by(DGE_cluster, Time) %>% 
  summarize(median = round(median(Values),1))

###
## Rev_1 - Figure 2 - (D) NMD relevance -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/NMD_relevance/DESeq2_DGE_combined_DGE_cluster_NMD_relevance.rds")

### Plot NMD relevance - barplot --------------------------------------------------------------------
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(DGE_cluster, UpDown, NMD_n_sig_perc) %>% 
  mutate(DGE_cluster = fct_rev(DGE_cluster)) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  # mutate(n = case_when(UpDown == "up" ~ n,
  #                      UpDown == "down" ~ -n)) %>% 
  ggplot(aes(y=DGE_cluster,
             x=n,
             fill=NMD_n_sig_perc)) +
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
        legend.margin=margin(c(0.1,0.1,0.1,0.1)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_col(position="fill") +
  # geom_vline(xintercept=0,
  #            linewidth=0.25,
  #            linetype="dashed") +
  scale_fill_gradientn(colours = c("#CACFD0",
                                   "#647588",
                                   # "#624B27",
                                   "#B09771",
                                   "#E5BF86")) +
  # scale_fill_viridis(direction = -1,
  #                    option="mako") +
  labs(fill="NMD\nrelevance\n(%)") +
  guides(size=guide_legend(nrow=2,byrow=TRUE),
         fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6))) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_DGE_cluster_NMD_confidence_perc_fill.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### NMD relevance bins ------------------------------------------------------

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late")) %>% 
  mutate(NMD_bin = fct_rev(NMD_bin)) %>% 
  dplyr::count(UpDown,NMD_bin) %>% 
  mutate(n=case_when(UpDown == "up" ~ n,
                     UpDown == "down" ~ -n)) %>% 
  mutate(NMD_bin = fct_rev(NMD_bin)) %>% 
  ggplot(aes(x=n,
             y=NMD_bin,
             fill=NMD_bin)) +
  theme_minimal() + 	
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
  geom_col(color="black",
           position=position_dodge(),
           linewidth=0.2) +
  scale_fill_manual(values = c("[0,25]" = "#CACFD0",
                               "(25,50]" =  "#647588",
                               "(50,75]" =  "#B09771",
                               "(75,100]" =  "#E5BF86")) +
  # scale_fill_viridis_d(direction = -1,
  #                      begin = 0.1,
  #                      end = 0.9,
  #                      option="mako") +
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="dashed") +
  geom_text(
    aes(label = case_when(n>0 ~ as.character(n),
                          TRUE ~ "")), 
    ## make labels left-aligned
    hjust = 0, nudge_x = +200,
    size = 5*0.36,
  ) +
  geom_text(
    aes(label = case_when(n>0 ~ "",
                          TRUE ~ as.character(-n))), 
    ## make labels left-aligned
    hjust = 1, nudge_x = -200,
    size = 5*0.36,
  ) +
  labs(fill="NMD\nrelevance\nbins") +
  scale_x_continuous(expand = expansion(c(0.25, 0.25))) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_DGE_cluster_NMD_relevance_bin_numbers.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

# Number per DGE_cluster
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  dplyr::count(NMD_bin, UpDown) %>% 
  group_by(UpDown) %>% 
  mutate(n_per = round(100*(n/sum(n)),1))


### Plot NMD relevance / BIN - barplot --------------------------------------------------------------------
DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(DGE_cluster, UpDown, NMD_bin) %>% 
  mutate(DGE_cluster = fct_rev(DGE_cluster)) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(NMD_bin = fct_rev(fct_relevel(NMD_bin,
                                   "up_(75,100]",
                                   "up_(50,75]",
                                   "up_(25,50]",
                                   "up_[0,25]",
                                   "down_[0,25]",
                                   "down_(25,50]",
                                   "down_(50,75]",
                                   "down_(75,100]",
  ))) %>% 
  # mutate(n = case_when(UpDown == "up" ~ n,
  #                      UpDown == "down" ~ -n)) %>% 
  ggplot(aes(y=DGE_cluster,
             x=n,
             fill=NMD_bin)) +
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
        legend.margin=margin(c(0.1,0.1,0.1,0.1)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_col(position="fill") +
  # geom_vline(xintercept=0,
  #            linewidth=0.25,
  #            linetype="dashed") +
  scale_fill_manual(values = c("[0,25]" = "#CACFD0",
                               "(25,50]" =  "#647588",
                               "(50,75]" =  "#B09771",
                               "(75,100]" =  "#E5BF86")) +
  # scale_fill_viridis_d(direction = 1,
  #                    option="mako") +
  labs(fill="NMD\nrelevance\n(%)") +
  guides(size=guide_legend(nrow=2,byrow=TRUE),
         fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6))) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_DGE_cluster_NMD_relevance_perc_fill_bin.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### GO bins -----------------------------------------------------------------

# Obtain background of expressed genes from degradation/recovery HCT116 data
Rev_1_F2_GO_DGE_bg_full <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c(
    "HCT116_UPF1_AID_degradation_this_Study",
    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Make lists of DGE cluster - excluding the non-informative sets
Rev_1_F2_GO_DGE_NMD_bin_UpDown <- split(DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
                                          separate(gene_id, c("gene_id", "Version")) %>% 
                                          pull(gene_id),
                                        DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
                                          mutate(UpDown_NMD_bin = paste0(UpDown,"_",NMD_bin), .after=NMD_bin) %>% 
                                          pull(UpDown_NMD_bin))

# Perform GO analysis
Rev_1_F2_GO_DGE_NMD_gostres <- gprofiler2::gost(query = Rev_1_F2_GO_DGE_NMD_bin_UpDown,
                                    custom_bg = Rev_1_F2_GO_DGE_bg_full,
                                    # multi_query = TRUE,
                                    sources = "GO:BP",
                                    domain_scope = "custom",
                                    organism = "hsapiens",
                                    correction_method = c("gSCS"),
                                    # evcodes = TRUE,
                                    # as_short_link = TRUE,
                                    significant = FALSE)

# Results as dataframe
Rev_1_F2_GO_DGE_NMD_gostres_result <- Rev_1_F2_GO_DGE_NMD_gostres$result

# Plot significance of GO terms per NMD relevance bin
Rev_1_F2_GO_DGE_NMD_gostres_result %>% 
  mutate(query = fct_rev(fct_relevel(query,
                                     "up_(75,100]",
                                     "up_(50,75]",
                                     "up_(25,50]",
                                     "up_[0,25]",
                                     "down_[0,25]",
                                     "down_(25,50]",
                                     "down_(50,75]",
                                     "down_(75,100]",
  ))) %>% 
  # filter(p_value != 1) %>% 
  ggplot(aes(x=-log10(p_value),
             y=query)) +
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
  geom_vline(xintercept = -log10(0.05)) +
  guides(fill="none") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_DGE_cluster_NMD_relevance_bin_GO_absolute.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### UPF1 log2FC bin ---------------------------------------------------------

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late")) %>% 
  dplyr::select(gene_id, NMD_n_sig_perc, NMD_bin, UpDown) %>% 
  left_join(DESeq2_DGE_combined %>% 
              filter(experimentSet %in% c("HeLa_UPF1_KD_cytoplasmic_Longman_2020",
                                          "HCT116_UPF1_AID_degradation_this_Study",
                                          "HCT116_UPF1_FKBP_degradation_this_Study",
                                          "HEK293_UPF1_FKBP_degradation_this_Study"
              )) %>% 
              filter(!condition_2 %in% c("control_48h",
                                         "UPF1_Nter_0h",
                                         "UPF1_Nter_2h",
                                         "UPF1_Nter_4h",
                                         "UPF1_Nter_8h",
                                         "UPF1_FKBP_HEK_0h",
                                         "UPF1_FKBP_HCT_0h")) %>% 
              mutate(condition_2 = case_when(condition_2 == "UPF1" ~ publicationName,
                                             condition_2 == "SMG1i" ~ experimentSet,
                                             TRUE ~ condition_2)) %>% 
              mutate(experiment_set = fct_relevel(experimentSet,
                                                  "HeLa_UPF1_KD_cytoplasmic_Longman_2020"))) %>% 
  arrange(experiment_set) %>% 
  filter(!is.na(condition_2)) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  # group_by(condition_2, bin) %>% 
  # summarise_at(vars(log2FoldChange),
  #              list(Q1=~quantile(., probs = 0.25, na.rm = TRUE),
  #                   median=median,
  #                   Q3=~quantile(., probs = 0.75, na.rm = TRUE)), na.rm = TRUE) %>% 
  ggplot(aes(x=log2FoldChange,
             y=condition_2,
             fill=NMD_bin)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
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
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="dashed") +
  geom_vline(xintercept=1,
             linewidth=0.25,
             color="darkred",
             linetype="dashed") +
  # geom_pointrange(aes(x=median,
  #                     xmin=Q1,
  #                     xmax=Q3),
  #                 position = position_dodge(width = 0.95),
  #                 linewidth = 0.25,
  #                 size = 0.4,
  #                 stroke=0.1,
  #                 alpha=0.75,
  #                 shape=21) +
  geom_boxplot(outliers =  FALSE,
               # coef = 0,
               color="black",
               linewidth=0.1) +
  scale_fill_manual(values = c("[0,25]" = "#CACFD0",
                               "(25,50]" =  "#647588",
                               "(50,75]" =  "#B09771",
                               "(75,100]" =  "#E5BF86")) +
  # scale_fill_viridis_d(direction = -1,
  #                      begin = 0.1,
  #                      end = 0.9,
  #                      option="mako") +
  scale_x_continuous(name = "log2FC",
                     breaks = c(0,2,4),
                     labels = c(0,2,4),
                     limits = c(-0.5,4)) +
  force_panelsizes(rows = unit(6*8+0.4, "mm"),
                   cols = unit(20+0.4, "mm")) 


ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_DGE_cluster_NMD_relevance_bin_log2FC_UPF1.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### UPF1 log2FC bin heatmap ---------------------------------------------------------

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  filter(DGE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late")) %>% 
  dplyr::select(gene_id,NMD_n_sig_perc,NMD_bin, UpDown) %>% 
  left_join(DESeq2_DGE_combined %>% 
              filter(experimentSet %in% c("HeLa_UPF1_KD_cytoplasmic_Longman_2020",
                                          "HEK293_SMG567_KO_KD_Boehm_2021",
                                          "HEK293_SMG67_KD_Britto_Borges_2024",
                                          "HeLa_SMG67_KD_Britto_Borges_2024",
                                          "MCF7_SMG67_KD_Britto_Borges_2024",
                                          "U2OS_SMG67_KD_Britto_Borges_2024",
                                          "HEK293_UPF3_dKO_Wallmeroth_2022",
                                          "HEK293_CASC3_KO_Gerbracht_2020",
                                          "HCT116_UPF1_AID_degradation_this_Study",
                                          "HCT116_UPF1_AID_recovery_this_Study",
                                          "HCT116_UPF1_FKBP_degradation_this_Study",
                                          "HEK293_UPF1_FKBP_degradation_this_Study",
                                          "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                                          "HFF_SMG1i_this_Study",
                                          "HUVEC_SMG1i_this_Study"
              )) %>% 
              filter(!condition_2 %in% c("SMG8_delKID_0uM",
                                         "SMG8_delKID_01uM",
                                         "SMG8_delKID_1uM")) %>% 
              mutate(condition_2 = case_when(condition_2 == "UPF1" ~ publicationName,
                                             condition_2 == "SMG1i" ~ experimentSet,
                                             TRUE ~ condition_2)) %>% 
              mutate(experiment_set = fct_relevel(experimentSet,
                                                  "HEK293_SMG567_KO_KD_Boehm_2021",
                                                  "HEK293_SMG67_KD_Britto_Borges_2024",
                                                  "HeLa_SMG67_KD_Britto_Borges_2024",
                                                  "MCF7_SMG67_KD_Britto_Borges_2024",
                                                  "U2OS_SMG67_KD_Britto_Borges_2024",
                                                  "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                                                  "HFF_SMG1i_this_Study",
                                                  "HUVEC_SMG1i_this_Study",
                                                  "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
                                                  "HEK293_UPF3_dKO_Wallmeroth_2022",
                                                  "HEK293_CASC3_KO_Gerbracht_2020",
                                                  "HCT116_UPF1_AID_degradation_this_Study",
                                                  "HCT116_UPF1_AID_recovery_this_Study",
                                                  "HCT116_UPF1_FKBP_degradation_this_Study",
                                                  "HEK293_UPF1_FKBP_degradation_this_Study"))) %>% 
  arrange(experiment_set) %>% 
  filter(!is.na(condition_2)) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>% 
  mutate(condition_2 = (fct_inorder(as_factor(condition_2)))) %>% 
  # distinct(condition_2)
  group_by(condition_2, NMD_bin) %>%
  summarise_at(vars(log2FoldChange),
               list(Q1=~quantile(., probs = 0.25, na.rm = TRUE),
                    median=median,
                    Q3=~quantile(., probs = 0.75, na.rm = TRUE)), na.rm = TRUE) %>%
  ggplot(aes(y=NMD_bin,
             x=condition_2,
             fill=median)) +
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
  geom_tile(color="black",
            linewidth=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey90") +
  force_panelsizes(rows = unit(8+0.4, "mm"),
                   cols = unit(50*2+0.4, "mm"))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y="",
       x="",
       fill="median\nDGE\nlog2FC") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.vjust = 0.5))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_D_median_log2FC_perBin.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

## Rev_1 - Figure 2 - (E) Colombo meta -----------------------------------------------------------
# Idea: check how "significant" the DGE cluster genes are in the Colombo et al 2017 meta analysis (PMID: 27864472)

### Colombo 2017 data -------------------------------------------------------
# Import Excel file (slightly modified to ease import into R)
Colombo_Supplemental_Table_S2_forR <- readxl::read_excel("Resources/External/Colombo2017_PMID_27864472/Supplemental_Table_S2_forR.xlsx")

# Join with GENCODE and determine significance (lower cutoffs due to single-end, etc.)
Colombo_Supplemental_Table_S2_forR_annotated <- Colombo_Supplemental_Table_S2_forR %>% 
  left_join(gtf_gencode_df_short %>% 
              filter(type == "gene") %>% 
              separate(gene_id,
                       into=c("gene_id_plain",NA),
                       remove = FALSE),
            by=c("gene" = "gene_id_plain")) %>% 
  mutate(sign_meta_meta = case_when(meta_meta < 0.01 ~ "sig.",
                                    TRUE ~ "Not Sig."),
         sign_SMGs = case_when(meta_SMGs < 0.01 ~ "sig.",
                               TRUE ~ "Not Sig."))

#### Plot --------------------------------------------------------------------

# Plot -log10(padj) per DGE_cluster
gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin %>% 
  left_join(Colombo_Supplemental_Table_S2_forR_annotated %>% dplyr::select(gene_id, UPF1_FDR, dKD_SMG6_FDR, dKD_SMG7_FDR, meta_SMG6, meta_SMG7, meta_SMGs, meta_meta, sign_meta_meta)) %>% 
  pivot_longer(cols=c(UPF1_FDR,  meta_SMGs, meta_meta),
               names_to = "meta_analyses",
               values_to = "pval") %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  mutate(meta_analyses = fct_relevel(meta_analyses,
                                     "meta_meta",
                                     "meta_SMGs",
                                     "UPF1_FDR")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  ggplot(aes(y=-log10(pval),
             x=DGE_cluster,
             fill=meta_analyses)) +
  geom_boxplot(outliers = FALSE,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.2) + 
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major.y = element_line(colour = 'gray80', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
  scale_fill_manual(values=c("meta_meta" = "#B8321A",
                             "meta_SMGs" = "#CBAC7C",
                             "meta_SMG6" = "#9ebcda",
                             "UPF1_FDR" = "#66B2B2")) +
  # scale_fill_manual(values=c("fast_up" = "#67001f",
  #                            "medium_up" = "#b2182b",
  #                            "slow_up" = "#d6604d",
  #                            "peak_fast_up" = "#f4a582",
  #                            "peak_medium_up" = "#fddbc7",
  #                            "expressed" = "gray30",
  #                            "complex" = "#6A3F6C",
  #                            "peak_medium_down" = "#d1e5f0",
  #                            "peak_fast_down" = "#92c5de",
  #                            "slow_down" = "#4393c3",
  #                            "medium_down" = "#2166ac",
#                            "fast_down" = "#053061")) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # coord_cartesian(xlim = c(0,50)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(42, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_E_Colombo_meta_DGE_cluster.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

# Plot -log10(padj) per NMD_bin
gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin %>% 
  left_join(Colombo_Supplemental_Table_S2_forR_annotated %>% dplyr::select(gene_id, UPF1_FDR, dKD_SMG6_FDR, dKD_SMG7_FDR, meta_SMG6, meta_SMG7, meta_SMGs, meta_meta, sign_meta_meta)) %>% 
  pivot_longer(cols=c(UPF1_FDR,  meta_SMGs, meta_meta),
               names_to = "meta_analyses",
               values_to = "pval") %>% 
  # filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(NMD_bin)) %>% 
  mutate(meta_analyses = fct_relevel(meta_analyses,
                                     "meta_meta",
                                     "meta_SMGs",
                                     "UPF1_FDR")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  mutate(NMD_bin_updown = paste0(UpDown,"_",NMD_bin)) %>% 
  mutate(NMD_bin_updown = fct_rev(fct_relevel(NMD_bin_updown,
                                              "up_(75,100]",
                                              "up_(50,75]",
                                              "up_(25,50]",
                                              "up_[0,25]",
                                              "n.s._n.s.",
                                              "down_[0,25]",
                                              "down_(25,50]",
                                              "down_(50,75]",
                                              "down_(75,100]",
  ))) %>% 
  filter(NMD_bin_updown %in% c("up_(75,100]",
                               "up_(50,75]",
                               "up_(25,50]",
                               "up_[0,25]",
                               "n.s._n.s.")) %>% 
  ggplot(aes(y=-log10(pval),
             x=NMD_bin_updown,
             fill=meta_analyses)) +
  geom_boxplot(outliers = FALSE,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.2) + 
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major.y = element_line(colour = 'gray80', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
  scale_fill_manual(values=c("meta_meta" = "#B8321A",
                             "meta_SMGs" = "#CBAC7C",
                             "meta_SMG6" = "#9ebcda",
                             "UPF1_FDR" = "#66B2B2")) +
  # scale_fill_manual(values=c("fast_up" = "#67001f",
  #                            "medium_up" = "#b2182b",
  #                            "slow_up" = "#d6604d",
  #                            "peak_fast_up" = "#f4a582",
  #                            "peak_medium_up" = "#fddbc7",
  #                            "expressed" = "gray30",
  #                            "complex" = "#6A3F6C",
  #                            "peak_medium_down" = "#d1e5f0",
  #                            "peak_fast_down" = "#92c5de",
  #                            "slow_down" = "#4393c3",
  #                            "medium_down" = "#2166ac",
#                            "fast_down" = "#053061")) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # coord_cartesian(xlim = c(0,50)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(30, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_E_Colombo_meta_NMD_bin.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  


## Rev_1 - Figure 2 - (F) phospho-UPF1 Kurosaki 2014 -----------------------------------------------------

# Load essential data for these analyses
load("Resources/External/Kurosaki2014_PMID_25184677/Kurosaki2014_phospho_UPF1.rds")

##
### Stats --------------------------------------------------------------------
##

# How many significant? lower significance due to duplicates, etc.
Kurosaki2014_tbl %>% 
  dplyr::count(padj < 0.01, condition_2)

##
### Plot --------------------------------------------------------------------
##


#### Boxplot -----------------------------------------------------------------

Kurosaki2014_tbl %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  mutate(DGE_cluster = fct_rev(DGE_cluster)) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  ggplot(aes(y=log2FoldChange,
             x=DGE_cluster,
             fill=DGE_cluster)) +
  geom_hline(yintercept = 0,
             linewidth=0.25,
             color="black") +
  geom_boxplot(outliers = FALSE,
               fatten=2,
               linewidth=0.25) +
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                             "up 2:delayed" = "#b2182b",
                             "up 3:late" = "#d6604d",
                             "down 3:late" = "#92c5de",
                             "down 2:delayed" = "#4393c3",
                             "down 1:early" = "#053061",
                             "expressed" = "gray30")) +
  facet_wrap(~condition_2) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(28, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_F_Kurosaki_pUPF1_DGE_cluster.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  


#### Boxplot NMD bin -----------------------------------------------------------------

Kurosaki2014_tbl %>% 
  # filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(DGE_cluster %in% c("up 1:early", "up 2:delayed", "up 3:late", "expressed")) %>% 
  # filter(!is.na(DGE_cluster)) %>% 
  mutate(NMD_bin_updown = paste0(UpDown,"_",NMD_bin)) %>% 
  mutate(NMD_bin_updown = fct_rev(fct_relevel(NMD_bin_updown,
                                              "up_(75,100]",
                                              "up_(50,75]",
                                              "up_(25,50]",
                                              "up_[0,25]",
                                              "n.s._n.s.",
                                              "down_[0,25]",
                                              "down_(25,50]",
                                              "down_(50,75]",
                                              "down_(75,100]",
  ))) %>% 
  filter(NMD_bin_updown %in% c("up_(75,100]",
                               "up_(50,75]",
                               "up_(25,50]",
                               "up_[0,25]",
                               "n.s._n.s.")) %>% 
  ggplot(aes(y=log2FoldChange,
             x=NMD_bin_updown,
             fill=NMD_bin)) +
  geom_hline(yintercept = 0,
             linewidth=0.25,
             color="black") +
  geom_boxplot(outliers = FALSE,
               fatten=2,
               linewidth=0.25) +
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("n.s." = "white",
                               "[0,25]" = "#CACFD0",
                               "(25,50]" =  "#647588",
                               "(50,75]" =  "#B09771",
                               "(75,100]" =  "#E5BF86")) +
  facet_wrap(~condition_2, nrow=1) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(16, "mm"))

ggsave(file.path("Plots", "Figure2", "Rev_1_F2_F_Kurosaki_pUPF1_NMD_bin.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  
