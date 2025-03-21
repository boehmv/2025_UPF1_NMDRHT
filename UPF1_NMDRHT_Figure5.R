#!/usr/bin/env Rscript

# Title: UPF1_Revision_Figure5
# Objective: Code for generating Panels of Figure 5 + Figure S5 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_Revision_Analysis"

###
## Rev_1 - Figure 5 - (A1) edgeR DTE SRSF2 -----------------------------------------------------------
###

edgeR_DTE_NMDRHT_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  filter(gene_name == "SRSF2") %>% 
  dplyr::rename("transcript_type" = "transcript_biotype") %>% 
  mutate(annotation = "NMDHRT") %>% 
  bind_rows(edgeR_DTE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              filter(gene_name == "SRSF2") %>% 
              mutate(annotation = "GENCODE")) %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "protein_coding_CDS_not_defined", "lncRNA", "nonsense_mediated_decay", "retained_intron") ~ "other",
                          transcript_type == "protein_coding" ~ "coding",
                          transcript_type == "protein_coding_CDS_not_defined" ~ "coding_noCDS",
                          transcript_type == "lncRNA" ~ "lncRNA",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "retained_intron" ~ "RI")) %>% 
  ggplot(aes(x=logCPM,
             y=logFC,
             size=Overdispersion)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_line(color = "gray80",
                                        linewidth = 0.1,
                                        linetype = 1),
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
  geom_point(aes(fill=type),
             shape=21,
             stroke=0.2,
             color="black") +
  scale_size(range = c(1,2.5)) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "NMD" = "#E5BF86",
                             "RI" = "#B9D5EB")) +
  coord_cartesian(clip = "off") +
  ggrepel::geom_text_repel(aes(label=transcript_name),
                           size = 4*0.36,
                           # fill="lightgray",
                           max.overlaps = Inf,
                           # label.padding = unit(0.1, "lines"),
                           # label.size = 0,
                           parse=F,
                           xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  facet_wrap(~annotation) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_A1_NMDHRT_GENCODE_SRSF2_edgeR.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### SRSF2 12h N-AID-UPF1 Overdispersion -------------------------------------

edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  filter(gene_name == "SRSF2") %>% 
  # dplyr::rename("transcript_type" = "transcript_biotype") %>% 
  mutate(annotation = "NMDHRT") %>% 
  summarize(n=sum(Overdispersion))

edgeR_DTE_NMDRHT_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  filter(gene_name == "SRSF2") %>% 
  dplyr::rename("transcript_type" = "transcript_biotype") %>% 
  mutate(annotation = "NMDHRT") %>% 
  summarize(n=sum(Overdispersion))

###
## Rev_1 - Figure 5 - (A2) Salmon mapping -----------------------------------------------------------
###

### Load GENCODE-based Salmon QC data ---------------------------------------------------------------

# Load essential data
GENCODE_UPF1_Salmon_QC_forComparison <- read_csv("Resources/QC/GENCODE_UPF1_Salmon_QC_forComparison.csv")

# Calculate mean total counts and percent mapped per experimentSet
GENCODE_UPF1_Salmon_QC_forComparison_summary <- GENCODE_UPF1_Salmon_QC_forComparison %>% 
  group_by(publicationName) %>% 
  summarize(reads = mean(general.num_processed),
            percent_mapped = mean(general.percent_mapped),
            eq_classes = mean(general.num_eq_classes),
            libType = paste(unique(general.library_types), collapse = ', ')) %>% 
  mutate(stranded = case_when(libType %in% c("ISR", "ISF") ~ "S",
                              libType %in% c("IU") ~ "U",
                              TRUE ~ "other")) %>% 
  filter(!is.na(publicationName))

### Load NMDRHT-based Salmon QC data ---------------------------------------------------------------

# Load essential data
NMDRHT_UPF1_Salmon_QC_forComparison <- read_csv("Resources/QC/NMDRHT_UPF1_Salmon_QC_forComparison.csv")

# Calculate mean total counts and percent mapped per experimentSet
NMDRHT_UPF1_Salmon_QC_forComparison_summary <- NMDRHT_UPF1_Salmon_QC_forComparison %>% 
  group_by(publicationName) %>% 
  summarize(reads = mean(general.num_processed),
            percent_mapped = mean(general.percent_mapped),
            eq_classes = mean(general.num_eq_classes),
            libType = paste(unique(general.library_types), collapse = ', ')) %>% 
  mutate(stranded = case_when(libType %in% c("ISR", "ISF") ~ "S",
                              libType %in% c("IU") ~ "U",
                              TRUE ~ "other")) %>% 
  filter(!is.na(publicationName))

#### Combined plots ----------------------------------------------------------
# Percent mapped
GENCODE_UPF1_Salmon_QC_forComparison_summary %>% 
  mutate(annotation = "GENCODE") %>% 
  bind_rows(NMDRHT_UPF1_Salmon_QC_forComparison_summary %>% 
              mutate(annotation = "NMDRHT")) %>% 
  pivot_longer(cols=c(eq_classes, percent_mapped),
               names_to = "stat",
               values_to = "value") %>% 
  mutate(stat = fct_rev(stat)) %>% 
  mutate(publicationName = fct_rev(publicationName)) %>% 
  mutate(annotation = fct_rev(fct_relevel(annotation,
                                          "GENCODE",
                                          "NMDRHT"))) %>% 
  filter(stat == "percent_mapped") %>% 
  ggplot(aes(x=value,
             y=publicationName,
             fill=annotation)) +
  geom_col(position = position_dodge(),
           color="black",
           linewidth=0.25) +
  geom_text(aes(label=paste(round(value,0), "%")),
            hjust = 1.25,
            size = 5*0.36,
            color="white",
            position = position_dodge(width = 0.9)
  ) +
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'black', linewidth = 0.1),
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
  scale_fill_manual(values=c("GENCODE" = "#B80C66",
                             "NMDRHT" = "#B09771")) +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_A2_Annotation_overdispersion_gene_wise_percMapped.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (A3) Eq classes -----------------------------------------------------------
###

GENCODE_UPF1_Salmon_QC_forComparison_summary %>% 
  mutate(annotation = "GENCODE") %>% 
  bind_rows(NMDRHT_UPF1_Salmon_QC_forComparison_summary %>% 
              mutate(annotation = "NMDRHT")) %>% 
  pivot_longer(cols=c(eq_classes, percent_mapped),
               names_to = "stat",
               values_to = "value") %>% 
  mutate(stat = fct_rev(stat)) %>% 
  mutate(publicationName = fct_rev(publicationName)) %>% 
  mutate(annotation = fct_rev(fct_relevel(annotation,
                                          "GENCODE",
                                          "NMDRHT"))) %>% 
  filter(stat == "eq_classes") %>% 
  ggplot(aes(x=value,
             y=publicationName,
             fill=annotation)) +
  geom_col(position = position_dodge(),
           color="black",
           linewidth=0.25) +
  geom_text(aes(label=paste(round(value/10^5,0), " x 10^5")),
            hjust = 1.25,
            size = 5*0.36,
            color="white",
            position = position_dodge(width = 0.9)
  ) +
  theme_minimal() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        strip.text=element_text(size=6), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'black', linewidth = 0.1),
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
  scale_fill_manual(values=c("GENCODE" = "#B80C66",
                             "NMDRHT" = "#B09771")) +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_A3_Annotation_overdispersion_gene_wise_eqClasses.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


###
## Rev_1 - Figure 5 - (A4) edgeR gene-wise overdispersion -----------------------------------------------------------
###

edgeR_DTE_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                              "HCT116_UPF1_AID_recovery_this_Study",
                              "HCT116_UPF1_AID_degradation_riboMinus_this_Study",
                              "HCT116_UPF1_FKBP_degradation_this_Study",
                              "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
  group_by(gene_id, publicationName) %>% 
  mutate(cumul_overdisp = sum(Overdispersion)) %>% 
  ungroup() %>% 
  group_by(publicationName) %>% 
  mutate(cumul_overdisp = cumul_overdisp/length(unique(condition_2))) %>% 
  ungroup() %>% 
  mutate(annotation = "GENCODE") %>% 
  bind_rows(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                          "HCT116_UPF1_AID_recovery_this_Study",
                                          "HCT116_UPF1_AID_degradation_riboMinus_this_Study",
                                          "HCT116_UPF1_FKBP_degradation_this_Study",
                                          "HEK293_UPF1_FKBP_degradation_this_Study")) %>% 
              group_by(gene_id, publicationName) %>% 
              mutate(cumul_overdisp = sum(Overdispersion)) %>% 
              ungroup() %>% 
              group_by(publicationName) %>% 
              mutate(cumul_overdisp = cumul_overdisp/length(unique(condition_2))) %>% 
              ungroup() %>% 
              mutate(annotation = "NMDRHT")) %>% 
  mutate(publicationName = fct_rev(publicationName)) %>% 
  mutate(annotation = fct_rev(fct_relevel(annotation,
                                          "GENCODE",
                                          "NMDRHT"))) %>% 
  ggplot(aes(x=cumul_overdisp,
             y=publicationName,
             fill=annotation)) +
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
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'black', linewidth = 0.1),
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
  scale_fill_manual(values=c("GENCODE" = "#B80C66",
                             "NMDRHT" = "#B09771")) +
  guides(fill=guide_legend(reverse=T)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_A4_Annotation_overdispersion_gene_wise.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (B) edgeR DTE NMDRHT -----------------------------------------------------------
###

# Prepare data
edgeR_DTE_NMDRHT_combined_forPlot <- edgeR_DTE_NMDRHT_combined %>% 
  left_join(NMDRHT.v1.2_tbl_length_GC_mfe %>% dplyr::select(transcript_id, 
                                                            NMD_tx_status,
                                                            NMD_tx_reason)) %>% 
  filter(experimentSet %in% c(
    # "HEK293_UPF1_KD_Kishor_2019",
    "HeLa_UPF1_KD_cytoplasmic_Longman_2020",
    # "HeLa_UPF1_KD_Longman_2020",
    # "HEK293_UPF1_KD_Fritz_2022",
    # "K562_UPF1_KD_Hug_2022",
    # "HUH7_UPF1_KD_Lee_2022",
    # "HepG2_UPF1_KD_He_2023",
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
  filter(!condition_2 %in% c("SMG5_KD",
                             "SMG7_KO_2",
                             "SMG7_KO_34",
                             "SMG8_KO_0uM",
                             "SMG9_KO_0uM",
                             "control_01uM",
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
                                      "HCT116_UPF1_AID_degradation_this_Study",
                                      "HCT116_UPF1_AID_recovery_this_Study",
                                      "HCT116_UPF1_FKBP_degradation_this_Study",
                                      "HEK293_UPF1_FKBP_degradation_this_Study",
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
  mutate(significant = case_when(FDR < 0.0001 & abs(logFC) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(logFC>0 ~ "up",
                             logFC<0 ~ "down")) %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_transcripts = length(unique(transcript_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, NMD_tx_status) %>% 
  mutate(detected_transcripts_type = length(unique(transcript_id))) %>% 
  ungroup()

### Heatmap -----------------------------------------------------------------

edgeR_DTE_NMDRHT_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, NMD_tx_status, up_down, detected_transcripts_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(n_perTotalType = n/detected_transcripts_type) %>% 
  mutate(type_up_down = paste0(NMD_tx_status, " - ", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding - down",
                                     "coding - up",
                                     "predicted_coding - down",
                                     "predicted_coding - up",
                                     "NMD - down",
                                     "NMD - up",
                                     "predicted_NMD - down",
                                     "predicted_NMD - up",
                                     "mixed - down",
                                     "mixed - up",
                                     "lncRNA - down",
                                     "lncRNA - up"))) %>% 
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
  force_panelsizes(rows = unit(length(unique(edgeR_DTE_NMDRHT_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(24+0.6, "mm"))  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y="",
       x="",
       fill="Fraction of sig.\nregulated transcripts\nper NMD status") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.vjust = 0.5))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_B1_edgeR_DTE_Heatmap_perStatus.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Barplot -----------------------------------------------------------

edgeR_DTE_NMDRHT_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, NMD_tx_status, up_down, detected_transcripts_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(type_up_down = paste0(NMD_tx_status, " - ", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding - up",
                                     "coding - down",
                                     "predicted_coding - up",
                                     "predicted_coding - down",
                                     "NMD - up",
                                     "NMD - down",
                                     "predicted_NMD - up",
                                     "predicted_NMD - down",
                                     "mixed - up",
                                     "mixed - down",
                                     "lncRNA - up",
                                     "lncRNA - down"))) %>% 
  mutate(NMD_tx_status = (fct_relevel(NMD_tx_status,
                                      "coding",
                                      "predicted_coding",
                                      "NMD",
                                      "predicted_NMD",
                                      "mixed",
                                      "lncRNA"
  ))) %>% 
  mutate(n = case_when(up_down == "down" ~ -n,
                       up_down == "up" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n,
             fill=NMD_tx_status)) +
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
             color="darkgray",
             linewidth=0.25) +
  geom_col(color="black",
           linewidth = 0.2) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD"="#E5BF86",
                             "mixed" = "#624B27",
                             "predicted_NMD" = "#B09771",
                             "predicted_coding" = "#287DAB")) +
  # scale_fill_manual(values=c("coding" = "#7C7189",
  #                            "lncRNA" = "#BC8E7D",
  #                            "NMD"="#D04E59",
  #                            "mixed" = "#85BEDC",
  #                            "predicted_NMD" = "#F08683",
  #                            "predicted_coding" = "#CABEE9")) +
  labs(y="",
       x="Number of \nsig. regulated transcripts",
       fill="NMD status",
       size="Sig. DTE events") +
  force_panelsizes(rows = unit(length(unique(edgeR_DTE_NMDRHT_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(20, "mm")) +
  guides(fill=guide_legend(nrow=3,byrow=FALSE,
                           title.position = "top"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_B2_edgeR_DTE_Barplot_perStatus.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Examples -----------------------------------------------------------

# Define representative NMD-targeted genes plus UPF1
standard_heatmap_transcripts <- c("SRSF2-204_FSM",
                                  "GABARAPL1-208_FSM",
                                  "RBM3-204_FSM",
                                  "CCNB1IP1-214_FSM",
                                  "SNHG12-207_FSM")

edgeR_DTE_NMDRHT_combined_forPlot %>%  
  filter(transcript_name %in% standard_heatmap_transcripts) %>% 
  dplyr::select(condition_2, experimentSet, logFC, FDR, transcript_name) %>% 
  mutate(transcript_name = fct_relevel(transcript_name,
                                       "GABARAPL1-208_FSM",
                                       "RBM3-204_FSM",
                                       "SRSF2-204_FSM",
                                       "CCNB1IP1-214_FSM",
                                       "SNHG12-207_FSM"
  )) %>% 
  mutate(FDR = replace(FDR, FDR == 0, 1e-320)) %>% 
  ggplot(aes(x=transcript_name,
             y=condition_2,
             fill=logFC
  )) +
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
       fill="DTE log2FC") +
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
  force_panelsizes(rows = unit(length(unique(edgeR_DTE_NMDRHT_combined_forPlot$condition_2))*2+0.6, "mm"),
                   cols = unit(5*2+0.6, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_B3_edgeR_DTE_Examples.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (C) NMDRHT DTE cluster -----------------------------------------------------------
###

load("Resources/Cluster/Rev_1_F5_transcript_cluster_information.rds")

### Plot Heatmaps -----------------------------------------------------------

#### Up ----------------------------------------------------------------------
Rev_1_F5_transcript_hclust_split_up <- factor(Rev_1_F5_transcript_hclust_up_transcript_cluster_fix, levels=c("up 1:early",
                                                                                                             "up 2:delayed",
                                                                                                             "up 3:late",
                                                                                                             "up 4:inverse"))

Rev_1_F5_transcript_split_col = rep(1:17, each = 3)

ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 0,
                             col = 0),
                   labels = unique(Rev_1_F5_coldata$condition),
                   labels_rot = 90,
                   labels_just = "right",
                   labels_offset = unit(1, "npc"),
                   labels_gp = gpar(fontsize = 6),
                   height = unit(1, "cm"),
                   which = c("column"))
)

# Color setting for Heatmap
Rev_1_F5_transcript_col_up = colorRamp2(c(-6,-4,-2, -1, 0, 1, 2, 4, 6), rev(c("#88002D","#D4332A", "#F9834C", "#FFC08E", "white", "#B7C2DE", "#7797C9", "#0C6DA8", "#0E3F5C")))

# log2-transformed *Z-SCALED* heatmap of upregulated DGE transcripts
Rev_1_F5_transcript_hclust_up_heat <- Heatmap(matrix = Rev_1_F5_transcript_test_up_fix,
                                              col =  Rev_1_F5_transcript_col_up,
                                              cluster_rows = FALSE,
                                              cluster_columns = FALSE,
                                              name="Z-score\n(log2-normalized\ncounts)",
                                              row_split=Rev_1_F5_transcript_hclust_split_up, 
                                              column_split = Rev_1_F5_transcript_split_col,
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
                                              width = unit((length(unique(Rev_1_F5_coldata$condition))*2)+(-1+length(unique(Rev_1_F5_coldata$condition))*1), "mm"),
                                              height = unit(6, "cm"),
                                              column_title = "NMDRHT - edgeR DTE sig. upregulated")

Rev_1_F5_transcript_hclust_up_heat = draw(Rev_1_F5_transcript_hclust_up_heat)

pdf(file = file.path("Plots", "Figure5", "Rev_1_F5_C_only_NMDRHT_edgeR_DTE_sig_upregulated_hclust_heatmap.pdf"),
    useDingbats = FALSE,
    family="Arial",
    fonts="Arial")

draw(Rev_1_F5_transcript_hclust_up_heat)

dev.off()

#### Down ----------------------------------------------------------------------
Rev_1_F5_transcript_hclust_split_down <- factor(Rev_1_F5_transcript_hclust_down_transcript_cluster_fix, levels=c("down 4:inverse",
                                                                                                                 "down 3:late",
                                                                                                                 "down 2:delayed",
                                                                                                                 "down 1:early"))

Rev_1_F5_transcript_split_col = rep(1:17, each = 3)

ha = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 0,
                             col = 0),
                   labels = unique(Rev_1_F5_coldata$condition),
                   labels_rot = 90,
                   labels_just = "right",
                   labels_offset = unit(1, "npc"),
                   labels_gp = gpar(fontsize = 6),
                   height = unit(1, "cm"),
                   which = c("column"))
)

# Color setting for Heatmap
Rev_1_F5_transcript_col_down = colorRamp2(c(-6,-4,-2, -1, 0, 1, 2, 4, 6), rev(c("#88002D","#D4332A", "#F9834C", "#FFC08E", "white", "#B7C2DE", "#7797C9", "#0C6DA8", "#0E3F5C")))

# log2-transformed *Z-SCALED* heatmap of downregulated DGE transcripts
Rev_1_F5_transcript_hclust_down_heat <- Heatmap(matrix = Rev_1_F5_transcript_test_down_fix,
                                                col =  Rev_1_F5_transcript_col_down,
                                                cluster_rows = FALSE,
                                                cluster_columns = FALSE,
                                                name="Z-score\n(log2-normalized\ncounts)",
                                                row_split=Rev_1_F5_transcript_hclust_split_down, 
                                                column_split = Rev_1_F5_transcript_split_col,
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
                                                width = unit((length(unique(Rev_1_F5_coldata$condition))*2)+(-1+length(unique(Rev_1_F5_coldata$condition))*1), "mm"),
                                                height = unit(6, "cm"),
                                                column_title = "NMDRHT - edgeR DTE sig. downregulated")

Rev_1_F5_transcript_hclust_down_heat = draw(Rev_1_F5_transcript_hclust_down_heat)

pdf(file = file.path("Plots", "Figure5", "Rev_1_F5_C_only_NMDRHT_edgeR_DTE_sig_downregulated_hclust_heatmap.pdf"),
    useDingbats = FALSE,
    family="Arial",
    fonts="Arial")

draw(Rev_1_F5_transcript_hclust_down_heat)

dev.off()

### Fit the log2-norm counts ------------------------------------------------

#### Up ------------------------------------------------

# Prepare for plotting
Rev_1_F5_transcript_test_up_klus_cond_forPlot <- Rev_1_F5_transcript_test_up_klus_cond %>% 
  pivot_longer(cols=-c(transcript_id, DTE_cluster_up),
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

Rev_1_F5_transcript_test_up_klus_cond_n <- Rev_1_F5_transcript_test_up_klus_cond_forPlot %>% 
  group_by(DTE_cluster_up) %>% 
  distinct(transcript_id,.keep_all = TRUE) %>% 
  summarize(n=n()) %>% 
  dplyr::rename("DTE_cluster" = "DTE_cluster_up")

#### Down ------------------------------------------------

# Prepare for plotting
Rev_1_F5_transcript_test_down_klus_cond_forPlot <- Rev_1_F5_transcript_test_down_klus_cond %>% 
  pivot_longer(cols=-c(transcript_id, DTE_cluster_down),
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

Rev_1_F5_transcript_test_down_klus_cond_n <- Rev_1_F5_transcript_test_down_klus_cond_forPlot %>% 
  group_by(DTE_cluster_down) %>% 
  distinct(transcript_id,.keep_all = TRUE) %>% 
  summarize(n=n()) %>% 
  dplyr::rename("DTE_cluster" = "DTE_cluster_down")

 #### Combined Up/Down ------------------------------------------------
Rev_1_F5_transcript_test_combined_klus <- Rev_1_F5_transcript_test_down_klus_cond_forPlot %>% 
  # filter(exp_cond %in% c("UPF1_AID",
  #                        "UPF1_AID_Rec")) %>% 
  filter(!condition %in% c("UPF1_Nter_48h",
                           "UPF1_Nter_24h_R0h")) %>% 
  dplyr::rename("DTE_cluster" = "DTE_cluster_down") %>% 
  bind_rows(Rev_1_F5_transcript_test_up_klus_cond_forPlot %>% 
              # filter(exp_cond %in% c("UPF1_AID",
              #                        "UPF1_AID_Rec")) %>% 
              filter(!condition %in% c("UPF1_Nter_48h",
                                       "UPF1_Nter_24h_R0h")) %>% 
              dplyr::rename("DTE_cluster" = "DTE_cluster_up")) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster,
                                    "up 1:early",
                                    "down 1:early",
                                    "up 2:delayed",
                                    "down 2:delayed",
                                    "up 3:late",
                                    "down 3:late",
                                    "up 4:inverse",
                                    "down 4:inverse"))) %>% 
  # group_by(timepoint,DTE_cluster_down) %>% 
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
    # data=Rev_1_F5_transcript_test_down_klus_cond_forPlot %>% 
    #                              filter(exp_cond %in% c("UPF1_AID",
    #                                                     "UPF1_AID_Rec")) %>%
    #                              filter(!condition %in% c("UPF1_Nter_48h")) %>% 
    #                              dplyr::rename("DTE_cluster" = "DTE_cluster_down") %>% 
    #                              bind_rows(Rev_1_F5_transcript_test_up_klus_cond_forPlot %>% 
    #                                          filter(exp_cond %in% c("UPF1_AID",
    #                                                                 "UPF1_AID_Rec")) %>%
    #                                          filter(condition != "UPF1_Nter_48h") %>% 
    #                                          dplyr::rename("DTE_cluster" = "DTE_cluster_up")),
    aes(group=transcript_id,
        color=DTE_cluster),
    alpha=0.5,
    linewidth=0.1),
    dpi=600, dev = "cairo") +
  stat_smooth(
    # data=Rev_1_F5_transcript_test_down_klus_cond_forPlot %>% 
    #             filter(exp_cond %in% c("UPF1_AID",
    #                                    "UPF1_AID_Rec")) %>%
    #             filter(!condition %in% c("UPF1_Nter_48h")) %>% 
    #             dplyr::rename("DTE_cluster" = "DTE_cluster_down") %>% 
    #             bind_rows(Rev_1_F5_transcript_test_up_klus_cond_forPlot %>% 
    #                         filter(exp_cond %in% c("UPF1_AID",
    #                                                "UPF1_AID_Rec")) %>%
    #                         filter(condition != "UPF1_Nter_48h") %>% 
    #                         dplyr::rename("DTE_cluster" = "DTE_cluster_up")),
    color="black",method = "gam", formula = y ~ s(x,k=7, bs="cs"), size = 1) +
  geom_boxplot(aes(group=interaction(timepoint,exp_cond),
                   y=score,
                   fill=exp_cond),
               outlier.shape = NA,
               fatten=2,
               linewidth=0.25) +
  geom_text(data=bind_rows(Rev_1_F5_transcript_test_up_klus_cond_n,
                           Rev_1_F5_transcript_test_down_klus_cond_n),
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
  facet_wrap(~DTE_cluster,ncol=2) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(total_width = unit(65, "mm"),
                   total_height = unit(100, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_C_transcript_test_combined_klus.pdf"),
       Rev_1_F5_transcript_test_combined_klus,
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (D) NMDRHT ImpulseDE2 -----------------------------------------------------------
###

load("Resources/ImpulseDE2/Rev_1_F5_ImpulseDE2.rds")

### Overarching analyses ----------------------------------------------------

# How many significant genes per biotype
Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  mutate(type = case_when(!transcript_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          transcript_biotype == "protein_coding" ~ "coding",
                          transcript_biotype == "nonsense_mediated_decay" ~ "NMD",
                          transcript_biotype == "lncRNA" ~ "lncRNA")) %>% 
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
  labs(x="Number of significant transcripts",
       y="NMDRHT\ntranscript biotype",
       fill="Best fit") +
  guides(fill=guide_legend(reverse=T))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_D_ImpulseDE2_biotypes_sig.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Frequency of onset/offset time

Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  mutate(timecourse=fct_rev(fct_relevel(timecourse,
                                        "combined"))) %>% 
  filter(bestFit != "noFit") %>% 
  dplyr::count(timecourse,bestFit)

# How many sig. genes per DTE_cluster
Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  # filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  mutate(timecourse=(fct_relevel(timecourse,
                                 "combined"))) %>% 
  mutate(bestFit = case_when(sigImpulseDE2 == "n.s." ~ "n.s.",
                             TRUE ~ bestFit)) %>% 
  filter(!DTE_cluster %in% c("not_expressed", "complex", "expressed")) %>% 
  dplyr::count(DTE_cluster, timecourse, bestFit) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  group_by(DTE_cluster, timecourse) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  mutate(bestFit = fct_relevel(bestFit,
                               "Sigmoid",
                               "Impulse",
                               "noFit",
                               "n.s.")) %>%
  ggplot(aes(y=DTE_cluster,
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
  labs(x="Fraction of significant transcripts",
       y="Sig. regulated DTE cluster",
       fill="Best fit") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, reverse=T))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_D_ImpulseDE2_DTE_cluster_percent.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Boxplot padj per DTE cluster
Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  # filter(sigImpulseDE2 == "sig.") %>%
  mutate(timecourse="combined") %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(!DTE_cluster %in% c("not_expressed", "complex", "expressed")) %>% 
  filter(bestFit %in% c("Impulse", "Sigmoid")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
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
             y=DTE_cluster,
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
  geom_boxplot(outlier.shape = NA,
               color="black",
               linewidth=0.2) +
  scale_fill_manual(values=c("Impulse" = "#4D5F8E",
                             "Sigmoid" = "#C582B2")) +
  facet_wrap(~timecourse) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_D_ImpulseDE2_DTE_cluster_boxplot_padj.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  


# Boxplot per DTE cluster

Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late")) %>% 
  pivot_longer(cols=c(t1,t2),
               names_to = "Time",
               values_to = "Values") %>% 
  mutate(Time = fct_rev(fct_relevel(Time,
                                    "t1",
                                    "t2"))) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  filter(bestFit == "Impulse") %>% 
  ggplot(aes(x=Values,
             y=DTE_cluster,
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
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_D_ImpulseDE2_DTE_cluster_boxplot_times.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

Rev_1_F5_ImpulseDE2_Combined_df_Cluster_Param %>% 
  filter(sigImpulseDE2 == "sig.") %>%
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late")) %>% 
  pivot_longer(cols=c(t1,t2),
               names_to = "Time",
               values_to = "Values") %>% 
  filter(bestFit == "Impulse") %>% 
  group_by(DTE_cluster, Time) %>% 
  summarize(median = round(median(Values),1))

###
## Rev_1 - Figure 5 - (E) NMD relevance -----------------------------------------------------------
###

load("Resources/NMD_relevance/edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance.rds")

### Plot NMD relevance - barplot --------------------------------------------------------------------
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(DTE_cluster, UpDown, NMD_n_sig_tx_perc) %>% 
  mutate(DTE_cluster = fct_rev(DTE_cluster)) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  # mutate(n = case_when(UpDown == "up" ~ n,
  #                      UpDown == "down" ~ -n)) %>% 
  ggplot(aes(y=DTE_cluster,
             x=n,
             fill=NMD_n_sig_tx_perc)) +
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

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_DTE_cluster_NMD_relevance_perc_fill.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### Plot NMD relevance bins - barplot --------------------------------------------------------------------
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late")) %>% 
  mutate(NMD_bin_tx = fct_rev(NMD_bin_tx)) %>% 
  dplyr::count(UpDown,NMD_bin_tx) %>% 
  mutate(n=case_when(UpDown == "up" ~ n,
                     UpDown == "down" ~ -n)) %>% 
  mutate(NMD_bin_tx = fct_rev(NMD_bin_tx)) %>% 
  ggplot(aes(x=n,
             y=NMD_bin_tx,
             fill=NMD_bin_tx)) +
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
  labs(fill="NMD\nconfidence\nbins") +
  scale_x_continuous(expand = expansion(c(0.25, 0.25))) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_DTE_cluster_NMD_confidence_bin_numbers.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

# Number per DT_cluster
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin %>% 
  dplyr::count(NMD_bin_tx, UpDown) %>% 
  group_by(UpDown) %>% 
  mutate(n_per = round(100*(n/sum(n)),1))

### Plot NMD relevance / BIN - barplot --------------------------------------------------------------------
edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "down 1:early",
                            "up 2:delayed",
                            "down 2:delayed",
                            "up 3:late",
                            "down 3:late"
  )) %>% 
  group_by(DTE_cluster, UpDown, NMD_bin_tx) %>% 
  mutate(DTE_cluster = fct_rev(DTE_cluster)) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(NMD_bin_tx = fct_rev(fct_relevel(NMD_bin_tx,
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
  ggplot(aes(y=DTE_cluster,
             x=n,
             fill=NMD_bin_tx)) +
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

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_DTE_cluster_NMD_confidence_perc_fill_bin.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### UPF1 log2FC bin ---------------------------------------------------------

edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late")) %>% 
  dplyr::select(transcript_id,NMD_n_sig_tx_perc,NMD_bin_tx, UpDown) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(experimentSet %in% c("HeLa_UPF1_KD_cytoplasmic_Longman_2020",
                                          # "HEK293_SMG567_KO_KD_Boehm_2021",
                                          # "HEK293_SMG67_KD_Britto_Borges_2024",
                                          # "HeLa_SMG67_KD_Britto_Borges_2024",
                                          # "MCF7_SMG67_KD_Britto_Borges_2024",
                                          # "U2OS_SMG67_KD_Britto_Borges_2024",
                                          # "HEK293_UPF3_dKO_Wallmeroth_2022",
                                          # "HEK293_CASC3_KO_Gerbracht_2020",
                                          "HCT116_UPF1_AID_degradation_this_Study",
                                          # "HCT116_UPF1_AID_recovery_this_Study",
                                          "HCT116_UPF1_FKBP_degradation_this_Study",
                                          "HEK293_UPF1_FKBP_degradation_this_Study"
                                          # "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                                          # "HFF_SMG1i_this_Study",
                                          # "HUVEC_SMG1i_this_Study"
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
                                                  "HEK293_SMG567_KO_KD_Boehm_2021",
                                                  "HEK293_SMG67_KD_Britto_Borges_2024",
                                                  "HeLa_SMG67_KD_Britto_Borges_2024",
                                                  "MCF7_SMG67_KD_Britto_Borges_2024",
                                                  "U2OS_SMG67_KD_Britto_Borges_2024",
                                                  "HCT116_SMG89KO_SMG1i_Kueckelmann_2024",
                                                  "HFF_SMG1i_this_Study",
                                                  "HUVEC_SMG1i_this_Study"))) %>% 
  arrange(experiment_set) %>% 
  filter(!is.na(condition_2)) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  # group_by(condition_2, bin) %>% 
  # summarise_at(vars(log2FoldChange),
  #              list(Q1=~quantile(., probs = 0.25, na.rm = TRUE),
  #                   median=median,
  #                   Q3=~quantile(., probs = 0.75, na.rm = TRUE)), na.rm = TRUE) %>% 
  ggplot(aes(x=logFC,
             y=condition_2,
             fill=NMD_bin_tx)) +
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
  geom_boxplot(outliers = FALSE,
               # coef = 0,
               color="black",
               linewidth=0.1) +
  scale_fill_manual(values = c("[0,25]" = "#CACFD0",
                               "(25,50]" =  "#647588",
                               "(50,75]" =  "#B09771",
                               "(75,100]" =  "#E5BF86")) +
  scale_x_continuous(name = "log2FC",
                     breaks = c(0,2,4),
                     labels = c(0,2,4),
                     limits = c(-2,10)) +
  force_panelsizes(rows = unit(6*8+0.4, "mm"),
                   cols = unit(20+0.4, "mm")) 


ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_DGE_cluster_NMD_confidence_bin_log2FC.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### UPF1 log2FC bin heatmap ---------------------------------------------------------

edgeR_DTE_NMDRHT_combined_DTE_cluster_NMD_relevance_complete_bin %>% 
  filter(DTE_cluster %in% c("up 1:early",
                            "up 2:delayed",
                            "up 3:late")) %>% 
  dplyr::select(transcript_id, NMD_n_sig_tx_perc,NMD_bin_tx, UpDown) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
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
                                                  "HUVEC_SMG1i_this_Study"))) %>% 
  arrange(experiment_set) %>% 
  filter(!is.na(condition_2)) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>% 
  mutate(condition_2 = (fct_inorder(as_factor(condition_2)))) %>% 
  # distinct(condition_2)
  group_by(condition_2, NMD_bin_tx) %>%
  summarise_at(vars(logFC),
               list(Q1=~quantile(., probs = 0.25, na.rm = TRUE),
                    median=median,
                    Q3=~quantile(., probs = 0.75, na.rm = TRUE)), na.rm = TRUE) %>%
  ggplot(aes(y=NMD_bin_tx,
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
       fill="median\nDTE\nlog2FC") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                direction = "vertical",
                                label.hjust = 0.5,
                                label.vjust = 0.5))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_median_log2FC_perBin.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

### DTE cluster per NMDRHT biotype --------------------------------------------------------------
NMDRHT.v1.2_MainTable %>%  
  filter(!DTE_cluster %in% c("not_expressed")) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster,
                                    "up 1:early",
                                    "up 2:delayed",
                                    "up 3:late",
                                    "up 4:inverse",  
                                    "expressed",
                                    "not_expressed",
                                    "down 4:inverse",
                                    "down 3:late",
                                    "down 2:delayed",
                                    "down 1:early"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  dplyr::count(NMD_tx_status, DTE_cluster) %>% 
  group_by(NMD_tx_status) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  print(n=60)

NMDRHT.v1.2_MainTable %>%  
  filter(!DTE_cluster %in% c("not_expressed")) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",  
                                           "expressed",
                                           "not_expressed",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  dplyr::count(NMD_tx_status, DTE_cluster) %>% 
  group_by(NMD_tx_status) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  ggplot(aes(x=NMD_tx_status,
             y=n,
             fill=DTE_cluster)) +
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
           position="fill",
           linewidth=0.2) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                             "up 2:delayed" = "#b2182b",
                             "up 3:late" = "#d6604d",
                             "up 4:inverse" = "#cfbeb4",
                             "down 4:inverse" = "#b9c3c8",
                             "down 3:late" = "#92c5de",
                             "down 2:delayed" = "#4393c3",
                             "down 1:early" = "#053061",
                             "expressed" = "lightgray",
                             "not_expressed" = "gray5")) +
  guides(fill=guide_legend(reverse=F)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(7*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_NMD_status_DTE_cluster_fraction_woNonExpressed.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### 12h IAA log2FC "expressed "NMD" --------------------------------------------------------------


#### log2FC ------------------------------------------------------------------

NMDRHT.v1.2_MainTable %>%  
  filter(DTE_cluster %in% c("expressed")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 == "UPF1_Nter_12h")) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  mutate(significant = case_when(FDR < 0.0001 & abs(logFC) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  ggplot(aes(x=NMD_tx_status,
             y=logFC,
             fill=significant)) +
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
  geom_boxplot(outliers = FALSE,
               # coef = 0,
               color="black",
               linewidth=0.1) +
  geom_hline(yintercept=1,
             color="#963B5A",
             linewidth=0.25,
             linetype="solid") +
  scale_fill_manual(values=c("TRUE" = "#963B5A",
                             "FALSE" = "lightgray")) +
  guides(fill=guide_legend(reverse=F)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(7*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_NMD_status_Expressed_log2FC_12h.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### FDR ------------------------------------------------------------------
NMDRHT.v1.2_MainTable %>%  
  filter(DTE_cluster %in% c("expressed")) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 == "UPF1_Nter_12h")) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  mutate(significant = case_when(FDR < 0.0001 & abs(logFC) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  ggplot(aes(x=NMD_tx_status,
             y=-log10(FDR),
             fill=significant)) +
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
  geom_boxplot(outliers = FALSE,
               # coef = 0,
               color="black",
               linewidth=0.1) +
  geom_hline(yintercept=-log10(0.0001),
             color="#963B5A",
             linewidth=0.25,
             linetype="solid") +
  scale_fill_manual(values=c("TRUE" = "#67001f",
                             "FALSE" = "lightgray")) +
  guides(fill=guide_legend(reverse=F)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(7*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_E_NMD_status_Expressed_FDR_12h.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (F) PARE analysis -----------------------------------------------------------
###

PARE_edgeR_combined <- read_csv("Resources/PARE/PARE_edgeR_combined.csv")

#### NMD status --------------------------------------------------------------

PARE_edgeR_combined %>% 
  # filter(DTE_cluster %in% c("up 1:early", "up 2:delayed", "up 3:late", "expressed")) %>% 
  pivot_longer(cols=c(L2FC_XS6, L2FC_XU1),
               names_to="condition_2",
               values_to="logFC") %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
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
  scale_fill_manual(values=c("L2FC_XS6" = "#624B27",
                             "L2FC_XU1" = "#006666")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(0,2,4),
                     # labels = c(0,2,4),
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(5*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_F_log2FC_PARE_seq_NMD_status_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### NMD reason --------------------------------------------------------------

PARE_edgeR_combined %>% 
  # filter(DTE_cluster %in% c("up 1:early", "up 2:delayed", "up 3:late", "expressed")) %>% 
  pivot_longer(cols=c(L2FC_XS6, L2FC_XU1),
               names_to="condition_2",
               values_to="logFC") %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=NMD_tx_reason_simple,
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
  scale_fill_manual(values=c("L2FC_XS6" = "#624B27",
                             "L2FC_XU1" = "#006666")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(0,2,4),
                     # labels = c(0,2,4),
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(7*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_F_log2FC_PARE_seq_NMD_reason_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### NMD relevance bin --------------------------------------------------------------

PARE_edgeR_combined %>% 
  # filter(DTE_cluster %in% c("up 1:early", "up 2:delayed", "up 3:late", "expressed")) %>% 
  pivot_longer(cols=c(L2FC_XS6, L2FC_XU1),
               names_to="condition_2",
               values_to="logFC") %>% 
  # filter(UpDown == "up") %>% 
  mutate(NMD_bin_updown = paste0(UpDown,"_",NMD_bin_tx)) %>% 
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
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=NMD_bin_updown,
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
  scale_fill_manual(values=c("L2FC_XS6" = "#624B27",
                             "L2FC_XU1" = "#006666")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(0,2,4),
                     # labels = c(0,2,4),
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(9*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_F_log2FC_PARE_seq_NMD_relevance_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### DTE_cluster --------------------------------------------------------------

PARE_edgeR_combined %>% 
  # filter(DTE_cluster %in% c("up 1:early", "up 2:delayed", "up 3:late", "expressed")) %>% 
  pivot_longer(cols=c(L2FC_XS6, L2FC_XU1),
               names_to="condition_2",
               values_to="logFC") %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=DTE_cluster,
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
  scale_fill_manual(values=c("L2FC_XS6" = "#624B27",
                             "L2FC_XU1" = "#006666")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(0,2,4),
                     # labels = c(0,2,4),
  ) +
  coord_cartesian(ylim = c(-2,2)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~DTE_cluster,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(7*4+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_F_log2FC_PARE_seq_DTE_cluster_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (G) Overarching analyses -----------------------------------------------------------
###


### log2FC 12h UPF1 depletion -------------------------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=NMD_tx_reason_simple,
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
             linewidth=0.25,
             linetype="dashed") +
  geom_hline(yintercept=1,
             linewidth=0.25,
             color="#963B5A",
             linetype="dashed") +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_tx_status,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(20+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_UPF1_12h_NMD_reason_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Structural Category  ----------------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  filter(structural_category_simple %in% c("FSM", "NIC", "NNIC")) %>% 
  ggplot(aes(x=structural_category_simple,
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
             linewidth=0.25,
             linetype="dashed") +
  geom_hline(yintercept=1,
             linewidth=0.25,
             color="#963B5A",
             linetype="dashed") +
  geom_boxplot(outliers = FALSE,
               color="black",
               alpha=0.5,
               linewidth=0.1) +
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_tx_status,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(10+0.4, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_UPF1_12h_structuralCategory_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### NMD bin log2FC  ----------------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  filter(UpDown == "up") %>% 
  mutate(NMD_bin_updown = paste0(UpDown,"_",NMD_bin_tx)) %>% 
  mutate(NMD_bin_updown = fct_rev(fct_relevel(NMD_bin_updown,
                                       "up_(75,100]",
                                       "up_(50,75]",
                                       "up_(25,50]",
                                       "up_[0,25]",
                                       "down_[0,25]",
                                       "down_(25,50]",
                                       "down_(50,75]",
                                       "down_(75,100]",
  ))) %>% 
  ggplot(aes(x=NMD_bin_updown,
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
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_tx_status,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_UPF1_12h_NMD_bin_Up_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### DTE cluster log2FC  ----------------------------------------------------

NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  ggplot(aes(x=DTE_cluster,
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
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_tx_status,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(21, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_UPF1_12h_DTE_cluster_Up_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### bakR gene & DTE tx -------------------------------------------------------------------

GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

NMDRHT.v1.2_MainTable %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id, kdeg_conclusion, RNA_conclusion, Mech_conclusion)) %>% 
  mutate(overall_conclusion = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ "Stabilized",
                                        kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ "Destabilized",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ "Increased Syn.",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ "Decreased Syn.",
                                        TRUE ~ "n.s.")) %>% 
  mutate(overall_conclusion = fct_rev(fct_relevel(overall_conclusion,
                                          "Stabilized",
                                          "Increased Syn.",
                                          "n.s.",
                                          "Decreased Syn.",
                                          "Destabilized"))) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h",
                                        "UPF1_FKBP_HCT_12h",
                                        "UPF1_FKBP_HEK_12h"))) %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  ggplot(aes(x=overall_conclusion,
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
  scale_fill_manual(values=c("UPF1_Nter_12h" = "#00696D",
                             "UPF1_FKBP_HCT_12h" = "#6B56A7",
                             "UPF1_FKBP_HEK_12h" = "#963B5A")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_tx_status,
             nrow=1,
             scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_UPF1_12h_DTE_gene_Conclusion_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### tx -vs- gene-level ----------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                      DGE_cluster)) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
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
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "expressed",
                                           "not_expressed",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  dplyr::count(DGE_cluster, DTE_cluster, NMD_50nt_rule) %>% 
  ggplot(aes(x=DTE_cluster,
             y=n,
             fill=DGE_cluster)) +
  geom_col(position="fill") +
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
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~NMD_50nt_rule)


##### Correlation Scatter -----------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(GENCODE_v42_MainTable %>% dplyr::select(gene_id,
                                                      DGE_cluster)) %>% 
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, log2FoldChange)) %>% 
  ggplot(aes(x=logFC,
             y=log2FoldChange)) +
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
  geom_hline(yintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  geom_vline(xintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                 alpha=0.75,
                 color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "eq"), sep = "\n"),
               size = 5*0.30,
               color="gray20",
               label.y = "bottom", label.x = "right") +
  facet_wrap(~NMD_50nt_rule) +
  scale_color_viridis_c(option="mako") +
  coord_fixed(ratio=1,
              clip = "off",
              xlim = c(-10,10),
              ylim = c(-10,10)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_log2FC_DGE_DTE_scatter_12h_AID_UPF1.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### ISAR NMDRHT - all ----------------------------------------------------------
ISAR_DTU_NMDRHT_combined %>% 
  filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study")) %>% 
  dplyr::rename("transcript_id" = "isoform_id") %>% 
  # filter(isoform_switch_q_value < 0.0001) %>% 
  dplyr::select(transcript_id, IF1, IF2, dIF, condition_2) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  # mutate(condition_2 = fct_rev(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "control_48h",
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
                                           "UPF1_FKBP_HEK_12h"))) %>%
  # filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed",
  #                            "down 3:late",
  #                            "down 2:delayed",
  #                            "down 1:early")) %>%
  # filter(NMD_tx_reason == "none") %>%
  # filter(utr3_len != 0) %>%
  mutate(DTE_cluster = (fct_relevel(DTE_cluster,
                                    "up 1:early",
                                    "up 2:delayed",
                                    "up 3:late",
                                    "expressed"))) %>%
  group_by(condition_2, NMD_tx_reason_simple) %>% 
  summarize(median_dIF = median(dIF)) %>% 
  ungroup() %>% 
  ggplot(aes(y=condition_2,
             x=NMD_tx_reason_simple,
             fill=median_dIF)) +
  theme(legend.position="right", 
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       # limits = c(-1,4),
                       na.value = "grey90") +
  # facet_wrap(~DTE_cluster, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))


ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_ISAR_all_NMD_reason_heatmap.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### ISAR NMDRHT - scatter ----------------------------------------------------------
ISAR_DTU_NMDRHT_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  dplyr::rename("transcript_id" = "isoform_id") %>% 
  mutate(isoform_switch_q_value = replace(isoform_switch_q_value, isoform_switch_q_value == 0, 1e-300)) %>% 
  # filter(isoform_switch_q_value < 0.0001) %>% 
  dplyr::select(transcript_id, IF1, IF2, dIF, condition_2, isoform_switch_q_value) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  mutate(DTE_cluster = (fct_relevel(DTE_cluster,
                                    "up 1:early",
                                    "up 2:delayed",
                                    "up 3:late",
                                    "expressed"))) %>%
  ggplot(aes(x=dIF,
             y=-log10(isoform_switch_q_value))) +
  theme(legend.position="right", 
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
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  # scale_color_manual(values=c("none" = "#214D65",
  #                            "lncRNA" = "#76716E",
  #                            "AS_NMD"="#E5BF85",
  #                            "AS_NMD_UTR3" = "#C1AC65",
  #                            "novel" = "#9A9A49",
  #                            "uORF" = "#708831",
  #                            "overl_uORF" = "#43761E",
  #                            "other" = "#D9DADA")) +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  scale_color_viridis_c(option="mako") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))


ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_ISAR_volcano_50nt_rule.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### ISAR NMDRHT - dIF vs log2FC ----------------------------------------------------------
ISAR_DTU_NMDRHT_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  dplyr::rename("transcript_id" = "isoform_id") %>% 
  # filter(isoform_switch_q_value < 0.0001) %>% 
  dplyr::select(transcript_id, IF1, IF2, dIF, condition_2, isoform_switch_q_value) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  ggplot(aes(x=logFC,
             y=dIF)) +
  geom_hline(yintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  geom_vline(xintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  theme(legend.position="right", 
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
  ggpmisc::stat_poly_line(show.legend=FALSE,
                 alpha=0.75,
                 color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "eq"), sep = "\n"),
               size = 5*0.30,
               color="gray20",
               label.y = "bottom", label.x = "right") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  scale_color_viridis_c(option="mako") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_ISAR_DTU_vs_DTE.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### ISAR NMDRHT - DGE vs dIF ----------------------------------------------------------
ISAR_DTU_NMDRHT_combined %>% 
  filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
  dplyr::rename("transcript_id" = "isoform_id") %>% 
  # filter(isoform_switch_q_value < 0.0001) %>% 
  dplyr::select(transcript_id, IF1, IF2, dIF, condition_2, isoform_switch_q_value) %>% 
  left_join(NMDRHT.v1.2_MainTable) %>% 
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, log2FoldChange)) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  ggplot(aes(y=log2FoldChange,
             x=dIF)) +
  geom_hline(yintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  geom_vline(xintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  theme(legend.position="right", 
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
  ggpmisc::stat_poly_line(show.legend=FALSE,
                 alpha=0.75,
                 color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "eq"), sep = "\n"),
               size = 5*0.30,
               color="gray20",
               label.y = "bottom", label.x = "right") +
  facet_wrap(~NMD_50nt_rule, nrow=1) +
  scale_color_viridis_c(option="mako") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_ISAR_DGE_vs_DTU.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


##### NMD_tx_status per DTE_cluster ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("expressed", "not_expressed")) %>% 
  mutate(NMD_tx_status = (fct_relevel(NMD_tx_status,
                                      "coding",
                                      "predicted_coding",
                                      "NMD",
                                      "predicted_NMD",
                                      "mixed",
                                      "lncRNA"))) %>%
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "expressed",
                                           "not_expressed",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  dplyr::count(DTE_cluster, NMD_tx_status) %>% 
  group_by(DTE_cluster) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  ggplot(aes(x=DTE_cluster,
             y=n_per,
             fill=NMD_tx_status)) +
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
  geom_col() +
  geom_text(aes(label = case_when(n_per > 25 & NMD_tx_status == "coding" ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 25 & NMD_tx_status != "coding" ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD"="#E5BF86",
                             "mixed" = "#624B27",
                             "predicted_NMD" = "#B09771",
                             "predicted_coding" = "#287DAB")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(8*3, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_NMD_tx_status_DTE_cluster.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### NMD_reason per DTE_cluster ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("expressed", "not_expressed")) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "expressed",
                                           "not_expressed",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  dplyr::count(DTE_cluster, NMD_tx_reason_simple) %>% 
  group_by(DTE_cluster) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  ggplot(aes(x=DTE_cluster,
             y=n_per,
             fill=NMD_tx_reason_simple)) +
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
  geom_col() +
  geom_text(aes(label = case_when(n_per > 25 & NMD_tx_reason_simple == "none" ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 25 & NMD_tx_reason_simple != "none" ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "other" = "#D9DADA")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(8*3, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_NMD_reason_DTE_cluster.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### StructuralCategory per DTE_cluster ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(!DTE_cluster %in% c("expressed", "not_expressed")) %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "up 4:inverse",
                                           "expressed",
                                           "not_expressed",
                                           "down 4:inverse",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  dplyr::count(DTE_cluster, structural_category_simple_2) %>% 
  group_by(DTE_cluster) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  ggplot(aes(x=DTE_cluster,
             y=n_per,
             fill=structural_category_simple_2)) +
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
  geom_col() +
  geom_text(aes(label = case_when(n_per > 25 & structural_category_simple_2 == "FSM" ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 25 & structural_category_simple_2 != "FSM"  ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("FSM" = "#414487",
                             "NIC" = "#7AD151",
                             "NNIC" = "#BBDF27",
                             "other" = "gray80")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(8*3, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_StructuralCategory_DTE_cluster.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### DTE NMD bin NMD reason ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  dplyr::count(NMD_bin_tx, UpDown, NMD_tx_status) %>% 
  arrange(desc(n)) %>% 
  group_by(UpDown, NMD_bin_tx) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>% 
  mutate(NMD_bin = paste0(UpDown,"_",NMD_bin_tx)) %>% 
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
  mutate(NMD_tx_status = (fct_relevel(NMD_tx_status,
                                      "coding",
                                      "predicted_coding",
                                      "NMD",
                                      "predicted_NMD",
                                      "mixed",
                                      "lncRNA"))) %>%
  ggplot(aes(x=NMD_bin,
             y=n_per,
             fill=NMD_tx_status)) +
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
  geom_col(position="stack") +
  geom_text(aes(label = case_when(n_per > 40 & NMD_tx_status == "coding" ~ paste0(n,
                                                                                  " (",
                                                                                  round(n_per,0),
                                                                                  "%)"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            angle=90,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 40 & NMD_tx_status != "coding" ~ paste0(n,
                                                                                  " (",
                                                                                  round(n_per,0),
                                                                                  "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            angle=90,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD"="#E5BF86",
                             "mixed" = "#624B27",
                             "predicted_NMD" = "#B09771",
                             "predicted_coding" = "#287DAB")) +
  facet_wrap(~UpDown,
             scales="free_x") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_DTE_NMD_bin_NMD_tx_status.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### High NMD relevant - NMD reason ------------------------------------------------------------
NMDRHT.v1.2_MainTable %>% 
  filter(NMD_bin_tx == "(75,100]" & UpDown == "up") %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  dplyr::count(NMD_tx_reason_simple) %>% 
  arrange(desc(n)) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ggplot(aes(x="all",
             y=n_per,
             fill=NMD_tx_reason_simple)) +
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
  geom_col(position="stack") +
  # geom_text(aes(label = case_when(n_per > 40 & NMD_tx_status == "coding" ~ paste0(n,
  #                                                                                 " (",
  #                                                                                 round(n_per,0),
  #                                                                                 "%)"),
  #                                 TRUE ~ "")),
  #           color="white",
  #           size = 4*0.36,
  #           angle=90,
  #           position = position_stack(vjust = 0.5)
  # ) +
  geom_text(aes(label = case_when(n_per > 40 ~ paste0(n,
                                                      " (",
                                                      round(n_per,0),
                                                      "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            angle=90,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "other" = "#D9DADA")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(3, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_G_DTE_High_up_NMD_bin_NMD_reason.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 5 - (H) EZbakR tx-level -----------------------------------------------------------
###

# Load EZbakR analysis of transcript-level NMDRHT-based
UPF1_NMDRHT_EZbakR_TEC_combined <- read_csv(file=file.path("Resources/EZbakR_NMDRHT/UPF1_NMDRHT_EZbakR_TEC_combined_TPM02.csv"))

##### gene-vs-tx-level bakR -----------------------------------------------------
UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  # filter(condition=="UPF1_12h") %>% 
  ggplot(aes(x=L2FC_kdeg_tx,
             y=L2FC_kdeg)) +
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
  geom_hline(yintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  geom_vline(xintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 4,size = 0.25),
                     dpi=1200, dev = "cairo") +
  ggpmisc::stat_poly_line(show.legend=FALSE,
                          alpha=0.75,
                          color="gray") +
  ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label(c("adj.R2", "eq"), sep = "\n"),
                        size = 5*0.30,
                        color="gray20",
                        label.y = "bottom", label.x = "right") +
  facet_wrap(~condition, nrow=1) +
  scale_color_viridis_c(option="mako") +
  coord_fixed(ratio=1,
              clip = "off",
              xlim = c(-4,4),
              ylim = c(-4,4)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_log2FC_kdeg_gene_tx_scatter_12h_AID_UPF1.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR DTE cluster  ----------------------------------------------------

UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  mutate(NMD_tx_status = fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA")) %>% 
  filter(!DTE_cluster %in% c("up 4:inverse", "down 4:inverse", "not_expressed", "down 3:late",
                             "down 2:delayed",
                             "down 1:early")) %>% 
  mutate(DTE_cluster = fct_rev(fct_relevel(DTE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed"))) %>% 
  ggplot(aes(x=DTE_cluster,
             y=L2FC_kdeg_tx,
             fill=DTE_cluster)) +
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
               alpha=0.75,
               linewidth=0.1) +
  scale_fill_manual(values=c("up 1:early" = "#67001f",
                              "up 2:delayed" = "#b2182b",
                              "up 3:late" = "#d6604d",
                              "up 4:inverse" = "#cfbeb4",
                              "down 4:inverse" = "#b9c3c8",
                              "down 3:late" = "#92c5de",
                              "down 2:delayed" = "#4393c3",
                              "down 1:early" = "#053061")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(-3,-2,-1,0,1),
                     breaks = c(-3,-2,-1,0,1),
                     limits = c(-3.5,1.5)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~NMD_tx_status,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_DTE_cluster_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR NMD relevance bin ----------------------------------------------------

UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  filter(UpDown != "down") %>%
  mutate(NMD_bin_updown = paste0(UpDown,"_",NMD_bin_tx)) %>% 
  mutate(NMD_bin_updown = fct_rev(fct_relevel(NMD_bin_updown,
                                              "up_(75,100]",
                                              "up_(50,75]",
                                              "up_(25,50]",
                                              "up_[0,25]",
                                              "n.s._n.s."
  ))) %>% 
  ggplot(aes(x=NMD_bin_updown,
             y=L2FC_kdeg_tx,
             fill=NMD_bin_updown)) +
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
               alpha=0.75,
               linewidth=0.1) +
  scale_fill_manual(values = c("up_[0,25]" = "#CACFD0",
                               "up_(25,50]" =  "#647588",
                               "up_(50,75]" =  "#B09771",
                               "up_(75,100]" =  "#E5BF86",
                               "n.s._n.s." = "white",
                               "down_[0,25]" = "#CACFD0",
                               "down_(25,50]" =  "#647588",
                               "down_(50,75]" =  "#B09771",
                               "down_(75,100]" =  "#E5BF86")) +
  scale_y_continuous(name = "log2FC",
                     # breaks = c(-3,-2,-1,0,1),
                     breaks = c(-3,-2,-1,0,1),
                     limits = c(-3.5,1.5)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~NMD_tx_status,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(15, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_NMD_bin_boxplot.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR NMD status ----------------------------------------------------

UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  # filter(UpDown == "up") %>% 
  mutate(NMD_tx_status = (fct_relevel(NMD_tx_status,
                                     "coding",
                                     "predicted_coding",
                                     "NMD",
                                     "predicted_NMD",
                                     "mixed",
                                     "lncRNA"))) %>%
  ggplot(aes(x=NMD_tx_status,
             y=L2FC_kdeg_tx,
             fill=NMD_tx_status)) +
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
               alpha=0.75,
               linewidth=0.1) +
  scale_fill_manual(values=c("coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD"="#E5BF86",
                             "mixed" = "#624B27",
                             "predicted_NMD" = "#B09771",
                             "predicted_coding" = "#287DAB")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~NMD_tx_status,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(18, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_NMD_status.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR NMD status ----------------------------------------------------

UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  # filter(UpDown == "up") %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  ggplot(aes(x=NMD_tx_reason_simple,
             y=L2FC_kdeg_tx,
             fill=NMD_tx_reason_simple)) +
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
               alpha=0.75,
               linewidth=0.1) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "other" = "#D9DADA")) +
  # scale_y_continuous(name = "log2FC",
  #                    # breaks = c(0,2,4),
  #                    # labels = c(0,2,4),
  #                    limits = c(-8,8)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # facet_wrap(~NMD_tx_status,
  #            nrow=1,
  #            scales="free_x") +
  force_panelsizes(rows = unit(20+0.4, "mm"),
                   cols = unit(18, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_NMD_reason.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR kdeg difference tx-vs-gene ------------------------------------------------------------
UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  filter(kdeg_tx_conclusion == "Stabilized") %>% 
  mutate(diff_kdeg = L2FC_kdeg_tx-L2FC_kdeg) %>% 
  arrange(diff_kdeg) %>% 
  relocate(diff_kdeg) %>% 
  filter(!is.na(NMD_50nt_rule)) %>% 
  ggplot(aes(x=diff_kdeg,
             fill=NMD_50nt_rule)) +
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
  geom_density(alpha=0.75,
               size=0.2) +
  xlim(-3.5,3.5) +
  scale_fill_manual(values=c("FALSE" = "#214D65",
                             "TRUE"="#E5BF86")) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_diff_kdeg.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR kdeg difference tx-vs-gene ------------------------------------------------------------
UPF1_NMDRHT_EZbakR_TEC_combined %>% 
  filter(condition == "UPF1_12h") %>% 
  filter(kdeg_tx_conclusion == "Stabilized") %>% 
  mutate(diff_kdeg = L2FC_kdeg_tx-L2FC_kdeg) %>% 
  relocate(diff_kdeg) %>% 
  mutate(gene_tx = paste0(gene_name," | ", transcript_id)) %>% 
  slice_min(diff_kdeg, n=10) %>% 
  arrange(L2FC_kdeg) %>% 
  mutate(gene_tx = fct_inorder(gene_tx)) %>% 
  pivot_longer(cols=c(L2FC_kdeg_tx, L2FC_kdeg, logFC),
               names_to = "name",
               values_to = "log2FC") %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other")) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "none",
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other",
                                            "lncRNA"
  )) %>% 
  ggplot(aes(x=name,
             y=gene_tx,
             fill=log2FC)) +
  theme(legend.position="right", 
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
  geom_tile() +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       # limits = c(-1,4),
                       na.value = "grey90") +
  # facet_wrap(~DTE_cluster, nrow=1) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggside::geom_ysidetile(aes(x = "NMD reason", yfill = NMD_tx_reason_simple))  +
  ggside::scale_yfill_manual(values=c("none" = "#214D65",
                                      "lncRNA" = "#76716E",
                                      "AS_NMD"="#E5BF85",
                                      "AS_NMD_UTR3" = "#C1AC65",
                                      "novel" = "#9A9A49",
                                      "uORF" = "#708831",
                                      "overl_uORF" = "#43761E",
                                      "multiple" = "gray30",
                                      "other" = "#D9DADA")) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_diff_kdeg_top10.pdf"),
       width = cw3,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### 12h N-AID-UPF1 EZbakR kdeg difference tx-vs-gene ------------------------------------------------------------
# Join with gene-level analyses for EZbakR parameter
NMDRHT.v1.2_MainTable_EZbakR_tx_all <- NMDRHT.v1.2_MainTable %>% 
  # filter(DTE_cluster == "up 1:early") %>% 
  # filter(!is.na(NMD_50nt_rule)) %>% 
  left_join(edgeR_DTE_NMDRHT_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h"))) %>%
  left_join(DESeq2_DGE_combined %>% 
              filter(condition_2 %in% c("UPF1_Nter_12h")) %>% 
              dplyr::select(gene_id, baseMean)) %>% 
  left_join(UPF1_NMDRHT_EZbakR_TEC_combined %>% 
            filter(condition == "UPF1_12h") %>% dplyr::select(transcript_id,
                                                              L2FC_kdeg_tx,
                                                              padj_kdeg_tx,
                                                              kdeg_tx_conclusion))

NMDRHT.v1.2_MainTable_EZbakR_tx_all %>% 
  mutate(EZbakR_detected = case_when(is.na(kdeg_tx_conclusion) ~ FALSE,
                                     !is.na(kdeg_tx_conclusion) ~ TRUE)) %>% 
  dplyr::count(NMD_50nt_rule, EZbakR_detected) %>% 
  group_by(NMD_50nt_rule) %>% 
  mutate(n_per = n/sum(n)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x=NMD_50nt_rule,
             y=n_per,
             fill=EZbakR_detected)) +
  theme(legend.position="right", 
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
             color="darkgray",
             linewidth=0.25) +
  geom_col(color="black",
           linewidth = 0.2) +
  geom_text(aes(label = case_when(n_per > 20 & EZbakR_detected == TRUE ~ paste0("n=",
                                                                                n,
                                                                                "\n(",
    round(n_per,0),
    "%)"),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle=90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 20 & EZbakR_detected == FALSE ~ paste0("n=",
                                                                                  n,
                                                                                  "\n(",
                                                                                  round(n_per,0),
                                                                                  "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            angle=90,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("TRUE" = "#963B5A",
                             "FALSE" = "lightgray")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(15, "mm")) 

ggsave(file.path("Plots", "Figure5", "Rev_1_F5_H_EZbakR_UPF1_12h_Barplot_Detected.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
