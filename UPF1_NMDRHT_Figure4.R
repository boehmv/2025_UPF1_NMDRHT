#!/usr/bin/env Rscript

# Title: UPF1_Revision_Figure4
# Objective: Code for generating Panels of Figure 4 + Figure S4 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_Revision_Analysis"

###
## Rev_1 - Figure 4 - (A1) GENCODE -----------------------------------------------------------
###

# Comment: Plot A1 = Number of genes and transcript in GENCODE release 42 annotation per biotype

# Plot of gene/transcript level relevant biotypes
gtf_gencode_df_short %>% 
  filter(type == "gene") %>% 
  distinct(gene_id, .keep_all = TRUE) %>% 
  dplyr::select(-type) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other",
                          gene_type == "protein_coding" ~ "protein-coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  dplyr::count(type) %>% 
  mutate(level = "gene") %>% 
  bind_rows(gtf_gencode_df_short %>% 
              filter(type == "transcript") %>% 
              distinct(transcript_id, .keep_all = TRUE) %>% 
              dplyr::select(-type) %>% 
              mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                                      transcript_type == "protein_coding" ~ "protein-coding",
                                      transcript_type == "nonsense_mediated_decay" ~ "NMD",
                                      transcript_type == "lncRNA" ~ "lncRNA")) %>% 
              dplyr::count(type) %>% 
              mutate(level = "transcript")) %>% 
  mutate(type = fct_rev(fct_relevel(type,
                                    "protein-coding",
                                    "lncRNA",
                                    "NMD",
                                    "other"))) %>% 
  ggplot(aes(x=level,
             y=n,
             fill=type)) +
  geom_col(color="black",
           linewidth=0.2) +
  geom_text(aes(label = case_when(type %in% c("other", "lncRNA", "NMD") ~ as.character(n),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(!type %in% c("other", "lncRNA", "NMD") ~ as.character(n),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(values=c("protein-coding" = "#214D65",
                             "lncRNA" = "#76716E",
                             "NMD" = "#E5BF86",
                             "other" = "#D9DADA")) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
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
  labs(fill="GENCODE\nbiotype",
       y="number of genes\nor transcripts") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure4", "Rev_1_F4_A1_GENCODE_overview_level.pdf"),
       width = 10,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###
## Rev_1 - Figure 4 - (A2) GENCODE TSL -----------------------------------------------------------
###

gtf_gencode_df_short %>% 
  # filter(gene_id %in% (gtf_gencode_df_short_DGE_cluster %>% dplyr::filter(DGE_cluster != "not_expressed") %>% pull(gene_id))) %>% 
  filter(type == "transcript") %>% 
  dplyr::select(-type) %>% 
  # filter(transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA")) %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          transcript_type == "protein_coding" ~ "protein-coding",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "lncRNA" ~ "lncRNA")) %>% 
  filter(type != "other") %>% 
  mutate(transcript_support_level = as.character(transcript_support_level)) %>%
  mutate(type = fct_rev(fct_relevel(type,
                                    "protein-coding",
                                    "NMD",
                                    "lncRNA"
                                    ))) %>% 
  group_by(transcript_support_level, type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ggplot(aes(y=n_per,
             x=type,
             fill=transcript_support_level)) +
  geom_col(position="stack",
           color="black",
           linewidth=0.1) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
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
  scale_fill_brewer(palette = "Purples", direction=-1,
                    na.value="gray50") +
  geom_text(aes(label = case_when(n_per > 20 & !transcript_support_level %in% c("1", "2") ~ paste0(n,
                                                                                                   "\n(",
                                                                                                   round(n_per,0),
                                                                                                   "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(n_per > 20 & (transcript_support_level %in% c("1", "2") | is.na(transcript_support_level)) ~ paste0(n,
                                                                                                                                      "\n(",
                                                                                                                                      round(n_per,0),
                                                                                                                                      "%)"),
                                  TRUE ~ "")),
            color="white",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="GENCODE\ntranscript\nbiotype",
       y="Fraction of transripts (%)",
       fill="Transcript\nsupport\nlevel") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure4", "Rev_1_F4_A2_GENCODE_TSL.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A3) Longread Stats -----------------------------------------------------------
###
longread_QC <- TidyMultiqc::load_multiqc(tibble(path="Resources/Longread/QC_longread/BAM/multiqc_data/multiqc_data.json") %>% 
                                           pull(path),
                                         sections = "raw") %>% 
  mutate(sample=case_when(str_detect(metadata.sample_id, "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_1") ~ "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_1",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_2") ~ "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_2",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_1") ~ "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_1",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_2") ~ "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_2",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_0h_IAA_PacBio_1") ~ "HCT_N_AID_UPF1_0h_IAA_PacBio_1",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_0h_IAA_PacBio_2") ~ "HCT_N_AID_UPF1_0h_IAA_PacBio_2",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_12h_IAA_PacBio_1") ~ "HCT_N_AID_UPF1_12h_IAA_PacBio_1",
                          str_detect(metadata.sample_id, "HCT_N_AID_UPF1_12h_IAA_PacBio_2") ~ "HCT_N_AID_UPF1_12h_IAA_PacBio_2"
  )) 

longread_QC %>% 
  dplyr::select(sample,raw.nanostat.median_read_length_aligned,raw.nanostat.median_read_quality_aligned,raw.nanostat.number_of_reads_aligned) %>% 
  dplyr::rename("median_length" = "raw.nanostat.median_read_length_aligned",
                "median_quality" = "raw.nanostat.median_read_quality_aligned",
                "reads" = "raw.nanostat.number_of_reads_aligned",
  ) %>% 
  mutate(sample = fct_rev(fct_relevel(sample,
                                      "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_1",
                                      "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_2",
                                      "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_1",
                                      "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_2",
                                      "HCT_N_AID_UPF1_0h_IAA_PacBio_1",
                                      "HCT_N_AID_UPF1_0h_IAA_PacBio_2",
                                      "HCT_N_AID_UPF1_12h_IAA_PacBio_1",
                                      "HCT_N_AID_UPF1_12h_IAA_PacBio_2"))) %>% 
  mutate(condition = case_when(sample %in% c("HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_1",
                                             "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA_2") ~ "HCT_N_AID_UPF1_0h_IAA_ONT_dRNA",
                               sample %in% c("HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_1",
                                             "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA_2") ~ "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA",
                               sample %in% c("HCT_N_AID_UPF1_0h_IAA_PacBio_1",
                                             "HCT_N_AID_UPF1_0h_IAA_PacBio_2") ~ "HCT_N_AID_UPF1_0h_IAA_PacBio",
                               sample %in% c("HCT_N_AID_UPF1_12h_IAA_PacBio_1",
                                             "HCT_N_AID_UPF1_12h_IAA_PacBio_2") ~ "HCT_N_AID_UPF1_12h_IAA_PacBio")) %>% 
  pivot_longer(cols=-c(sample, condition),
               names_to="parameter",
               values_to="value") %>% 
  ggplot(aes(x=value,
             y=sample,
             fill=condition)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
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
  geom_col(color="black",
           linewidth=0.1) +
  geom_text(aes(label = value),
            color="black",
            size = 4*0.36,
            hjust = 1.25
  ) +
  scale_fill_manual(values = c("HCT_N_AID_UPF1_0h_IAA_ONT_dRNA" = "#A2D9F7",
                               "HCT_N_AID_UPF1_12h_IAA_ONT_dRNA" = "#66BAE8",
                               "HCT_N_AID_UPF1_0h_IAA_PacBio" = "#D2CDE7",
                               "HCT_N_AID_UPF1_12h_IAA_PacBio" = "#A597DB")) +
  guides(color="none",
         fill="none") +
  labs(x="",
       y="",
       fill="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~parameter,
             scales = "free_x") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots", "Figure4", "Rev_1_F4_A3_LongRead_statistics.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A4) NMDRHT SQANTI stats -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_1_datasources.rds")

NMDRHT_SQANTI3_GENCODE_fromGTF %>% 
  dplyr::count(structural_category) %>% 
  arrange(desc(n)) %>%
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  mutate(structural_category = fct_rev(structural_category)) %>% 
  ggplot(aes(x=n,
             y=structural_category,
             fill=structural_category)) +
  geom_col(color="black",
           linewidth = 0.2) +
  geom_text(aes(label = case_when(structural_category %in% c("novel_not_in_catalog", "novel_in_catalog", "intergenic") ~ paste0(
    round(n),
    ""),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(!structural_category %in% c("novel_not_in_catalog", "novel_in_catalog", "intergenic") ~ paste0(
    round(n),
    ""),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_viridis_d(direction = -1,
                       begin=0.1,
                       end=0.9) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
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
  labs(x="Number of transcript isoforms",
       y="",
       fill="SQANTI3\nstructural category") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A4_NMDRHT_SQANTI3_total_overview_absolute.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A5) NMDRHT UIC stats -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_2_datasources.rds")

NMDRHT_tracking_cleaned_sum %>% 
  mutate(UIC_bin = cut_width(UIC_total_support, 5, boundary = 0), .before=UIC_total_support) %>% 
  dplyr::count(UIC_bin) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  arrange(UIC_bin) %>% 
  mutate(UIC_bin_numbers = paste0(UIC_bin," n=", round(n), " (",round(n_per,0),"%)")) %>% 
  mutate(UIC_bin_numbers = fct_inorder(UIC_bin_numbers)) %>% 
  ggplot(aes(y=n_per,
             x="UIC",
             fill=UIC_bin_numbers)) +
  geom_col() +
  scale_fill_viridis_d(begin=0.1,
                       end=0.9,
                       option="inferno") +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
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
  labs(x="",
       y="",
       fill="Unique intron chain\ntotal support") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(10, "mm")) 

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A5_NMDRegHumanTxome_UIC_bin_fraction.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A6) SRSF2 UIC stats -----------------------------------------------------------
###

NMDRHT_tracking_cleaned_sum %>% 
  filter(gffcompare_gene_name == "SRSF2") %>% 
  mutate(mean_UIC_support = mean(UIC_total_support)) %>% 
  # mutate(name_class = paste0(gffcompare_transcript_name,"_",class_code)) %>% 
  ggplot(aes(x=gffcompare_transcript_name,
             y=UIC_total_support)) +
  geom_point(aes(fill=class_code),
             size=1.5,
             shape=21,
             alpha=0.75,
             color="black",
             stroke=0.2) +
  geom_hline(aes(yintercept=mean_UIC_support),
             color="darkred",
             linewidth=0.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label=paste("cutoff=",round(mean_UIC_support,0))),
            x=2,
            color="darkred",
            size = 5*0.36,
            y=20) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray80",
                                        linewidth = 0.1,
                                        linetype = 1),
        panel.grid.minor = element_blank(),
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
  scale_fill_brewer(palette="Dark2",
                    direction=1) +
  labs(x="",
       y="Total UIC support",
       fill="Gffcompare\nclass code") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) 

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A6_NMDRegHumanTxome_SRSF2_UIC.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A7) UIC filtering -----------------------------------------------------------
###

NMDRHT_tracking_filtered %>% 
  ggplot(aes(x=UIC_filter_level,
             y=UIC_total_support,
             fill=UIC_filter_level)) +
  geom_boxplot(outlier.shape = NA,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.1) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "gray80",
                                        linewidth = 0.1,
                                        linetype = 1),
        panel.grid.minor = element_blank(),
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
  scale_fill_manual(values=c("high" = "#A7473A",
                             "low" = "#565052")) +
  guides(fill="none") +
  labs(x="",
       y="Total UIC support",
       fill="Gffcompare\nclass code") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) 

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A7_NMDRegHumanTxome_UIC_support_filtered.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# How many loci are left?
NMDRHT_tracking_filtered %>% 
  distinct(gene_id,UIC_filter_level) %>%
  dplyr::count(UIC_filter_level)

# How many transcripts are left?
NMDRHT_tracking_filtered %>% 
  # distinct(NMDRHT_transcript_id) %>% 
  dplyr::count()

###
## Rev_1 - Figure 4 - (A8) NMDRHT final structural category -----------------------------------------------------------
###

NMDRHT.v1.1_tbl <- read_csv("Resources/NMDRHT/NMDRHT.v1.1_short.csv")

NMDRHT.v1.1_tbl %>% 
  dplyr::count(structural_category) %>% 
  mutate(n_per = n/sum(n))

NMDRHT.v1.1_tbl %>% 
  dplyr::count(structural_category) %>% 
  ggplot(aes(x="NMDRHT",
             y=n,
             fill=structural_category)) +
  geom_col(color="black",
           linewidth=0.2) +
  scale_fill_viridis_d(begin=0.1,
                       end=0.9) +
  geom_text(aes(label = case_when(structural_category %in% c("novel_not_in_catalog", "novel_in_catalog") ~ paste0(
    round(n),
    ""),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(structural_category %in% c("full-splice_match") ~ paste0(
    round(n),
    ""),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    position = position_stack(vjust = 0.5)
  ) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
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
  labs(fill="SQANTI3\ncategory",
       y="number of transcripts",
       y="") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A8_NMDRHT_SQANTI3_overview_level.pdf"),
       width = 10,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###
## Rev_1 - Figure 4 - (A9) NMDRHT ORFanage ORF type -----------------------------------------------------------
###

NMDRHT.v1.1_tbl %>% 
  filter(structural_category %in% c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>% 
  dplyr::count(ORF_type, structural_category) %>% 
  # mutate(structural_category = fct_rev(structural_category)) %>% 
  mutate(ORF_type = fct_rev(ORF_type)) %>%
  ggplot(aes(x=structural_category,
             y=n,
             fill=ORF_type)) +
  geom_col(color="black",
           linewidth=0.2) +
  # facet_wrap(~structural_category) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
  scale_fill_brewer(palette="Greens") +
  labs(x="",
       y="Number of transcripts",
       fill="ORFanage\npredicted\nORF type") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A9_NMDRegHumanTxome_structural_ORF.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A10) NMDRHT ORFanage NMD status -----------------------------------------------------------
###

NMDRHT.v1.1_tbl %>% 
  filter(structural_category %in% c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>% 
  dplyr::count(NMD_status, structural_category) %>% 
  mutate(NMD_status = case_when(is.na(NMD_status) ~ "NA",
                                TRUE ~ as.character(NMD_status))) %>% 
  mutate(NMD_status = fct_rev(fct_relevel(as.character(NMD_status),
                                          "NA",
                                          "FALSE",
                                          "TRUE"
  ))) %>% 
  mutate(structural_category = fct_rev(structural_category)) %>% 
  ggplot(aes(y=structural_category,
             x=n,
             fill=NMD_status)) +
  geom_col(color="black",
           linewidth=0.2) +
  # facet_wrap(~structural_category) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
  scale_fill_manual(values=c("FALSE" = "#287DAB",
                             "TRUE" = "#B09771",
                             "NA" = "#D9DADA")) +
  labs(y="",
       x="Number of transcripts",
       fill="ORFanage\npredicted\nNMD status") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A10_NMDRegHumanTxome_structural_NMD.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A11) NMDRHT ORFanage NMD status (fraction) -----------------------------------------------------------
###

NMDRHT.v1.1_tbl %>% 
  filter(structural_category %in% c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog")) %>% 
  dplyr::count(NMD_status, structural_category) %>% 
  group_by(structural_category) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ungroup() %>%
  mutate(NMD_status = case_when(is.na(NMD_status) ~ "NA",
                                TRUE ~ as.character(NMD_status))) %>% 
  mutate(NMD_status = fct_rev(fct_relevel(as.character(NMD_status),
                                          "NA",
                                          "FALSE",
                                          "TRUE"
  ))) %>% 
  mutate(structural_category = fct_rev(structural_category)) %>% 
  ggplot(aes(y=structural_category,
             x=n_per,
             fill=NMD_status)) +
  geom_col(color="black",
           linewidth=0.2) +
  geom_text(aes(label = case_when(n_per > 25 ~ paste0(n,
                                                      "\n(",
                                                      round(n_per,0),
                                                      "%)"),
                                  TRUE ~ "")),
            color="black",
            size = 4*0.36,
            position = position_stack(vjust = 0.5)
  ) +
  # facet_wrap(~structural_category) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
  scale_fill_manual(values=c("FALSE" = "#287DAB",
                             "TRUE" = "#B09771",
                             "NA" = "#D9DADA")) +
  labs(y="",
       x="Number of transcripts",
       fill="ORFanage\npredicted\nNMD status") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A11_NMDRegHumanTxome_structural_NMD_percent.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A12) ORFquant RiboTISH overlap -----------------------------------------------------------
###

# Load essential data for these analyses
HCT_N_AID_UPF1_ORFquant <- read_csv("Resources/ORFquant/HCT_N_AID_UPF1_ORFquant.csv")
HCT_N_AID_UPF1_RiboTISH_combined <- read_csv("Resources/RiboTISH/HCT_N_AID_UPF1_RiboTISH.csv")
HCT_N_AID_UPF1_ORFquant_RiboTish <- read_csv("Resources/NMDRHT/HCT_N_AID_UPF1_ORFquant_RiboTish.csv")

library(nVennR) # For overlaps with nVenns

# Select ORFquant and Ribo-TISH sets per condition
R1_ORFquant_ORFs_detected_0h <- HCT_N_AID_UPF1_ORFquant %>% 
  dplyr::filter(condition == "HCT_N_AID_UPF1_0h") %>% 
  pull(ORF_region)

R1_ORFquant_ORFs_detected_12h <- HCT_N_AID_UPF1_ORFquant %>% 
  dplyr::filter(condition == "HCT_N_AID_UPF1_12h") %>% 
  pull(ORF_region)

R1_RiboTISH_ORFs_detected_0h <- HCT_N_AID_UPF1_RiboTISH_combined %>% 
  dplyr::filter(condition == "HCT_N_AID_UPF1_0h") %>% 
  pull(ORF_region_ORFquantComp)

R1_RiboTISH_ORFs_detected_12h <- HCT_N_AID_UPF1_RiboTISH_combined %>% 
  dplyr::filter(condition == "HCT_N_AID_UPF1_12h") %>% 
  pull(ORF_region_ORFquantComp)

R1_ORFquant_RiboTISH_ORFs_detected_myNV <- plotVenn(list(R1_ORFquant_ORFs_detected_0h,
                                                         R1_ORFquant_ORFs_detected_12h,
                                                         R1_RiboTISH_ORFs_detected_0h,
                                                         R1_RiboTISH_ORFs_detected_12h), 
                                                    sNames=c("ORFquant_0h", 
                                                             "ORFquant_12h",
                                                             "RiboTISH_0h", 
                                                             "RiboTISH_12h"), showPlot = FALSE)

showSVG(R1_ORFquant_RiboTISH_ORFs_detected_myNV,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#420A68FF", "#932667FF", "#DD513AFF", "#FCA50AFF"))

# This output was saved and converted to PDF -> imported into CorelDraw for further visualization

###
## Rev_1 - Figure 4 - (A13) ORFquant RiboTISH ORF length -----------------------------------------------------------
###

HCT_N_AID_UPF1_ORFquant %>% 
  dplyr::select(ORF_region, ORF_length, condition) %>% 
  mutate(method = "ORFquant") %>% 
  bind_rows(HCT_N_AID_UPF1_RiboTISH_combined %>% 
              dplyr::select(ORF_region_ORFquantComp, AALen, condition) %>% 
              dplyr::rename("ORF_region" = "ORF_region_ORFquantComp",
                            "ORF_length" = "AALen") %>% 
              mutate(method = "RiboTISH")) %>% 
  left_join(HCT_N_AID_UPF1_ORFquant_RiboTish %>% dplyr::select(ORF_region, RiboTISH),
            relationship = "many-to-many") %>% 
  distinct() %>% 
  mutate(RiboTISH = case_when(is.na(RiboTISH) ~ FALSE,
                              TRUE ~ RiboTISH)) %>% 
  ggplot(aes(x=ORF_length,
             y=method,
             fill=condition)) +
  theme(legend.position="top", 
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
               # coef = 0,
               color="black",
               linewidth=0.1) +
  scale_fill_manual(values=c("HCT_N_AID_UPF1_0h" = "#65BB8C",
                             "HCT_N_AID_UPF1_12h" = "#00696D")) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_wrap(~RiboTISH) +
  labs(x="ORF length (aa)",
       y="") +
  force_panelsizes(rows = unit(10+0.4, "mm"),
                   cols = unit(40+0.4, "mm")) 

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A13_ORF_length_ORFprediction.pdf"),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")  

###
## Rev_1 - Figure 4 - (A14) ORFquant ORF status -----------------------------------------------------------
###

# Load updated NMDRHT information (contains Ribo-Seq ORFs and transcript properties)
NMDRHT.v1.2_tbl_length_GC_mfe <- as_tibble(read_csv(file.path("Resources/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe.csv")))

NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  mutate(NMD_tx_status = fct_rev(fct_relevel(NMD_tx_status,
                                             "coding",
                                             "mixed",
                                             "NMD",
                                             "predicted_coding",
                                             "predicted_NMD",
                                             "lncRNA"))) %>% 
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
                                            "lncRNA"
                                            
  )) %>% 
  dplyr::count(NMD_tx_status,NMD_tx_reason_simple) %>% 
  arrange(desc(n)) %>% 
  group_by(NMD_tx_status) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ggplot(aes(y=NMD_tx_status,
             x=n,
             fill=NMD_tx_reason_simple)) +
  geom_col(color="black",
           linewidth=0.1,
           position="stack") +
  theme(legend.position="right", 
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
  geom_text(aes(label = case_when(NMD_tx_reason_simple %in% c("lncRNA", "none") & n > 1500 ~ paste0(
    round(n,0)),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle = 90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(!NMD_tx_reason_simple %in% c("lncRNA", "none") & n > 1500 ~ paste0(
    round(n,0)),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    angle = 90,
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
  force_panelsizes(rows = unit(32, "mm"),
                   cols = unit(20, "mm")) +
  labs(x="Number of transcripts",
       y="ORFquant NMD status",
       fill="Reason") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(reverse=F))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A14_ORFquant_NMD_reason_stack.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A15) ORFquant ORF status fraction -----------------------------------------------------------
###

NMDRHT.v1.2_tbl_length_GC_mfe %>% 
  mutate(NMD_tx_status = fct_rev(fct_relevel(NMD_tx_status,
                                             "coding",
                                             "mixed",
                                             "NMD",
                                             "predicted_coding",
                                             "predicted_NMD",
                                             "lncRNA"))) %>% 
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
                                            "lncRNA"
                                            
  )) %>% 
  dplyr::count(NMD_tx_status,NMD_tx_reason_simple) %>% 
  arrange(desc(n)) %>% 
  group_by(NMD_tx_status) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ggplot(aes(y=NMD_tx_status,
             x=n_per,
             fill=NMD_tx_reason_simple)) +
  geom_col(color="black",
           linewidth=0.1,
           position="stack") +
  theme(legend.position="right", 
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
  geom_text(aes(label = case_when(NMD_tx_reason_simple %in% c("lncRNA", "none", "overl_uORF") & n_per > 10 ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="white",
    size = 4*0.36,
    angle = 90,
    position = position_stack(vjust = 0.5)
  ) +
  geom_text(aes(label = case_when(!NMD_tx_reason_simple %in% c("lncRNA", "none", "overl_uORF") & n_per > 10 ~ paste0(
    round(n_per,0),
    "%"),
    TRUE ~ "")),
    color="black",
    size = 4*0.36,
    angle = 90,
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
  force_panelsizes(rows = unit(32, "mm"),
                   cols = unit(20, "mm")) +
  labs(x="Number of transcripts",
       y="ORFquant NMD status",
       fill="Reason") +
  guides(fill=guide_legend(reverse=F))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A15_ORFquant_NMD_reason_fill.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A16) ORFquant alluvial plot -----------------------------------------------------------
###

library(ggalluvial)

# Load essential data for these analyses
NMDRHT.v1.2_tbl_NMD <- read_csv("Resources/NMDRHT/NMDRHT.v1.2_tbl_NMD.csv")

NMDRHT.v1.2_tbl_NMD %>% 
  mutate(structural_category_simple_2 = case_when(structural_category_simple %in% c("FSM", "NIC", "NNIC") ~ structural_category_simple,
                                                  TRUE  ~ "other")) %>% 
  filter(!NMD_tx_reason %in% c("none", "lncRNA")) %>% 
  mutate(NMD_tx_reason_simple = case_when(NMD_tx_reason %in% c("none",
                                                               "AS_NMD",
                                                               "AS_NMD_UTR3",
                                                               "novel",
                                                               "uORF",
                                                               "overl_uORF",
                                                               "lncRNA") ~ NMD_tx_reason,
                                          TRUE ~ "other_NMD")) %>% 
  dplyr::count(structural_category_simple_2,NMD_tx_status,NMD_tx_reason_simple, NMD_50nt_rule) %>% 
  ungroup() %>% 
  mutate(NMD_tx_status = (fct_relevel(NMD_tx_status,
                                      "mixed",
                                      "NMD",
                                      "predicted_NMD"))) %>% 
  mutate(NMD_50nt_rule = case_when(is.na(NMD_50nt_rule) ~ "lncRNA",
                                   TRUE ~ as.character(NMD_50nt_rule))) %>% 
  mutate(NMD_tx_reason_simple = fct_relevel(as_factor(NMD_tx_reason_simple),
                                            "AS_NMD",
                                            "AS_NMD_UTR3",
                                            "novel",
                                            "uORF",
                                            "overl_uORF",
                                            "other_NMD"
                                            
  )) %>% 
  ggplot(aes(
    axis1 = structural_category_simple_2,
    axis2 = NMD_tx_status,
    axis3 = NMD_tx_reason_simple,
    y = n)) +
  theme(legend.position="right", 
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
  geom_flow(aes(fill = NMD_tx_reason_simple),
            color="gray20",
            linewidth=0.2) +
  scale_fill_manual(values=c("none" = "#214D65",
                             "lncRNA" = "#76716E",
                             "AS_NMD"="#E5BF85",
                             "AS_NMD_UTR3" = "#C1AC65",
                             "novel" = "#9A9A49",
                             "uORF" = "#708831",
                             "overl_uORF" = "#43761E",
                             "other_NMD" = "#D9DADA")) +
  geom_stratum(color="gray20",
               linewidth=0.2) +
  geom_text(stat = "stratum",
            size = 4*0.36,
            aes(label = after_stat(stratum))) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(80, "mm"))


ggsave(file.path("Plots/Figure4", "Rev_1_F4_A16_ORFquant_NMD-only_reason_alluvial.pdf"),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A17) splam scores -----------------------------------------------------------
###

NMDRHT.v1.2_tbl_length_GC_mfe_splam <- read_csv(file.path("Resources/splam/NMDRHT/NMDRHT.v1.2_tbl_length_GC_mfe_splam.csv"))

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>% 
  pivot_longer(cols=c("avg_donor_score", "avg_acceptor_score"),
               names_to = "spliceSite",
               values_to = "score") %>% 
  mutate(spliceSite = fct_relevel(spliceSite,
                                  "avg_donor_score", "avg_acceptor_score")) %>% 
  ggplot(aes(y=score,
             x=UIC_total_support)) +
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
  geom_hex(binwidth = c(2, 0.1)) +
  facet_wrap(~spliceSite,ncol=1) +
  scale_fill_viridis(option="cividis") +
  # scale_fill_manual(values=c("avg_donor_score" = "#4C4C53",
  #                            "avg_acceptor_score" = "#7C6C65")) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A17_splam_NMDRHT_avgScore.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A18) splam scores filter-level & structure -----------------------------------------------------------
###

# Load 
load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_2_datasources.rds")

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>% 
  left_join(NMDRHT_tracking_filtered %>% 
              dplyr::select(NMDRHT_transcript_id, UIC_filter_level) %>% 
              dplyr::rename("transcript_id" = "NMDRHT_transcript_id")) %>% 
  pivot_longer(cols=c("avg_donor_score", "avg_acceptor_score"),
               names_to = "spliceSite",
               values_to = "score") %>% 
  mutate(spliceSite = fct_relevel(spliceSite,
                                  "avg_donor_score", "avg_acceptor_score")) %>% 
  filter(structural_category_simple %in% c("FSM", "NIC", "NNIC")) %>% 
  ggplot(aes(x=structural_category_simple,
             y=score,
             fill=spliceSite)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  facet_wrap(~spliceSite+UIC_filter_level) +
  scale_fill_manual(values=c("avg_donor_score"="#93c6e1",
                             "avg_acceptor_score"  = "#5f93ac")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(6, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A18_splam_NMDRHT_avgScore_StructuralCategory.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A19) GENCODE splam scores filter-level & structure -----------------------------------------------------------
###

GENCODE_splam_junction_score <- read_delim("Resources/splam/GENCODE/junction_score.bed", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = c("chr", "start", "end", "junc_name", "intron_num", "stramd", "donor_score", "acceptor_score", "transcript_id"), trim_ws = TRUE)  %>%
  separate_rows(transcript_id, sep = ",") %>% 
  group_by(transcript_id) %>% 
  mutate(avg_donor_score = mean(donor_score, na.rm = TRUE),
         avg_acceptor_score = mean(acceptor_score, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct(transcript_id, avg_donor_score, avg_acceptor_score) %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="transcript"))

# Plot
GENCODE_splam_junction_score %>% 
  pivot_longer(cols=c("avg_donor_score", "avg_acceptor_score"),
               names_to = "spliceSite",
               values_to = "score") %>% 
  mutate(spliceSite = fct_relevel(spliceSite,
                                  "avg_donor_score", "avg_acceptor_score")) %>% 
  ggplot(aes(x=as.factor(transcript_support_level),
             y=score,
             fill=spliceSite)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  facet_wrap(~spliceSite,ncol=1) +
  scale_fill_manual(values=c("avg_donor_score"="#93c6e1",
                             "avg_acceptor_score"  = "#5f93ac")) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A19_splam_GENCODE_avgScore_StructaralCategory.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A20) NMDRHT CAGE & polyA -----------------------------------------------------------
###

load("Resources/NMDRHT/2025-02-26_NMDRHT_checkpoint_3_datasources.rds")

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>%  
  left_join(NMDRHT_final_selection %>% 
              dplyr::select(NMDRHT_transcript_id, dist_to_CAGE_peak, dist_to_polyA_site) %>% 
              dplyr::rename("transcript_id"="NMDRHT_transcript_id")) %>%
  left_join(NMDRHT_tracking_filtered %>% 
              dplyr::select(NMDRHT_transcript_id, UIC_filter_level) %>% 
              dplyr::rename("transcript_id" = "NMDRHT_transcript_id")) %>% 
  filter(structural_category_simple %in% c("FSM", "NIC", "NNIC")) %>% 
  pivot_longer(cols=c(dist_to_CAGE_peak,dist_to_polyA_site),
               names_to = "tx_end_dist",
               values_to = "distance") %>% 
  ggplot(aes(x=structural_category_simple,
             y=distance,
             fill=tx_end_dist)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  scale_fill_manual(values=c("dist_to_CAGE_peak"="#f8b150",
                             "dist_to_polyA_site"  = "#c17d17")) +
  coord_cartesian(ylim = c(-100,9000)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~tx_end_dist+UIC_filter_level) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(6, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A20_NMDRHT_CAGE_polyA_dist_global_StructuralCategory.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A21) NMDRHT CAGE & polyA - Zoom -----------------------------------------------------------
###

NMDRHT.v1.2_tbl_length_GC_mfe_splam %>%  
  left_join(NMDRHT_final_selection %>% 
              dplyr::select(NMDRHT_transcript_id, dist_to_CAGE_peak, dist_to_polyA_site) %>% 
              dplyr::rename("transcript_id"="NMDRHT_transcript_id")) %>%
  left_join(NMDRHT_tracking_filtered %>% 
              dplyr::select(NMDRHT_transcript_id, UIC_filter_level) %>% 
              dplyr::rename("transcript_id" = "NMDRHT_transcript_id")) %>% 
  filter(structural_category_simple %in% c("FSM", "NIC", "NNIC")) %>% 
  pivot_longer(cols=c(dist_to_CAGE_peak,dist_to_polyA_site),
               names_to = "tx_end_dist",
               values_to = "distance") %>% 
  ggplot(aes(x=structural_category_simple,
             y=distance,
             fill=tx_end_dist)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  scale_fill_manual(values=c("dist_to_CAGE_peak"="#f8b150",
                             "dist_to_polyA_site"  = "#c17d17")) +
  coord_cartesian(ylim=c(-100,100)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~tx_end_dist+UIC_filter_level) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(6, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A21_NMDRHT_CAGE_polyA_dist_zoom_StructuralCategory.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A22) GENCODE_SQANTI3 zoom -----------------------------------------------------------
###

# Define chromomes to keep
selected_chromosomes = c(paste0("chr",seq_len(22)), "chrX", "chrY", "chrM")

# Import SQANTI3
GENCODE_SQANTI3 <-readr::read_tsv("Resources/GENCODE/SQANTI3/SQANTI3_gencode.v42_classification.txt") %>% 
  dplyr::rename("transcript_id" = "isoform",
                "seqname" = "chrom") %>%
  dplyr::filter(seqname %in% selected_chromosomes) %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="transcript"))

GENCODE_SQANTI3 %>% 
  mutate(transcript_support_level = as.factor(transcript_support_level)) %>% 
  pivot_longer(cols=c(dist_to_CAGE_peak,dist_to_polyA_site),
               names_to = "tx_end_dist",
               values_to = "distance") %>% 
  ggplot(aes(x=transcript_support_level,
             y=distance,
             fill=tx_end_dist)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  scale_fill_manual(values=c("dist_to_CAGE_peak"="#f8b150",
                             "dist_to_polyA_site"  = "#c17d17")) +
  # coord_cartesian(ylim = c(-3000,3000)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~tx_end_dist,ncol=1) +
  coord_cartesian(ylim=c(-100,100)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A22_GENCODE_CAGE_polyA_dist_zoom_TSL.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 4 - (A23) GENCODE_SQANTI3 -----------------------------------------------------------
###

GENCODE_SQANTI3 %>% 
  mutate(transcript_support_level = as.factor(transcript_support_level)) %>% 
  pivot_longer(cols=c(dist_to_CAGE_peak,dist_to_polyA_site),
               names_to = "tx_end_dist",
               values_to = "distance") %>% 
  ggplot(aes(x=transcript_support_level,
             y=distance,
             fill=tx_end_dist)) +
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
  geom_boxplot(outliers=FALSE,
               # position=position_dodge(0.7),
               # width=0.15,
               linewidth=0.1,
               color="black") +
  scale_fill_manual(values=c("dist_to_CAGE_peak"="#f8b150",
                             "dist_to_polyA_site"  = "#c17d17")) +
  coord_cartesian(ylim = c(-100,9000)) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~tx_end_dist,ncol=1) +
  # coord_cartesian(ylim=c(-100,3000)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(12, "mm"))

ggsave(file.path("Plots/Figure4", "Rev_1_F4_A23_GENCODE_CAGE_polyA_dist_global_TSL.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")