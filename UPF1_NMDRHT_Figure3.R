#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_Figure3
# Objective: Code for generating Panels of Figure 3 + Figure S3 for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

# Note upfront: relies on libraries and data loaded in by the main script "UPF1_NMDRHT_Analysis"

# Load DESeq2 DGE *GENCODE*-annotation-based data
DESeq2_DGE_combined <- read_csv("Resources/DESeq2_DGE_combined.csv")

###
## Rev_1 - Figure 3 - (A) bakR -----------------------------------------------------------
###

# Load UPF1_4SU_Fit_stan
load(file = "Resources/bakR/UPF1_4SU_Fit_stan.csv")

# Load Mechs_UPF1_combined
Mechs_UPF1_combined <- read_csv(file=file.path("Resources/bakR/Mechs_UPF1_combined.csv"))

### Mutation rates -----------------------------------------------------------

# Raw Mutation rate
UPF1_4SU_Fit_stan_raw_mutrate <- UPF1_4SU_Fit_stan$Data_lists$Fast_df %>% 
  dplyr::group_by(sample, type) %>%
  dplyr::summarise(mutrate = sum(TC*n)/sum(nT*n)) %>% 
  separate(sample,
           c("sample", "replicate"),
           sep=-2) %>% 
  mutate(replicate = parse_number(replicate)) %>% 
  mutate(sample = fct_rev(fct_relevel(sample,
                                      "Ctrl_0h_IAA_4SU",
                                      "Ctrl_12h_IAA_4SU",
                                      "Ctrl_24h_IAA_4SU",
                                      "Nter_0h_IAA",
                                      "Nter_0h_IAA_4SU",
                                      "Nter_12h_IAA_4SU",
                                      "Nter_24h_IAA_4SU")))

# Raw Mutation rate Plot
UPF1_4SU_Fit_stan_raw_mutrate %>% 
  group_by(sample) %>% 
  summarize(mean_mutrate = mean(mutrate)) %>% 
  ggplot(aes(y=sample,
             x="T-to-C\nrate",
             fill=(mean_mutrate*100))) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  scale_fill_viridis(option="rocket",
                     direction=1,
                     begin=0.2,
                     end=0.8) +
  geom_text(aes(label=round(mean_mutrate*100,1)),
            size = 5*0.36,
            hjust=-0.8) +
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
       fill="Raw T-to-C\nmutation rate (%)") +
  coord_fixed(ratio=1,
              clip="off") +
  # theme(legend.direction = "horizontal", legend.box = "vertical") +
  # guides(fill = guide_colourbar(order = 1,
  #                               title.position = "top",
  #                               label.position = "bottom",
  #                               label.hjust = 0.5,
  #                               label.vjust = 0.5,
  #                               label.theme = element_text(angle = 90, size = 6)),
  #        size = guide_legend(order = 2,
  #                            title.position = "top",
  #                            label.position = "bottom")) +
  # theme(axis.text.y=element_blank(),
#       axis.ticks.y = element_blank()) +
theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(4, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_A_UPF1_4SU_Raw_mutrates_heat.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Scatter plot ------------------------------------------------------------

# Prepare numbers
Mechs_UPF1_combined_conclusion <- Mechs_UPF1_combined %>% 
  dplyr::count(condition,Mech_conclusion, kdeg_conclusion, RNA_conclusion) %>% 
  arrange(desc(n))

# Plot
Mechs_UPF1_combined %>% 
  add_row(condition="Prediction", .before = 0) %>% 
  mutate(condition = fct_inorder(condition)) %>% 
  # dplyr::select(condition) %>% 
  # distinct(condition) %>% 
  ggplot(aes()) +
  geom_hline(yintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  geom_vline(xintercept = 0,
             color = "gray80",
             linewidth = 0.5,
             linetype = 1) +
  ggrastr::rasterise(geom_point(aes(size=-log10(mech_padj),
                                    x = L2FC_kdeg,
                                    y = L2FC_RNA,
                                    fill = log10(abs(mech_stat) + 1)*sign(mech_stat)),
                                shape=21,
                                stroke=0,
                                # color="black",
                                alpha=0.75), dpi=1200, dev = "cairo") +
  geom_abline(intercept = 0, slope = -1, linewidth=0.5, color="black", linetype = "dotted") +
  geom_label(data=(Mechs_UPF1_combined_conclusion %>% 
                     add_row(condition="Prediction", .before = 0) %>% 
                     mutate(condition = fct_inorder(condition))),
             aes(x = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ -6.5,
                               kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ 6.5,
                               kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ 6.5,
                               kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ -6.5),
                 y = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ 6.5,
                               kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ -6.5,
                               kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ 6.5,
                               kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ -6.5),
                 label = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ as.character(n),
                                   kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ as.character(n),
                                   kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ as.character(n),
                                   kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ as.character(n),
                                   TRUE ~ "")),
             size= 5*0.36,
             label.padding = unit(0.10, "lines")) +
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
  scale_fill_distiller(palette="BrBG") +
  # scale_color_viridis_c() + 
  scale_size(range = c(0.1, 1.25)) +
  xlab("L2FC(kdeg) from bakR") + 
  ylab("L2FC(RNA) from DESeq2") +
  labs(color = "Mechanism",
       title = "Mechanistic dissection | L2FC(kdeg) vs. L2FC(RNA)",
       subtitle = "Mechanism score: Synthesis driven = positive numbers; degradation driven = negative numbers.") +
  coord_fixed(ratio=1) + 
  xlim(-7.5, 7.5) +
  ylim(-7.5, 7.5) +
  facet_wrap(~condition, nrow=2) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_A_UPF1_combined_Mechanistic_Dissection.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### kdeg_ksyn-DGE Cluster -----------------------------------------------

gtf_gencode_df_short_DGE_cluster <- read_csv("Resources/GENCODE/gtf_gencode_df_short_DGE_cluster.csv") %>% 
  mutate(DGE_cluster = (fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early")))

gtf_gencode_df_short_DGE_cluster_Mech_4SU <- gtf_gencode_df_short_DGE_cluster %>% 
  left_join(Mechs_UPF1_combined,
            by=c("gene_id" = "XF"))

gtf_gencode_df_short_DGE_cluster_Mech_4SU %>% 
  filter(!is.na(condition)) %>% 
  filter(!DGE_cluster %in% c("not_expressed", "complex", "up 4:inverse", "down 4:inverse")) %>% 
  mutate(DGE_cluster = fct_rev(DGE_cluster)) %>% 
  pivot_longer(cols= c(L2FC_kdeg,
                       L2FC_ksyn),
               names_to = "rate",
               values_to = "log2FC") %>% 
  ggplot(aes(x=DGE_cluster,
             y=log2FC,
             fill=rate)) +
  geom_boxplot(outliers = FALSE,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.1) + 
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
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d")) +
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
  scale_y_continuous(limits = c(-4,4)) +
  facet_wrap(~condition, nrow=1) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_A_UPF1_combined_Mechanistic_Dissection_Cluster_boxplot.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Combine with NMD relevance bin Info -----------------------------------------------

# Load essential data for these analyses
load("Resources/NMD_relevance/DESeq2_DGE_combined_DGE_cluster_NMD_relevance.rds")

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin_Mech_4SU <- DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin %>% 
  left_join(Mechs_UPF1_combined,
            by=c("gene_id" = "XF"))

DESeq2_DGE_combined_DGE_cluster_NMD_relevance_complete_bin_Mech_4SU %>% 
  mutate(UpDown_NMD_bin = paste0(UpDown,"_",NMD_bin), .after=NMD_bin) %>% 
  filter(!is.na(condition)) %>% 
  mutate(UpDown_NMD_bin = fct_rev(fct_relevel(UpDown_NMD_bin,
                                          "up_(75,100]",
                                          "up_(50,75]",
                                          "up_(25,50]",
                                          "up_[0,25]",
                                          "down_[0,25]",
                                          "down_(25,50]",
                                          "down_(50,75]",
                                          "down_(75,100]"))) %>% 
  pivot_longer(cols= c(L2FC_kdeg,
                       L2FC_ksyn),
               names_to = "rate",
               values_to = "log2FC") %>% 
  ggplot(aes(x=UpDown_NMD_bin,
             y=log2FC,
             fill=rate)) +
  geom_boxplot(outliers = FALSE,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.1) + 
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
  scale_fill_manual(values=c("L2FC_kdeg" = "#80cdc1",
                             "L2FC_ksyn" = "#dfc27d")) +
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
  scale_y_continuous(limits = c(-4,4)) +
  facet_wrap(~condition, nrow=1) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_A_UPF1_combined_Mechanistic_Dissection_NMD_bin_boxplot.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### GO of synthesis-driven genes -----------------------------------------------

# Load necessary data
gostres_Mechs_UPF1_combined_Syn_GO_result_reduced <- read_csv("Resources/bakR/gostres_Mechs_UPF1_combined_Syn_GO_result_reduced.csv")

# Plot
gostres_Mechs_UPF1_combined_Syn_GO_result_reduced %>% 
  # mutate(parentTerm = paste0(parentTerm,"\n","(",parent_id,")")) %>% 
  mutate(query = fct_relevel(query,
                             "Increased Syn.",
                             "Decreased Syn."
  )) %>% 
  mutate(score=-log10(p_value)) %>% 
  group_by(query, parentTerm) %>% 
  mutate(mean_score=mean(score)) %>% 
  ungroup() %>% 
  # mutate(score = case_when(UpDown == "up" ~ score,
  #                          UpDown == "down" ~ -score)) %>%  
  # distinct(parentTerm, query, group_score) %>% 
  arrange((query),desc(mean_score)) %>% 
  mutate(parentTerm = fct_inorder(parentTerm)) %>% 
  # filter(parentTerm %in% c("response to stress",
  #                          "metabolic process",
  #                          "regulation of programmed cell death",
  #                          "biological regulation")) %>%
  ggplot(aes(x=score,
             y=parentTerm,
             fill=query)) +
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
  geom_bar(position = position_dodge(width=0.9,
                                     preserve = "total"),
           stat = "summary", 
           color="black",
           linewidth= 0.2,
           fun='mean',
           alpha=0.5) +
  geom_point(aes(group = query),
             shape=21,
             size=0.75,
             alpha=0.75,
             color="gray20",
             stroke=0.2,
             position = position_dodge(width=0.9, preserve = 'total')) +
  scale_fill_manual(values=c("Decreased Syn." = "#7F9494",
                             "Increased Syn." = "#BC6145")) +
  theme(axis.text.y = element_text(size = 5)) +
  xlim(0,15) +
  # theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  guides(fill=guide_legend(reverse=F)) +
  facet_wrap(~query,ncol=1, scales = "free_y") +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_A_Mechs_UPF1_combined_Syn_GO_result_reduced.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 3 - (B) HMD -----------------------------------------------------------
###

### Heatmap - RiboMinus --------------------------------------------------------------------- 

# For Histone-mediated decay (HMD) - Histone gene set (from https://www.genenames.org/data/genegroup/#!/group/864 retrieved on 27.09.2023)
Histone_HGNC <- read_delim("Resources/HMD/Histone_HGNC_group864.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) 

### HeatMap HMD -------------------------------------------------------------

Histone_HGNC %>% 
  dplyr::select(`Ensembl gene ID`) %>% 
  dplyr::rename("gene_id" = `Ensembl gene ID`) %>% 
  left_join(DESeq2_DGE_combined %>% 
              dplyr::filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                                 "HCT116_UPF1_AID_degradation_riboMinus_this_Study")) %>% 
              separate(gene_id, c("gene_id", "gene_id_version"))) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "control_48h", 
                                           "UPF1_Nter_0h", 
                                           "UPF1_Nter_2h", 
                                           "UPF1_Nter_4h", 
                                           "UPF1_Nter_8h",
                                           "UPF1_Nter_12h",
                                           "UPF1_Nter_24h", 
                                           "UPF1_Nter_48h",
                                           "UPF1_Nter_12h_riboMinus",
                                           "UPF1_Nter_48h_riboMinus"))) %>% 
  # dplyr::filter(gene_id %in% Histone_HGNC$`Ensembl gene ID`) %>% 
  arrange(desc(baseMean)) %>% 
  mutate(gene_id = fct_inorder(gene_id)) %>% 
  ggplot(aes(y=condition_2,
             x=gene_id,
             fill=log10(baseMean))) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako",
                       begin = 0.1,
                       end = 0.9) +
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
        axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  # guides(fill = guide_colourbar(barwidth = 0.5,
  #                               barheight = 2.5,
  #                               title = "log2FC",
  #                               direction = "vertical", 
  #                               title.position = "top",
  #                               label.position="right", 
  #                               label.theme = element_text(angle = 0,
  #                                                          size = 6))) +
  labs(x="Histone genes",
       y="",
       fill="log10\n(baseMean)") +
  coord_fixed(ratio=1) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(11*2+0.4, "mm"),
                   cols = unit(40, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_B_HMD_HeatMap_all_UPF1_BaseMean.pdf"),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

### Point log2FC plot -------------------------------------------------------

DESeq2_DGE_combined %>% 
  dplyr::filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_riboMinus_this_Study")) %>% 
  separate(gene_id, c("gene_id", "gene_id_version")) %>% 
  dplyr::filter(gene_id %in% Histone_HGNC$`Ensembl gene ID`) %>% 
  ggplot(aes(x=log10(baseMean),
             y=log2FoldChange,
             size=-log10(padj),
             fill=log2FoldChange)) +
  theme_classic() +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size = 6),
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
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_hline(yintercept = 0,
             linewidth = 0.5) +
  geom_point(shape=21,color="black", alpha=0.75) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  coord_cartesian(clip = 'off') +
  facet_wrap(~condition_2) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(35, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_B_HMD_Exp_L2FC_RiboMinus_UPF1_DGE.pdf"),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 


### H3-3A | H2AC18 | H1-2 ---------------------------------------------------

DESeq2_DGE_combined %>% 
  dplyr::filter(experimentSet %in% c("HCT116_UPF1_AID_degradation_this_Study",
                                     "HCT116_UPF1_AID_degradation_riboMinus_this_Study")) %>% 
  separate(gene_id, c("gene_id", "gene_id_version")) %>% 
  dplyr::filter(gene_id %in% Histone_HGNC$`Ensembl gene ID`) %>% 
  filter(gene_name %in% c("H1-2",
                          "H3-3A",
                          "H2AC18")) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "control_48h", 
                                           "UPF1_Nter_0h", 
                                           "UPF1_Nter_2h", 
                                           "UPF1_Nter_4h", 
                                           "UPF1_Nter_8h",
                                           "UPF1_Nter_12h",
                                           "UPF1_Nter_24h", 
                                           "UPF1_Nter_48h",
                                           "UPF1_Nter_12h_riboMinus",
                                           "UPF1_Nter_48h_riboMinus"))) %>% 
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
  # theme(axis.text.y=element_blank(),
  #       axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(10*2+0.4, "mm"),
                   cols = unit(3*2+0.4, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_B_HMD_examples_Dotmap_UPF1_DGE.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


###
## Rev_1 - Figure 3 - (C) Metabolic Parameters - Mukherjee 2017 -----------------------------------------------------------
###

# Load essential data for these analyses
load("Resources/NMD_relevance/DESeq2_DGE_combined_DGE_cluster_NMD_relevance.rds")

### Import data --------------------------------------------------------------------- 
# From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84722
# PMID: 27870833

# Import data and join with GENCODE global table
GSE84722_ST3 <- read_delim("Resources/External/Mukherjee2017_PMID_27870833/GSE84722_ST3.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>% 
  dplyr::select(-c(Annotation, Cluster, Host, Complex)) %>% 
  dplyr::rename("gene_id" = "Gene") %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  left_join(gtf_gencode_df_short_DGE_cluster_NMD_relevance_fill_bin %>% separate(gene_id, c("gene_id_plain", NA), remove = FALSE))

# Save as csv
GSE84722_ST3 %>% write_csv("Resources/External/Mukherjee2017_PMID_27870833/GSE84722_ST3.csv")

### Plot rates --------------------------------------------------------------

GSE84722_ST3 %>% 
  pivot_longer(cols=c(Syn, Proc, Deg, CytNuc, PolyCyt),
               names_to = "rate_name",
               values_to = "rate_value") %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  mutate(rate_name = fct_relevel(rate_name,
                                 "Syn",
                                 "Proc",
                                 "Deg",
                                 "CytNuc",
                                 "PolyCyt")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  ggplot(aes(x=DGE_cluster,
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
  coord_cartesian(ylim = c(-2,2)) +
  scale_fill_manual(values=c("Deg"="#80CDC1",
                             "Proc" = "#c17ecd",
                             "Syn"  = "#DFC27D",
                             "CytNuc" = "#cd7e8a",
                             "PolyCyt"="#D3CDBF")) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~rate_name, nrow=1) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(7*1.8, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_C_MetabolicRates_global_boxplot.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Plot others -------------------------------------------------------------
GSE84722_ST3 %>% 
  pivot_longer(cols=c(CytNuc, PolyCyt),
               names_to = "parameter_name",
               values_to = "parameter_value") %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  mutate(rate_name = fct_relevel(parameter_name,
                                 "CytNuc",
                                 "PolyCyt")) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "expressed",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early")) %>% 
  ggplot(aes(x=DGE_cluster,
             y=parameter_value,
             fill=parameter_name)) +
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
  # scale_fill_manual(values=c("Deg"="#80CDC1",
  #                            "Proc" = "#c17ecd",
  #                            "Syn"  = "#DFC27D")) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~parameter_name, scales = "free_y") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(14, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_C_MetabolicRates_CytoNuc_global_boxplot.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Stats -------------------------------------------------------------------
GSE84722_ST3 %>% 
  # pivot_longer(cols=c(Syn, Proc, Deg),
  #              names_to = "rate_name",
  #              values_to = "rate_value") %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "expressed",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early")) %>% 
  group_by(DGE_cluster) %>% 
  rstatix::get_summary_stats(Syn) 

# Prepare data
GSE84722_ST3_forStats <- GSE84722_ST3 %>% 
  # pivot_longer(cols=c(Syn, Proc, Deg),
  #              names_to = "rate_name",
  #              values_to = "rate_value") %>% 
  filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  filter(!is.na(DGE_cluster)) %>% 
  mutate(DGE_cluster = fct_relevel(DGE_cluster, 
                                   "up 1:early",
                                   "up 2:delayed",
                                   "up 3:late",
                                   "expressed",
                                   "down 3:late",
                                   "down 2:delayed",
                                   "down 1:early"))

# Expressed as control
GSE84722_ST3_stats <- GSE84722_ST3_forStats %>% rstatix::dunn_test(Syn ~ DGE_cluster) %>% 
  filter(group1 == "expressed" | group2 == "expressed") %>% 
  bind_rows(GSE84722_ST3_forStats %>% rstatix::dunn_test(Proc ~ DGE_cluster) %>% 
              filter(group1 == "expressed" | group2 == "expressed")) %>% 
  bind_rows(GSE84722_ST3_forStats %>% rstatix::dunn_test(Deg ~ DGE_cluster) %>% 
              filter(group1 == "expressed" | group2 == "expressed")) %>% 
  bind_rows(GSE84722_ST3_forStats %>% rstatix::dunn_test(CytNuc ~ DGE_cluster) %>% 
              filter(group1 == "expressed" | group2 == "expressed")) %>% 
  bind_rows(GSE84722_ST3_forStats %>% rstatix::dunn_test(PolyCyt ~ DGE_cluster) %>% 
              filter(group1 == "expressed" | group2 == "expressed"))

## Plots
GSE84722_ST3_stats %>% 
  mutate(p.adj = replace(p.adj, p.adj < 1e-100, 1e-100)) %>% 
  mutate(DGE_cluster = case_when(group1 == "expressed" ~ group2,
                                 group1 != "expressed" ~ group1)) %>% 
  add_row(DGE_cluster = "expressed", .y. = c("Syn",
                                             "Proc",
                                             "Deg",
                                             "CytNuc",
                                             "PolyCyt")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster, 
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  ggplot(aes(x=DGE_cluster,
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
                       limits = c(0, 100),
                       na.value = "grey80") +
  # theme(axis.text = element_blank(),
  #       axis.ticks = element_blank()) +
  labs(x="",
       y="",
       fill="-log10(padj)") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(7*1.8, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_C_MetabolicRates_CytoNuc_DunnTest_padj.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")



### NMD_bin -----------------------------------------------------------------
GSE84722_ST3 %>% 
  pivot_longer(cols=c(Syn, Proc, Deg, CytNuc, PolyCyt),
               names_to = "rate_name",
               values_to = "rate_value") %>% 
  # filter(!DGE_cluster %in% c("not_expressed", "up 4:inverse", "down 4:inverse")) %>% 
  # filter(!is.na(DGE_cluster)) %>% 
  mutate(rate_name = fct_relevel(rate_name,
                                 "Syn",
                                 "Proc",
                                 "Deg",
                                 "CytNuc",
                                 "PolyCyt")) %>% 
  filter(!is.na(NMD_bin)) %>% 
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
  ggplot(aes(x=NMD_bin_updown,
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
  coord_cartesian(ylim = c(-2,2)) +
  scale_fill_manual(values=c("Deg"="#80CDC1",
                             "Proc" = "#c17ecd",
                             "Syn"  = "#DFC27D",
                             "CytNuc" = "#cd7e8a",
                             "PolyCyt"="#D3CDBF")) +
  # scale_fill_manual(values=c("Deg"="#80CDC1",
  #                            "Proc" = "#c17ecd",
  #                            "Syn"  = "#DFC27D")) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~rate_name, nrow=1) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(9*1.8, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_C_MetabolicRates_global_NMD_bin_boxplot.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Stats NMD_bin -------------------------------------------------------------------
GSE84722_ST3 %>% 
  # pivot_longer(cols=c(Syn, Proc, Deg),
  #              names_to = "rate_name",
  #              values_to = "rate_value") %>% 
  filter(!is.na(NMD_bin)) %>% 
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
  group_by(NMD_bin_updown) %>% 
  rstatix::get_summary_stats(Deg) 

# Prepare data
GSE84722_ST3_forStats_NMDbin <- GSE84722_ST3 %>% 
  # pivot_longer(cols=c(Syn, Proc, Deg),
  #              names_to = "rate_name",
  #              values_to = "rate_value") %>% 
  filter(!is.na(NMD_bin)) %>% 
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
  )))

# Expressed as control
GSE84722_ST3_Stats_NMDbin <- GSE84722_ST3_forStats_NMDbin %>% rstatix::dunn_test(Syn ~ NMD_bin_updown) %>% 
  filter(group1 == "n.s._n.s." | group2 == "n.s._n.s.") %>% 
  bind_rows(GSE84722_ST3_forStats_NMDbin %>% rstatix::dunn_test(Proc ~ NMD_bin_updown) %>% 
              filter(group1 == "n.s._n.s." | group2 == "n.s._n.s.")) %>% 
  bind_rows(GSE84722_ST3_forStats_NMDbin %>% rstatix::dunn_test(Deg ~ NMD_bin_updown) %>% 
              filter(group1 == "n.s._n.s." | group2 == "n.s._n.s.")) %>% 
  bind_rows(GSE84722_ST3_forStats_NMDbin %>% rstatix::dunn_test(CytNuc ~ NMD_bin_updown) %>% 
              filter(group1 == "n.s._n.s." | group2 == "n.s._n.s.")) %>% 
  bind_rows(GSE84722_ST3_forStats_NMDbin %>% rstatix::dunn_test(PolyCyt ~ NMD_bin_updown) %>% 
              filter(group1 == "n.s._n.s." | group2 == "n.s._n.s."))

# Plot
GSE84722_ST3_Stats_NMDbin %>% 
  mutate(NMD_bin_updown = case_when(group1 == "n.s._n.s." ~ group2,
                                    group1 != "n.s._n.s." ~ group1)) %>% 
  add_row(NMD_bin_updown = "n.s._n.s.", .y. = c("Syn",
                                                "Proc",
                                                "Deg",
                                                "CytNuc",
                                                "PolyCyt")) %>% 
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
  ggplot(aes(x=NMD_bin_updown,
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
                       limits = c(0, 100),
                       na.value = "grey80") +
  # theme(axis.text = element_blank(),
  #       axis.ticks = element_blank()) +
  labs(x="",
       y="",
       fill="-log10(padj)") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(10, "mm"),
                   cols = unit(9*1.8, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_C_MetabolicRates_CytoNuc_NMD_bin_DunnTest_padj.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

###
## Rev_1 - Figure 3 - (D) Ribo-Seq -----------------------------------------------------------
###

### Ribo-Seq counts --------------------------------------------------------------------- 
ReadCounts_RiboSeq_Merged <- readxl::read_excel("Resources/Translation/ReadCounts_RiboSeq.xlsx",
                                        sheet = "Merged")

ReadCounts_RiboSeq_Merged %>% 
  dplyr::select(-c(P_sites_total_NMDRegHumanTxome)) %>% 
  pivot_longer(cols=-c(condition),
               names_to = "reads",
               values_to = "counts") %>% 
  mutate(condition=fct_rev(fct_relevel(condition,
                                       "N_AID_UPF1_0h",
                                       "N_AID_UPF1_12h"))) %>% 
  mutate(reads = fct_rev(fct_relevel(reads,
                                     "raw",
                                     "filtered",
                                     "mapped",
                                     "deduplicated",
                                     "P_sites_total_GENCODE"))) %>% 
  ggplot(aes(x=counts,
             y=condition,
             fill=reads)) +
  geom_col(position="dodge",
           color="black",
           alpha=0.75,
           linewidth=0.1) +
  geom_text(aes(x=10^2.5,
                label=paste(round(counts/10^6,0), " x 10^6")),
            size = 5*0.36,
            color="black",
            position = position_dodge(width = 0.9)
  ) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray60",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "right",
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x="Number of reads",
       y="",
       fill="") +
  guides(fill=guide_legend(reverse=T)) +
  scale_fill_brewer(palette = "Greens",
                    direction = -1) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(25, "mm"))

ggsave(filename = "Plots/Figure3/Rev_1_F3_D_GENCODE_RiboSeq_ReadCounts.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")


### P-site frame preference -------------------------------------------------

HCT_N_AID_UPF1_0h_IAA_P_sites_calcs <- read_delim("Resources/Translation/HCT_N_AID_UPF1_0h_IAA_P_sites_calcs", 
                                                  delim = "\t", escape_double = FALSE, 
                                                  trim_ws = TRUE) %>% 
  filter(comp == "nucl") %>% 
  arrange(read_length) %>% 
  mutate(condition = "N_AID_UPF1_0h")

HCT_N_AID_UPF1_12h_IAA_P_sites_calcs <- read_delim("Resources/Translation/HCT_N_AID_UPF1_12h_IAA_P_sites_calcs", 
                                                   delim = "\t", escape_double = FALSE, 
                                                   trim_ws = TRUE) %>% 
  filter(comp == "nucl") %>% 
  arrange(read_length) %>% 
  mutate(condition = "N_AID_UPF1_12h")

HCT_N_AID_UPF1_0h_IAA_P_sites_calcs %>% 
  bind_rows(HCT_N_AID_UPF1_12h_IAA_P_sites_calcs) %>% 
  mutate(condition=fct_rev(fct_relevel(condition,
                                       "N_AID_UPF1_0h",
                                       "N_AID_UPF1_12h"))) %>% 
  ggplot(aes(x=read_length,
             y=condition,
             fill=frame_preference)) +
  geom_tile(color = "black",
            lwd = 0.1,
            linetype = 1) +
  scale_x_continuous(limits = c(15,40),
                     breaks = c(15,20,25,30,35,40)) +
  scale_size(range = c(0.25, 1)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray60",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "gray60",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6), 
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "right",
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_viridis_c(option="rocket",
                       begin=0.05,
                       end=0.95) +
  force_panelsizes(rows = unit(8, "mm"),
                   cols = unit(25.1, "mm"))

ggsave(filename = "Plots/Figure3/Rev_1_F3_D_GENCODE_RiboSeq_FramePreference.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")


### Scatter plot ------------------------------------------------------------

# Load essential data for these analyses
combined_TE_results  <-  read_csv(file.path("Resources/Translation", "Rev_1_combined_TE_results.csv"))

# Count per class
combined_TE_results_stats <- combined_TE_results %>% 
  filter(!is.na(type)) %>% 
  dplyr::count(class, type)

# Scatter plot
combined_TE_results %>% 
  arrange(desc(class)) %>% 
  filter(!is.na(type)) %>% 
  ggplot(aes(x=L2FC_RNA,
             y=L2FC_Ribo)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=6),
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
  geom_abline(intercept = 0, slope = 1, linewidth=0.25, color="black", linetype = "dotted") +
  coord_fixed(ratio=1) +
  ggrastr::rasterise(geom_point(aes(fill=class,
                                    alpha=case_when(class == "Not Sig." ~ 0.25,
                                                    TRUE ~ 0.75)),
                                size=1,
                                shape=21,
                                stroke=0.1),
                     dpi=600, dev = "cairo") +
  # scale_color_manual(values=c("upregulated" = "#b2182b",
  #                            "expressed" = "gray30",
  #                            "downregulated" = "#2166ac")) +
  geom_label(data=combined_TE_results_stats,
             aes(x = case_when(class == "Concordant_up" ~ 8.5,
                               class == "Concordant_down" ~ -8.5,
                               class == "Discordant_up" ~ -8.5,
                               class == "Discordant_down" ~ 8.5),
                 y = case_when(class == "Concordant_up" ~ 8.5,
                               class == "Concordant_down" ~ -8.5,
                               class == "Discordant_up" ~ 8.5,
                               class == "Discordant_down" ~ -8.5),
                 label = case_when(class == "Concordant_up" ~ as.character(n),
                                   class == "Concordant_down" ~ as.character(n),
                                   class == "Discordant_up" ~ as.character(n),
                                   class == "Discordant_down" ~ as.character(n),
                                   TRUE ~ ""),
                 fill=class),
             size= 5*0.36,
             color="black",
             label.size = 0,
             alpha=0.25,
             show.legend = FALSE,
             label.padding = unit(0.10, "lines")) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  scale_alpha_identity() +
  scale_fill_manual(values=c("Concordant_up" = "#A7473A",
                             "Concordant_down" = "#4B5F6C",
                             "Discordant_up" = "#B09B37",
                             "Discordant_down" = "#7d9fc2",
                             "Not Sig." = "gray")) +
  facet_wrap(~type) +
  force_panelsizes(rows = unit(25, "mm"),
                   cols = unit(25, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_D_RiboSeq_DESeq2_TE_scatter.pdf"),
       width = cw2,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

### Discordant TE - Stabilized? ---------------------------------------------

# Load essential data for these analyses
GENCODE_v42_MainTable <- read_csv("Resources/GENCODE/GENCODE_v42_MainTable.csv")

GENCODE_v42_MainTable %>% 
  filter(!is.na(class_RiboSeq)) %>% 
  mutate(overall_conclusion = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ "Stabilized",
                                        kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ "Destabilized",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ "Increased Syn.",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ "Decreased Syn.",
                                        TRUE ~ "n.s.")) %>% 
  mutate(overall_conclusion = fct_relevel(overall_conclusion,
                                          "Stabilized",
                                          "Increased Syn.",
                                          "n.s.",
                                          "Decreased Syn.",
                                          "Destabilized")) %>% 
  mutate(class_RiboSeq = fct_relevel(class_RiboSeq,
                             "Concordant_up",
                             "Concordant_down",
                             "Not Sig.",
                             "Discordant_down",
                             "Discordant_up"
  )) %>% 
  dplyr::count(type,class_RiboSeq,overall_conclusion) %>% 
  group_by(type,class_RiboSeq) %>% 
  mutate(n_per = 100*round(n / sum(n), 5)) %>% 
  ggplot(aes(y=n_per,
             x=class_RiboSeq,
             fill=overall_conclusion)) +
  geom_col(position = "stack",
           color="black",
           linewidth=0.1) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.minor = element_blank(),
        strip.text=element_text(size=6, color="black"),
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
  # geom_vline(xintercept=0,
  #            linewidth=0.25,
  #            linetype="dashed") +
  scale_fill_manual(values=c("Stabilized" = "#8C6767",
                             "Increased Syn." = "#BC6145",
                             "n.s." = "gray80",
                             "Decreased Syn." = "#7F9494",
                             "Destabilized" = "#4F9AB7")) +
  labs(fill="bakR conclusion") +
  guides(size=guide_legend(nrow=2,byrow=TRUE),
         fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "right",
                                barwidth = unit(2, "mm"),
                                barheight = unit(8, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6))) +
  facet_wrap(~type,ncol=3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_D_RiboTE_bakR_conclusion.pdf"),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent") 

### Cluster -----------------------------------------------------------------

combined_TE_results %>% 
  dplyr::select(gene_id,L2FC_TE,L2FC_RNA,L2FC_Ribo,gene_name,type,DGE_cluster,class,NMD_bin) %>% 
  filter(!DGE_cluster %in% c("not_expressed", "complex", "up 4:inverse", "down 4:inverse")) %>% 
  mutate(DGE_cluster = fct_rev(fct_relevel(DGE_cluster,
                                           "up 1:early",
                                           "up 2:delayed",
                                           "up 3:late",
                                           "expressed",
                                           "down 3:late",
                                           "down 2:delayed",
                                           "down 1:early"))) %>% 
  pivot_longer(cols= c(L2FC_TE,L2FC_RNA,L2FC_Ribo),
               names_to = "parameter",
               values_to = "log2FC") %>% 
  ggplot(aes(x=DGE_cluster,
             y=log2FC,
             fill=parameter)) +
  geom_boxplot(outlier.shape = NA,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.1) + 
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
  scale_fill_manual(values=c("L2FC_RNA" = "#FAE093",
                             "L2FC_Ribo" = "#A7473A",
                             "L2FC_TE" = "#7D9FC2")) +
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
  scale_y_continuous(limits = c(-4,4)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(30, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_D_UPF1_RiboSeq_Cluster_boxplot.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


###  NMD bin ----------------------------------------------------------------

combined_TE_results %>% 
  filter(!is.na(NMD_bin)) %>% 
  dplyr::select(gene_id,L2FC_TE,L2FC_RNA,L2FC_Ribo,gene_name,type,DGE_cluster,class,UpDown,NMD_bin) %>% 
  mutate(UpDown_bin = paste0(UpDown,"_",NMD_bin), .after=NMD_bin) %>% 
  filter(!DGE_cluster %in% c("not_expressed", "complex")) %>% 
  pivot_longer(cols= c(L2FC_TE,L2FC_RNA,L2FC_Ribo),
               names_to = "parameter",
               values_to = "log2FC") %>% 
  mutate(UpDown_bin = fct_rev(fct_relevel(UpDown_bin,
                                          "up_(75,100]",
                                          "up_(50,75]",
                                          "up_(25,50]",
                                          "up_[0,25]",
                                          "down_[0,25]",
                                          "down_(25,50]",
                                          "down_(50,75]",
                                          "down_(75,100]"))) %>% 
  ggplot(aes(x=UpDown_bin,
             y=log2FC,
             fill=parameter)) +
  geom_boxplot(outlier.shape = NA,
               # varwidth = TRUE,
               fatten=2,
               linewidth=0.1) + 
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
  scale_fill_manual(values=c("L2FC_RNA" = "#FAE093",
                             "L2FC_Ribo" = "#A7473A",
                             "L2FC_TE" = "#7D9FC2")) +
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
  scale_y_continuous(limits = c(-4,4)) +
  force_panelsizes(rows = unit(15, "mm"),
                   cols = unit(30, "mm"))

ggsave(file.path("Plots/Figure3", "Rev_1_F3_D_UPF1_RiboSeq_NMD_bin_boxplot.pdf"),
       width = 30,
       height = 30,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
