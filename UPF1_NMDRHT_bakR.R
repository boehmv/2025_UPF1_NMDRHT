#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_bakR
# Objective: Differential kinetic analysis via bakR for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(bakR)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(corrplot)
# library(ggside)

##
# Load data sources --------------------------------------------------------
##

# Load simplified gencode annotation
gtf_gencode_df_short <- read_csv("Resources/GENCODE/gencode.v42.gtf_df_short.csv")

# Load cB file from bam2bakR
# Comment: generated previously b< following the bam2bakR snakemake pipeline
UPF1_cB_csv <- read_csv("Resources/bakR/cB.csv.gz")

# How does it look?
head(UPF1_cB_csv)

# Which samples do we have?
print(unique(UPF1_cB_csv$sample))

# Load Metadata file
UPF1_metadf <- read_delim("Resources/bakR/metadf", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% 
  column_to_rownames(var="X1") %>% 
  dplyr::rename(tl = "X2",
         Exp_ID = "X3")

# Assemble bakR dataset
UPF1_bakRData <- bakRData(UPF1_cB_csv, UPF1_metadf)

##
# Fast fit analysis -------------------------------------------------------
##

# Fit
UPF1_4SU_Fit <- bakRFit(UPF1_bakRData)

# Do the QC -> access plots if needed
UPF1_QC_checks <- QC_checks(UPF1_4SU_Fit)

# Use alternative Stan-based pnew estimation
UPF1_4SU_Fit_stan <- bakRFit(UPF1_4SU_Fit, FastRerun = TRUE, StanRateEst = TRUE)

# Do the QC -> access plots if needed
UPF1_QC_checks_stan <- QC_checks(UPF1_4SU_Fit_stan)

# QC plot 1
UPF1_QC_checks_stan$raw_mutrates

ggsave(file.path("Resources/bakR", "UPF1_4SU_Raw_mutrates.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# QC plot 2
UPF1_QC_checks_stan$conversion_rates

ggsave(file.path("Resources/bakR", "UPF1_4SU_Conversion_rates.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Distribution of inferred mutation rates
ggplot(UPF1_4SU_Fit_stan$Fast_Fit$Fn_Estimates, aes(x = logit_fn, color = as.factor(sample))) + 
  geom_density() + 
  theme_classic() +
  scale_color_viridis_d() + 
  xlab("logit(fn) estimates") + 
  ylab("density") +
  labs(title = "Distribution of inferred mutation rates",
       subtitle = "logit x-axis scale")

ggsave(file.path("Resources/bakR", "UPF1_4SU_mutation_rates_distribution.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Using function from corrplot package
corrplot.mixed(UPF1_QC_checks_stan$correlation_matrix, 
               upper = "square", lower = "number", 
               addgrid.col = "black", tl.col = "black")


## Extract L2FC(kdeg) and padj
UPF1_4SU_L2FC_df <- UPF1_4SU_Fit_stan$Fast_Fit$Effects_df

## Add significance conclusion
UPF1_4SU_L2FC_df <- UPF1_4SU_L2FC_df %>% dplyr::mutate(conclusion = ifelse(padj < 0.05, ifelse(L2FC_kdeg < 0, "Stabilized", "Destabilized"), "Not Sig."))

# Add gene names and type | from GENCODE annotation
UPF1_4SU_L2FC_df_anot <- UPF1_4SU_L2FC_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id"))

# Plot volcano per condition
ggplot2::ggplot(UPF1_4SU_L2FC_df_anot, ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion )) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::theme_classic() +
  ggplot2::ylab(bquote(-log[10](p[adj]))) +
  ggplot2::xlab(bquote(L2FC(k[deg]))) +
  facet_wrap(~Exp_ID)

# Count up/down per condition
UPF1_4SU_L2FC_df_anot %>% 
  group_by(Exp_ID, conclusion) %>% 
  summarize(n = n())

# Export fast fit to CSV file
UPF1_4SU_L2FC_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/UPF1_4SU_L2FC_df_anot.csv"))

# Load in fast fit to CSV file
UPF1_4SU_L2FC_df_anot <- read_csv(file=file.path("Resources/bakR/UPF1_4SU_L2FC_df_anot.csv"))

##
# Hybrid fit analysis -----------------------------------------------------
##

# Load options that will make running models more efficient
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Run Full model (This might take ~10-30 minutes to run)
UPF1_4SU_Fit_stan <- bakRFit(UPF1_4SU_Fit_stan, HybridFit = TRUE)

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

ggsave(file.path("Resources/bakR", "UPF1_4SU_Raw_mutrates_heat.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Do the QC -> access plots if needed
UPF1_QC_checks_stan <- QC_checks(UPF1_4SU_Fit_stan)

# QC plot 1
UPF1_QC_checks_stan$raw_mutrates 

ggsave(file.path("Resources/bakR", "UPF1_4SU_Raw_mutrates.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# QC plot 2
UPF1_QC_checks_stan$conversion_rates

ggsave(file.path("Resources/bakR", "UPF1_4SU_Conversion_rates.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Distribution of inferred mutation rates
ggplot(UPF1_4SU_Fit_stan$Fast_Fit$Fn_Estimates, aes(x = logit_fn, color = as.factor(sample))) + 
  geom_density() + 
  theme_classic() +
  scale_color_viridis_d() + 
  xlab("logit(fn) estimates") + 
  ylab("density") +
  labs(title = "Distribution of inferred mutation rates",
       subtitle = "logit x-axis scale")

ggsave(file.path("Resources/bakR", "UPF1_4SU_mutation_rates_distribution.pdf"),
       width = 15,
       height = 10,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# PCA plot
FnPCA2(UPF1_4SU_Fit_stan, Model = "Hybrid")

ggsave(file.path("Resources/bakR", "UPF1_4SU_PCA_hybrid_fit.pdf"),
       width = 20,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

## Extract L2FC(kdeg) and padj
UPF1_4SU_L2FC_hybrid_stan_df <- UPF1_4SU_Fit_stan$Hybrid_Fit$Effects_df

## Add significance ID
UPF1_4SU_L2FC_hybrid_stan_df <- UPF1_4SU_L2FC_hybrid_stan_df %>% dplyr::mutate(conclusion = ifelse(padj < 0.05, ifelse(L2FC_kdeg < 0, "Stabilized", "Destabilized"), "Not Sig."))

# Add gene names and type
UPF1_4SU_L2FC_hybrid_stan_df_anot <- UPF1_4SU_L2FC_hybrid_stan_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id"))

# Plot volcano per condition
ggplot2::ggplot(UPF1_4SU_L2FC_hybrid_stan_df_anot, ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion )) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::theme_classic() +
  ggplot2::ylab(bquote(-log[10](p[adj]))) +
  ggplot2::xlab(bquote(L2FC(k[deg]))) +
  facet_wrap(~Exp_ID)

# Count up/down per condition
UPF1_4SU_L2FC_hybrid_stan_df_anot %>% 
  group_by(Exp_ID, conclusion) %>% 
  summarize(n = n())

# Export hybrid fit to CSV file
UPF1_4SU_L2FC_hybrid_stan_df_anot %>% 
  write_csv(file=file.path("/srv/2023_UPF1_4SU/bam2bakR/UPF1_4SU_L2FC_hybrid_stan_df_anot.csv"))

# Load hybrid fit from CSV file
UPF1_4SU_L2FC_hybrid_stan_df_anot <- read_csv(file=file.path("Resources/bakR/UPF1_4SU_L2FC_hybrid_stan_df_anot.csv"))

# Save Fit 
save(UPF1_4SU_Fit_stan, file = "Resources/bakR/UPF1_4SU_Fit_stan.csv")

# Load Fit if needed
load(file = "Resources/bakR/UPF1_4SU_Fit_stan.csv")

##
# Advanced analyses -------------------------------------------------------
##

# Get the count matrix from bakR
Counts <- UPF1_4SU_Fit_stan$Data_lists$Count_Matrix

# Define conditions -> needed for DESeq2
colData <- data.frame(conditions = c("Ctrl_0h",
                                     "Ctrl_0h",
                                     "Ctrl_0h",
                                     "Ctrl_12h",
                                     "Ctrl_12h",
                                     "Ctrl_12h",
                                     "Ctrl_24h",
                                     "Ctrl_24h",
                                     "Ctrl_24h",
                                     "Nter_0h",
                                     "Nter_0h",
                                     "Nter_0h",
                                     "Nter_0h",
                                     "Nter_0h",
                                     "Nter_0h",
                                     "Nter_12h",
                                     "Nter_12h",
                                     "Nter_12h",
                                     "Nter_24h",
                                     "Nter_24h",
                                     "Nter_24h"))

# Make the colData input for DESeq2
rownames(colData) <- colnames(Counts)

# Make DESeq2 data object
dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = colData,
                              design = ~conditions)

# Fit DESeq2 model
ddso <- DESeq(dds)

## Mechanism 12h_Ctrl ------------------------------------------------------

# Extract results of experimental vs. reference comparison
reso_Ctrl_12h_vs_0h <- results(ddso, contrast = c("conditions", "Ctrl_12h", "Ctrl_0h"))

# Look at the column names of the reso object
colnames(as.data.frame(reso_Ctrl_12h_vs_0h))

# Convert to data frame
reso_Ctrl_12h_vs_0h <- as.data.frame(reso_Ctrl_12h_vs_0h)

# Make data frame
DE_df_Ctrl_12h_vs_0h <- data.frame(XF = row.names(reso_Ctrl_12h_vs_0h),
                             L2FC_RNA = reso_Ctrl_12h_vs_0h$log2FoldChange,
                             DE_score = reso_Ctrl_12h_vs_0h$stat,
                             DE_se = reso_Ctrl_12h_vs_0h$lfcSE,
                             DE_pval = reso_Ctrl_12h_vs_0h$pval,
                             DE_padj = reso_Ctrl_12h_vs_0h$padj)

# mechanistic dissection
Mechs_Ctrl_12h <- DissectMechanism(UPF1_4SU_Fit_stan, 
                                  DE_df_Ctrl_12h_vs_0h,
                                  bakRModel = "Hybrid",
                                  Exp_ID = 2,
                                  sims = 1000000)

# Extract data frame
Mechs_Ctrl_12h_df <- Mechs_Ctrl_12h$Mechanism_df 

# Add gene names and type
Mechs_Ctrl_12h_df_anot <- Mechs_Ctrl_12h_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id")) %>% 
  mutate(condition="Ctrl_12h")

# Export mechanism analysis from hybrid fit to CSV file
Mechs_Ctrl_12h_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_12h_Ctrl_df_hybrid_stan_anot.csv"))

## Mechanism 24h_Ctrl ------------------------------------------------------

# Extract results of experimental vs. reference comparison
reso_Ctrl_24h_vs_0h <- results(ddso, contrast = c("conditions", "Ctrl_24h", "Ctrl_0h"))

# Look at the column names of the reso object
colnames(as.data.frame(reso_Ctrl_24h_vs_0h))

# Convert to data frame
reso_Ctrl_24h_vs_0h <- as.data.frame(reso_Ctrl_24h_vs_0h)

# Make data frame
DE_df_Ctrl_24h_vs_0h <- data.frame(XF = row.names(reso_Ctrl_24h_vs_0h),
                                   L2FC_RNA = reso_Ctrl_24h_vs_0h$log2FoldChange,
                                   DE_score = reso_Ctrl_24h_vs_0h$stat,
                                   DE_se = reso_Ctrl_24h_vs_0h$lfcSE,
                                   DE_pval = reso_Ctrl_24h_vs_0h$pval,
                                   DE_padj = reso_Ctrl_24h_vs_0h$padj)

# mechanistic dissection
Mechs_Ctrl_24h <- DissectMechanism(UPF1_4SU_Fit_stan, 
                                   DE_df_Ctrl_24h_vs_0h,
                                   bakRModel = "Hybrid",
                                   Exp_ID = 3,
                                   sims = 1000000)

# Extract data frame
Mechs_Ctrl_24h_df <- Mechs_Ctrl_24h$Mechanism_df 

# Add gene names and type
Mechs_Ctrl_24h_df_anot <- Mechs_Ctrl_24h_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id"))  %>% 
  mutate(condition="Ctrl_24h")

# Export mechanism analysis from hybrid fit to CSV file
Mechs_Ctrl_24h_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_24h_Ctrl_df_hybrid_stan_anot.csv"))

## Mechanism 0h_UPF1 ------------------------------------------------------

# Extract results of experimental vs. reference comparison
reso_0h_vs_0h <- results(ddso, contrast = c("conditions", "Nter_0h", "Ctrl_0h"))

# Look at the column names of the reso object
colnames(as.data.frame(reso_0h_vs_0h))

# Convert to data frame
reso_0h_vs_0h <- as.data.frame(reso_0h_vs_0h)

# Make data frame
DE_df_0h_vs_0h <- data.frame(XF = row.names(reso_0h_vs_0h),
                              L2FC_RNA = reso_0h_vs_0h$log2FoldChange,
                              DE_score = reso_0h_vs_0h$stat,
                              DE_se = reso_0h_vs_0h$lfcSE,
                              DE_pval = reso_0h_vs_0h$pval,
                              DE_padj = reso_0h_vs_0h$padj)

# mechanistic dissection
Mechs_0h_UPF1 <- DissectMechanism(UPF1_4SU_Fit_stan, 
                                   DE_df_0h_vs_0h,
                                   bakRModel = "Hybrid",
                                   Exp_ID = 4,
                                   sims = 1000000)

# Extract data frame
Mechs_0h_UPF1_df <- Mechs_0h_UPF1$Mechanism_df 

# Add gene names and type
Mechs_0h_UPF1_df_anot <- Mechs_0h_UPF1_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id")) %>% 
  mutate(condition="Nter_0h")

# Export mechanism analysis from hybrid fit to CSV file
Mechs_0h_UPF1_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_0h_UPF1_df_hybrid_stan_anot.csv"))

## Mechanism 12h_UPF1 ------------------------------------------------------

# Extract results of experimental vs. reference comparison
reso_12h_vs_0h <- results(ddso, contrast = c("conditions", "Nter_12h", "Ctrl_0h"))

# Look at the column names of the reso object
colnames(as.data.frame(reso_12h_vs_0h))

# Convert to data frame
reso_12h_vs_0h <- as.data.frame(reso_12h_vs_0h)

# Make data frame
DE_df_12h_vs_0h <- data.frame(XF = row.names(reso_12h_vs_0h),
                    L2FC_RNA = reso_12h_vs_0h$log2FoldChange,
                    DE_score = reso_12h_vs_0h$stat,
                    DE_se = reso_12h_vs_0h$lfcSE,
                    DE_pval = reso_12h_vs_0h$pval,
                    DE_padj = reso_12h_vs_0h$padj)

# mechanistic dissection
Mechs_12h_UPF1 <- DissectMechanism(UPF1_4SU_Fit_stan, 
                                   DE_df_12h_vs_0h,
                          bakRModel = "Hybrid",
                          Exp_ID = 5,
                          sims = 1000000)

# Extract data frame
Mechs_12h_UPF1_df <- Mechs_12h_UPF1$Mechanism_df 

# Add gene names and type
Mechs_12h_UPF1_df_anot <- Mechs_12h_UPF1_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id")) %>% 
  mutate(condition="Nter_12h")

# Export mechanism analysis from hybrid fit to CSV file
Mechs_12h_UPF1_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_12h_UPF1_df_hybrid_stan_anot.csv"))

## Mechanism 24h_UPF1 ------------------------------------------------------

# Extract results of experimental vs. reference comparison
reso_24h_vs_0h <- results(ddso, contrast = c("conditions", "Nter_24h", "Ctrl_0h"))

# Look at the column names of the reso object
colnames(as.data.frame(reso_24h_vs_0h))

# Convert to data frame
reso_24h_vs_0h <- as.data.frame(reso_24h_vs_0h)

# Make data frame
DE_df_24h_vs_0h <- data.frame(XF = row.names(reso_24h_vs_0h),
                              L2FC_RNA = reso_24h_vs_0h$log2FoldChange,
                              DE_score = reso_24h_vs_0h$stat,
                              DE_se = reso_24h_vs_0h$lfcSE,
                              DE_pval = reso_24h_vs_0h$pval,
                              DE_padj = reso_24h_vs_0h$padj)

# mechanistic dissection
Mechs_24h_UPF1 <- DissectMechanism(UPF1_4SU_Fit_stan, 
                                   DE_df_24h_vs_0h,
                                   bakRModel = "Hybrid",
                                   Exp_ID = 6,
                                   sims = 1000000)

# Extract data frame
Mechs_24h_UPF1_df <- Mechs_24h_UPF1$Mechanism_df 

# Add gene names and type
Mechs_24h_UPF1_df_anot <- Mechs_24h_UPF1_df %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="gene"),
            by = c("XF" = "gene_id")) %>% 
  mutate(condition="Nter_24h")

# Export mechanism analysis from hybrid fit to CSV file
Mechs_24h_UPF1_df_anot %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_24h_UPF1_df_hybrid_stan_anot.csv"))

## Combined Mechanism ------------------------------------------------------
Mechs_UPF1_combined <- bind_rows(Mechs_Ctrl_12h_df_anot,
          Mechs_Ctrl_24h_df_anot,
          Mechs_0h_UPF1_df_anot,
          Mechs_12h_UPF1_df_anot,
          Mechs_24h_UPF1_df_anot) %>% 
  mutate(Mech_score = log10(abs(mech_stat) + 1)*sign(mech_stat)) %>% 
  mutate(kdeg_conclusion = ifelse(bakR_padj < 0.01, case_when(L2FC_kdeg < -1 ~ "Stabilized",
                                                              L2FC_kdeg > 1 ~ "Destabilized",
                                                              TRUE ~ "Not Sig."), "Not Sig.")) %>% 
  mutate(RNA_conclusion = ifelse(DE_padj < 0.01, case_when(L2FC_RNA < -1 ~ "Downregulated",
                                                           L2FC_RNA > 1 ~ "Upregulated",
                                                           TRUE ~ "Not Sig."), "Not Sig.")) %>% 
  mutate(Mech_conclusion = ifelse(mech_padj < 0.01, ifelse(Mech_score < 0, "Degradation", "Synthesis"), "Not Sig."))

# Export mechanism analysis from hybrid fit to CSV file
Mechs_UPF1_combined %>% 
  write_csv(file=file.path("Resources/bakR/Mechs_UPF1_combined.csv"))

# Read file again if needed
Mechs_UPF1_combined <- read_csv(file=file.path("Resources/bakR/Mechs_UPF1_combined.csv"))

# Define general conclusion
Mechs_UPF1_combined_generalConclusion <- Mechs_UPF1_combined %>% 
  mutate(general_conclusion = case_when(kdeg_conclusion == "Stabilized" & RNA_conclusion == "Upregulated" & Mech_conclusion == "Degradation" ~ "Stabilized",
                                        kdeg_conclusion == "Destabilized" & RNA_conclusion == "Downregulated" & Mech_conclusion == "Degradation" ~ "Destabilized",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Upregulated" & Mech_conclusion == "Synthesis" ~ "Increased Syn.",
                                        kdeg_conclusion == "Not Sig." & RNA_conclusion == "Downregulated" & Mech_conclusion == "Synthesis" ~ "Decreased Syn.",
                                        is.na(kdeg_conclusion) ~ "NA",
                                        TRUE ~ "Not Sig."))


# GO of synthesis-driven genes -----------------------------------------------------

# Obtain background of expressed genes from degradation/recovery HCT116 data
Rev_1_F2_GO_DGE_bg_full <- DESeq2_DGE_combined %>% 
  filter(experimentSet %in% c(
    "HCT116_UPF1_AID_degradation_this_Study",
    "HCT116_UPF1_AID_recovery_this_Study")) %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Make lists
Mechs_UPF1_combined_Syn_GO_list <- split(Mechs_UPF1_combined_generalConclusion %>% 
                                                                 filter(general_conclusion %in% c("Decreased Syn.", "Increased Syn.")) %>% 
                                                                 filter(condition == "Nter_12h") %>% 
                                                                 dplyr::rename("gene_id" = "XF") %>% 
                                                                 arrange((ksyn_padj)) %>% 
                                                                 distinct(gene_id, general_conclusion) %>% 
                                                                 separate(gene_id, c("gene_id", "Version")) %>% 
                                                                 pull(gene_id),
                                                               Mechs_UPF1_combined_generalConclusion %>% 
                                                                 filter(general_conclusion %in% c("Decreased Syn.", "Increased Syn.")) %>% 
                                                                 filter(condition == "Nter_12h") %>% 
                                                                 dplyr::rename("gene_id" = "XF") %>% 
                                                                 arrange((ksyn_padj)) %>% 
                                                                 distinct(gene_id, general_conclusion) %>% 
                                                                 separate(gene_id, c("gene_id", "Version")) %>% 
                                                                 pull(general_conclusion))

# Perform GO analysis
# Ordered query set to true (more significant genes contribute more!)
gostres_Mechs_UPF1_combined_Syn_GO <- gprofiler2::gost(query = Mechs_UPF1_combined_Syn_GO_list,
                                                       custom_bg = Rev_1_F2_GO_DGE_bg_full,
                                                       # multi_query = TRUE,
                                                       sources = "GO:BP",
                                                       ordered_query = TRUE,
                                                       domain_scope = "custom",
                                                       organism = "hsapiens",
                                                       correction_method = c("gSCS"),
                                                       evcodes = TRUE,
                                                       # as_short_link = TRUE,
                                                       significant = TRUE)

# As tibble
gostres_Mechs_UPF1_combined_Syn_GO_result <- gostres_Mechs_UPF1_combined_Syn_GO$result %>% 
  arrange(p_value)

###### Enrichment visualization via rrvgo ----------------------------------------

##
# sim Matrix
##
simMatrix_Mechs_UPF1_combined_Syn_GO <- rrvgo::calculateSimMatrix(gostres_Mechs_UPF1_combined_Syn_GO_result %>% 
                                                          distinct(term_id) %>% 
                                                          pull(term_id),
                                                        orgdb="org.Hs.eg.db",
                                                        ont="BP",
                                                        method="Rel")

##
# scores
##

scores_Mechs_UPF1_combined_Syn_GO <- setNames(-log10(gostres_Mechs_UPF1_combined_Syn_GO_result %>% 
                                                    distinct(term_id, .keep_all = TRUE) %>% 
                                                    pull(p_value)), 
                                              gostres_Mechs_UPF1_combined_Syn_GO_result %>% 
                                             distinct(term_id) %>% 
                                             pull(term_id))

##
# reduced terms
##

reducedTerms_Mechs_UPF1_combined_Syn_GO <-  rrvgo::reduceSimMatrix(simMatrix_Mechs_UPF1_combined_Syn_GO,
                                                        scores_Mechs_UPF1_combined_Syn_GO,
                                                        threshold=0.7,
                                                        orgdb="org.Hs.eg.db")

# Combine G:profiler output with combined&reduced GO IDs
gostres_Mechs_UPF1_combined_Syn_GO_result_reduced <- gostres_Mechs_UPF1_combined_Syn_GO_result %>% 
  left_join(reducedTerms_Mechs_UPF1_combined_Syn_GO %>% 
              dplyr::select(go, parent, parentTerm) %>% 
              dplyr::rename("term_id" = "go",
                            "parent_id" = "parent")) %>% 
  filter(!is.na(parent_id)) %>% 
  group_by(parentTerm,query) %>% 
  mutate(group_score = mean(-log10(p_value))) %>% 
  ungroup()

# Save as csv
gostres_Mechs_UPF1_combined_Syn_GO_result_reduced %>% write_csv("Resources/bakR/gostres_Mechs_UPF1_combined_Syn_GO_result_reduced.csv")