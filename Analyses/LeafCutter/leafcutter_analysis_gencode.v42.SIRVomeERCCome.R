#!/usr/bin/env Rscript

# Title: leafcutter_analysis_gencode.v42.SIRVomeERCCome.R
# Objective: Prepare LeafCutter output for further analysis. 
#            Standard leafcutter pipeline giving alternative splicing (AS) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(extrafont))

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# Define current analysis
myexp=arguments$args[2]

# Load cluster significance and split cluster information properly
cluster_significance <- read_delim(file.path(mydir, "leafcutter", myexp, paste0(myexp,"_cluster_significance.txt")),
                                   delim = "\t", 
                                   escape_double = FALSE, 
                                   trim_ws = TRUE) %>% 
  separate(col = cluster,
           into = c(NA, "clusterID"), 
           sep = ":")

# Load effect sizes and split intron information properly
effect_size <- read_delim(file.path(mydir, "leafcutter", myexp, paste0(myexp,"_effect_sizes.txt")),
                                   delim = "\t", 
                                   escape_double = FALSE, 
                                   trim_ws = TRUE) %>% 
  separate(col = intron,
           into = c("chr", "start", "end", "clusterID"), 
           sep = ":")

# Merge both tables and separate gene id/names properly
merged <- left_join(effect_size, cluster_significance) %>% 
  separate(col = genes,
           into = c("gene_name", "gene_id"), 
           sep = " ",
           extra="merge") %>% 
  separate(col = gene_id,
           into = c("gene_id", "altgenes"), 
           sep = ",",
           extra="merge") %>% 
  separate(col = clusterID,
           into = c(NA, NA, "strand"), 
           sep = "_",
           remove = FALSE)

# Write results file
write.csv(merged, 
          file=(file.path(mydir, 
                          "leafcutter", 
                          myexp, 
                          paste0(myexp,"_final.csv"))))


# Create "Plots" folder if it does not exists
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)

# Create "leafcutter" folder if it does not exists
dir.create(file.path(mydir, "Plots", "leafcutter"), showWarnings = FALSE)

# Create "leafcutter" folder if it does not exists
dir.create(file.path(mydir, "Plots", "leafcutter", paste(myexp)), showWarnings = FALSE)

# Generate first volcano plot without filtering
ggplot(data = merged, aes(x=deltapsi, y=-log10(p.adjust))) + 
  theme_classic() + 
  geom_point(alpha = 0.5, aes(color=deltapsi)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 12), 
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12), 
        plot.caption = element_text(size = 10), 
        text=element_text(family="Arial")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0.1, linetype="dashed") +
  geom_vline(xintercept = -0.1, linetype="dashed") +
  xlim(-1, 1) +
  scale_color_gradientn(colors = c("#980043","#ca0020", "#f4a582", "lightgray", "#92c5de", "#0571b0", "#253494"), values = rescale(c(-1,-0.5, -0.2, 0, 0.2, 0.5, 1), to = c(0, 1)), name = "deltaPSI") + labs(title = paste0(paste(myexp)), subtitle = "LeafCutter deltaPSI analysis", x = "deltaPSI", y = "-log10 p.adjust", caption = "LeafCutter (0.2.9)")

# Save plot
ggsave(file.path(mydir, "Plots", "leafcutter", paste(myexp), paste0(paste(myexp), "_raw_Volcano.pdf")), width = 6, height=4, 
       device=cairo_pdf,
       bg = "transparent")


####
# Filtered_data
####
my_data_filt <- merged %>% 
  filter(abs(deltapsi) > 0.2 & p.adjust < 0.01) %>%
  mutate(unique = paste(gene_id,clusterID, sep="_")) %>% 
  mutate(score = ".") %>% 
  filter(!is.na(strand))

# Calculate strand-aware 5'SS start and end position (3 bases in exon, 6 bases in intron)
my_data_filt$fiveSS_start <- ifelse(my_data_filt$strand=="+", as.double(my_data_filt$start)-3, as.double(my_data_filt$end)-7)

my_data_filt$fiveSS_end <- ifelse(my_data_filt$strand=="+", as.double(my_data_filt$start)+6, as.double(my_data_filt$end)+2)

# Calculate strand-aware 3'SS start and end position (3 bases in exon, 6 bases in intron)
my_data_filt$threeSS_start <- ifelse(my_data_filt$strand=="+", as.double(my_data_filt$end)-21, as.double(my_data_filt$start)-3)

my_data_filt$threeSS_end <- ifelse(my_data_filt$strand=="+", as.double(my_data_filt$end)+2, as.double(my_data_filt$start)+20)

# Generate final 5SS table to be output as bed file
output.5SS.bed <- subset(my_data_filt, select = c("chr", "fiveSS_start", "fiveSS_end","unique","score", "strand"))

# Generate final 3SS table to be output as bed file
output.3SS.bed <- subset(my_data_filt, select = c("chr", "threeSS_start", "threeSS_end","unique","score", "strand"))

# Export bed files
write.table(output.5SS.bed, file=file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), ".5SS.bed")), sep="\t", quote = F, col.names = F, row.names = F)

write.table(output.3SS.bed, file=file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), ".3SS.bed")), sep="\t", quote = F, col.names = F, row.names = F)

# Invoke bedtools from system - get fasta file corresponding to bed files
system(paste0("bedtools getfasta -fi /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa -bed ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".5SS.bed"), " -fo ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".5SS.fa.out"), " -s -name"))

system(paste0("bedtools getfasta -fi /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa -bed ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".3SS.bed"), " -fo ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".3SS.fa.out"), " -s -name"))

# Run MaxEnt score5.pl script
setwd("/home/volker/Tools/MaxEnt")

system(paste0("perl score5.pl ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".5SS.fa.out"), " > " , mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".5SS.MaxEnt.out")))

system(paste0("perl score3.pl ", mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".3SS.fa.out"), " > " , mydir, "/", "leafcutter", "/", paste(myexp), "/", paste0(paste(myexp), ".3SS.MaxEnt.out")))


# Read in MaxEnd outputs
MaxEnt_5SS_out <- read_tsv(file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), ".5SS.MaxEnt.out")), col_names = F)

MaxEnt_3SS_out <- read_tsv(file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), ".3SS.MaxEnt.out")), col_names = F)

# Retrieve 5'SS sequence and MaxEnt score -> incorporate in leafcutter table
my_data_filt$fiveSS_seq <- MaxEnt_5SS_out$X1

my_data_filt$fiveSS_MaxEnt <- MaxEnt_5SS_out$X2

# Retrieve 3'SS sequence and MaxEnt score -> incorporate in leafcutter table
my_data_filt$threeSS_seq <- MaxEnt_3SS_out$X1

my_data_filt$threeSS_MaxEnt <- MaxEnt_3SS_out$X2

# Categorize in "up"/"down"
my_data_filt$updown <- ifelse(my_data_filt$deltapsi>0.1, "up", "down")

# Plot MaxEnt 5SS output
ggplot(data = my_data_filt, aes(x=deltapsi, y=-log10(p.adjust))) + 
  theme_classic() + 
  geom_point(alpha = 0.5, aes(color=fiveSS_MaxEnt)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) +
  xlim(-1, 1) +
  scale_color_gradientn(colors = c("#980043","#ca0020", "#f4a582", "lightgray", "#92c5de", "#0571b0", "#253494"), values = rescale(c(-10,-5, -2, 0, 2, 5, 10), to = c(0, 1)), name = "5'ss MaxEnt score") + labs(title = paste0(paste(myexp)), subtitle = "LeafCutter 5'ss MaxEnt analysis", x = "deltaPSI", y = "-log10 p.adjust", caption = "LeafCutter (0.2.9) cut-offs: |deltaPSI| > 0.1 & p.adjust < 0.05")

# Save plot
ggsave(file.path(mydir, "Plots", "leafcutter", paste(myexp), paste0(paste(myexp), "_MaxEnt_5SS_Volcano.pdf")), width = 6, height=4, device=cairo_pdf)

# Plot MaxEnt 3SS output
ggplot(data = my_data_filt, aes(x=deltapsi, y=-log10(p.adjust))) + 
  theme_classic() + 
  geom_point(alpha = 0.5, aes(color=threeSS_MaxEnt)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) +
  xlim(-1, 1) +
  scale_color_gradientn(colors = c("#980043","#ca0020", "#f4a582", "lightgray", "#92c5de", "#0571b0", "#253494"), values = rescale(c(-10,-5, -2, 0, 2, 5, 10), to = c(0, 1)), name = "3'ss MaxEnt score") + labs(title = paste0(paste(myexp)), subtitle = "LeafCutter 3'ss MaxEnt analysis", x = "deltaPSI", y = "-log10 p.adjust", caption = "LeafCutter (0.2.9) cut-offs: |deltaPSI| > 0.1 & p.adjust < 0.05")

# Save plot
ggsave(file.path(mydir, "Plots", "leafcutter", paste(myexp), paste0(paste(myexp), "_MaxEnt_3SS_Volcano.pdf")), width = 6, height=4, device=cairo_pdf)

# Kolmogorov-Smirnov test for boxplot
up <- subset(my_data_filt, updown == "up")
down <- subset(my_data_filt, updown == "down")

ks_5SS_out <- ks.test(up$fiveSS_MaxEnt, down$fiveSS_MaxEnt)
ks_3SS_out <- ks.test(up$threeSS_MaxEnt, down$threeSS_MaxEnt)

# Plot MaxEnt 5SS output as boxplot
ggplot(data=my_data_filt, aes(x=updown, y=fiveSS_MaxEnt, fill=updown)) + geom_boxplot(notch=TRUE) + theme_classic() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) + scale_fill_manual(values = c("#ca0020", "#0571b0"), name = "5'ss MaxEnt score") + theme(legend.position="none") + labs(title = paste0(paste(myexp)), subtitle = "LeafCutter 5'ss MaxEnt analysis", x = "deltaPSI Up-/Down", y = "5'ss MaxEnt score", caption = "LeafCutter (0.2.9) cut-offs: |deltaPSI| > 0.1 & p.adjust < 0.05") + annotate("text", x = 1.5, y = 16, size = 3, label = paste0("Kolmogorov-Smirnov test p-value = ", paste(signif(ks_5SS_out$p.value, digits = 4))))

# Save plot
ggsave(file.path(mydir, "Plots", "leafcutter", paste(myexp), paste0(paste(myexp), "_MaxEnt_5SS_Boxplot.pdf")), width = 4, height=4, device=cairo_pdf)

# Plot MaxEnt 3SS output as boxplot
ggplot(data=my_data_filt, aes(x=updown, y=threeSS_MaxEnt, fill=updown)) + geom_boxplot(notch=TRUE) + theme_classic() + scale_fill_manual(values = c("#ca0020", "#0571b0"), name = "3'ss MaxEnt score") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10), text=element_text(family="Arial")) + theme(legend.position="none") + labs(title = paste0(paste(myexp)), subtitle = "LeafCutter 3'ss MaxEnt analysis", x = "deltaPSI Up-/Down", y = "3'ss MaxEnt score", caption = "LeafCutter (0.2.9) cut-offs: |deltaPSI| > 0.1 & p.adjust < 0.05") + annotate("text", x = 1.5, y = 16, size = 3, label = paste0("Kolmogorov-Smirnov test p-value = ", paste(signif(ks_3SS_out$p.value, digits = 4))))

# Save plot
ggsave(file.path(mydir, "Plots", "leafcutter", paste(myexp), paste0(paste(myexp), "_MaxEnt_3SS_Boxplot.pdf")), width = 4, height=4, device=cairo_pdf)

# Save Rdata
#save.image(file=file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), "_session.rData")))

# Export final filtered table
write.csv(my_data_filt, file=file.path(mydir, "leafcutter", paste(myexp), paste0(paste(myexp), "_filtered_table.csv")))

writeLines(capture.output(sessionInfo()), paste0(mydir, "/leafcutter/leafcutter_analysis_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file

