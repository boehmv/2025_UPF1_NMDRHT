#!/usr/bin/env Rscript

# Combined analysis of Leafcutter output

# Created by: boehmv (Volker Böhm)

######
# Load libraries
######

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("readxl"))
suppressMessages(library(extrafont))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
suppressPackageStartupMessages(library(ComplexHeatmap))
library("RColorBrewer")
suppressPackageStartupMessages(library(circlize))
library(extrafont)
#suppressPackageStartupMessages(library(tidyHeatmap))
library(reshape2)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 1)

# mydir is the first argument
mydir=arguments$args[1]

# Leafcutter: read in ConditionTable
df <- read.table(file = file.path(mydir, "leafcutter", "ConditionTable.txt"), sep = '\t', header = TRUE)

# Access the LeafViz result files
df2 <- df %>% 
  mutate(basePath = dirname(path)) %>% 
  mutate(intronPath = paste0(basePath,"/control_vs_",cond,"_introns.csv")) %>% 
  mutate(clusterPath = paste0(basePath,"/control_vs_",cond,"_clusters.csv")) %>% 
  mutate(intronSummaryPath = paste0(basePath,"/control_vs_",cond,"_intron_summary.csv")) %>% 
  mutate(clusterSummaryPath = paste0(basePath,"/control_vs_",cond,"_cluster_summary.csv")) 

# Collect into cat files
introns_cat <- df2 %>% 
  pull(intronPath) %>% 
  map_dfr(read_csv, col_select = (-1)) 

clusters_cat <- df2 %>% 
  pull(clusterPath) %>% 
  map_dfr(read_csv, col_select = (-1)) 

intronSummary_cat <- df2 %>% 
  pull(intronSummaryPath) %>% 
  map_dfr(read_csv, col_select = (-1)) 

clusterSummary_cat <- df2 %>% 
  pull(clusterSummaryPath) %>% 
  map_dfr(read_csv, col_select = (-1)) 

# Output
introns_cat %>% 
  write_csv(file=file.path(mydir, "leafcutter", paste0("leafcutter_introns_cat",
                                                               ".csv")))

clusters_cat %>% 
  write_csv(file=file.path(mydir, "leafcutter", paste0("leafcutter_clusters_cat",
                                                       ".csv")))

intronSummary_cat %>% 
  write_csv(file=file.path(mydir, "leafcutter", paste0("leafcutter_intronSummary_cat",
                                                       ".csv")))

clusterSummary_cat %>% 
  write_csv(file=file.path(mydir, "leafcutter", paste0("leafcutter_clusterSummary_cat",
                                                       ".csv")))


# Generate empty data frame for concatenated file
AS_cat <- data.frame(matrix(ncol = 7, nrow = 0))
y <- c("gene_id","gene_name","chr","start","end","deltapsi","p.adjust")
colnames(AS_cat) <- y

# Iterate over conditions and generate data frames, generate combined df by merge
for (row in 1:nrow(df)) {
  cond=df[row, "cond"]
  path=df[row, "path"]
  
  # Read in Excel files, subset for relevant genes and columns, concatenate into AS_cat
  d <- read_csv(file.path(path))
  d_cat <- subset(d, !grepl("^SIRV", gene_id) & !grepl("^ERCC", gene_id), select = c("gene_id","gene_name","chr","start","end","deltapsi","p.adjust"))
  d_cat$condition_2 <- paste0(cond)
  AS_cat <- rbind(AS_cat, d_cat)
}

write.csv(AS_cat, file=file.path(mydir, "leafcutter", paste0("leafcutter_AS_cat",
                                                                 ".csv")))

# Produce Tibble, group by and summarise 
AS_Tbl <- as_tibble(AS_cat) %>% 
  filter(p.adjust < 0.01 & abs(deltapsi) > 0.2) %>% 
  mutate(gene_coordinates=paste0(gene_id,"-",chr,":",start,"-",end)) %>% 
  select(-c("gene_id", "gene_name", "p.adjust", "chr", "start", "end"))  
    
# Pivot wider in order to get matrix in the end
AS_Tbl_wide <- pivot_wider(AS_Tbl, names_from = condition_2, values_from = deltapsi)

# Remove NA (replace with 0)
AS_Tbl_wide[is.na(AS_Tbl_wide)] <- 0

combined_heat <- as.data.frame(AS_Tbl_wide)

combined_heat <- combined_heat %>% column_to_rownames(var = "gene_coordinates")

# Order data frame
HeatcolSums <- colSums(combined_heat != 0)

HeatcolSums <- HeatcolSums[order(-HeatcolSums)]

combined_heat <- combined_heat[,names(HeatcolSums)]

i <- length(HeatcolSums)

if (i>1) {
  while (i>0) {
    combined_heat <- combined_heat[order(-combined_heat[,i]),]
    i=i-1
  }
} else {
  combined_heat <- sort(combined_heat)
}


# Get Annotations of Up-/Downregulated genes
if (length(HeatcolSums)>1) {
  UpReg <- colSums(combined_heat > 0.1)
} else {
  UpReg <- sum(combined_heat > 0.1)
}

UpAnno = HeatmapAnnotation("dPSI Up" = anno_barplot(UpReg, gp = gpar(fill="#253494"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

if (length(HeatcolSums)>1) {
  DownReg <- colSums(combined_heat < -0.1)
} else {
  DownReg <- sum(combined_heat < -0.1)
}


DownAnno = HeatmapAnnotation("dPSI Down" = anno_barplot(DownReg, gp = gpar(fill="#980043"), axis_param = list(gp = gpar(fontsize=6)), height = unit(0.75, "cm")), annotation_name_gp = gpar(fontsize=6))

######
# Heatmap
######
# Set Color for Heatmap
col = colorRamp2(c(min(combined_heat),-0.5, -0.1, 0, 0.1, 0.5, max(combined_heat)), c("#980043","#ca0020", "#f4a582", "white", "#92c5de", "#0571b0", "#253494"))

# Generate Heatmap
heat <- Heatmap(as.matrix(combined_heat),
                name = "deltaPSI",
                heatmap_legend_param = list(title_position = "topleft", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)),
                col =  col,
                na_col = "white",
                cluster_rows = FALSE,
                cluster_columns = FALSE,	# Set to TRUE if clustering is desired
                column_dend_side = "top",
                border = TRUE,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                column_dend_reorder = FALSE,
                row_dend_reorder = FALSE,
                column_names_gp = gpar(fontsize=6),
                top_annotation = UpAnno,
                bottom_annotation = DownAnno,
                width = unit((length(HeatcolSums)/2), "cm"),
                height = unit(5, "cm")
)

# Print to pdf
pdf(file = file.path(mydir, "Plots", "leafcutter", "LeafCutter_all_heatmap_string.pdf"))

draw(heat)

dev.off()