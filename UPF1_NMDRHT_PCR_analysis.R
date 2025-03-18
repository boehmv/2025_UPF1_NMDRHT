#!/usr/bin/env Rscript

# Title: UPF1_NMDRHT_PCR_analysis
# Objective: Code for generating PCR analysis plots for "Rapid UPF1 depletion illuminates the temporal dynamics of the NMD-regulated human transcriptome"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(readxl)
library(janitor)
library(tidyverse)
library(ggh4x)

# Fig1 KD SRSF2 -----------------------------------------------------------

Fig1_KD_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                       sheet = "Fig1_KD_SRSF2") %>% 
  clean_names() %>% 
  mutate(knockdown_time = paste0(knockdown,"_",time)) %>% 
  mutate(knockdown_time = fct_inorder(knockdown_time))

Fig1_KD_SRSF2 %>% 
  ggplot(aes(x=knockdown_time,
             y=percent,
             fill=knockdown)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
               show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("Luc" = "#B7B7BA",
                             "SMG6+7" = "#CBAC7C",
                             "UPF1" = "#66B2B2",
                             "siPOOL" = "#006666")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_KD_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 KD GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_KD_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                            sheet = "Fig1_KD_qPCR") %>% 
  clean_names() %>% 
  mutate(knockdown_time = paste0(knockdown,"_",time)) %>% 
  mutate(knockdown_time = fct_inorder(knockdown_time))

Fig1_KD_qPCR %>% 
  ggplot(aes(x=knockdown_time,
             y=log2fc,
             fill=knockdown)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("Luc" = "#B7B7BA",
                             "SMG6+7" = "#CBAC7C",
                             "UPF1" = "#66B2B2",
                             "siPOOL" = "#006666")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_KD_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 AID SRSF2 -----------------------------------------------------------

Fig1_AID_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                            sheet = "Fig1_AID_SRSF2") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_SRSF2 %>% 
  ggplot(aes(x=IAA_time,
             y=percent,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0" = "#F9F9F9",
                             "FALSE_48" = "#636363",
                             "TRUE_0" = "#65BB8C",
                             "TRUE_2" = "#43A787",
                             "TRUE_4" = "#229380",
                             "TRUE_8" = "#007E78",
                             "TRUE_12" = "#00696D",
                             "TRUE_24" = "#005560",
                             "TRUE_48" = "#0B4151"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(45, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 AID GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_AID_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                           sheet = "Fig1_AID_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_qPCR %>% 
  ggplot(aes(x=IAA_time,
             y=log2fc,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0" = "#F9F9F9",
                             "FALSE_48" = "#636363",
                             "TRUE_0" = "#65BB8C",
                             "TRUE_2" = "#43A787",
                             "TRUE_4" = "#229380",
                             "TRUE_8" = "#007E78",
                             "TRUE_12" = "#00696D",
                             "TRUE_24" = "#005560",
                             "TRUE_48" = "#0B4151"
                             )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(45, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 AID recovery SRSF2 -----------------------------------------------------------

Fig1_AID_recovery_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                             sheet = "Fig1_AID_recovery_SRSF2") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time,"_",recovery)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_recovery_SRSF2 %>% 
  ggplot(aes(x=IAA_time,
             y=percent,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("TRUE_0_0" = "#65BB8C",
                             "TRUE_24_0" = "#682714",
                             "TRUE_24_2" = "#953C00",
                             "TRUE_24_4" = "#BD5500",
                             "TRUE_24_8" = "#E17100",
                             "TRUE_24_12" = "#EF9300",
                             "TRUE_24_24" = "#F5B300",
                             "TRUE_24_48" = "#F8CF69"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_recovery_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 AID recovery GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_AID_recovery_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                            sheet = "Fig1_AID_recovery_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time,"_",recovery)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_recovery_qPCR %>% 
  ggplot(aes(x=IAA_time,
             y=log2fc,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("TRUE_0_0" = "#65BB8C",
                             "TRUE_24_0" = "#682714",
                             "TRUE_24_2" = "#953C00",
                             "TRUE_24_4" = "#BD5500",
                             "TRUE_24_8" = "#E17100",
                             "TRUE_24_12" = "#EF9300",
                             "TRUE_24_24" = "#F5B300",
                             "TRUE_24_48" = "#F8CF69"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_recovery_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 FKBP SRSF2 -----------------------------------------------------------

Fig1_FKBP_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                             sheet = "Fig1_FKBP_SRSF2") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(d_tag,"_",time,"_",cell_line)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_FKBP_SRSF2 %>% 
  ggplot(aes(x=IAA_time,
             y=percent,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_HCT116" = "#F9F9F9",
                             "TRUE_0_HCT116" = "#C8C3E2",
                             "TRUE_24_HCT116" = "#6B56A7",
                             "FALSE_0_HEK293" = "#F9F9F9",
                             "TRUE_0_HEK293" = "#DD8492",
                             "TRUE_24_HEK293" = "#963B5A"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  # facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_FKBP_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig1 FKBP GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_FKBP_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                            sheet = "Fig1_FKBP_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(d_tag,"_",time,"_",cell_line)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_FKBP_qPCR %>% 
  ggplot(aes(x=IAA_time,
             y=log2fc,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_HCT116" = "#F9F9F9",
                             "TRUE_0_HCT116" = "#C8C3E2",
                             "TRUE_24_HCT116" = "#6B56A7",
                             "FALSE_0_HEK293" = "#F9F9F9",
                             "TRUE_0_HEK293" = "#DD8492",
                             "TRUE_24_HEK293" = "#963B5A"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_FKBP_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig3 HU-6h SRSF2 -----------------------------------------------------------

Fig3_HMD_HU6_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                              sheet = "Fig3_HMD_HU6_SRSF2") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HU6_SRSF2 %>% 
  ggplot(aes(x=IAA_time,
             y=percent,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_20" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "FALSE_0_60" = "#65BB8C",
                             "TRUE_6_0" = "#11897C",
                             "TRUE_6_20" = "#11897C",
                             "TRUE_6_40" = "#11897C",
                             "TRUE_6_60" = "#11897C"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits = c(0,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  # facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HU6_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig3 HU-6h H1-2 & H2AC18 -----------------------------------------------------------

Fig3_HMD_HU6_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                             sheet = "Fig3_HMD_HU6_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HU6_qPCR %>% 
  ggplot(aes(x=IAA_time,
             y=log2fc,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_20" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "FALSE_0_60" = "#65BB8C",
                             "TRUE_6_0" = "#11897C",
                             "TRUE_6_20" = "#11897C",
                             "TRUE_6_40" = "#11897C",
                             "TRUE_6_60" = "#11897C"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HU6_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig3 HU-12-48h SRSF2 -----------------------------------------------------------

Fig3_HMD_HUlong_SRSF2 <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                                 sheet = "Fig3_HMD_HUlong_SRSF2") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HUlong_SRSF2 %>% 
  ggplot(aes(x=IAA_time,
             y=percent,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "TRUE_12_0" = "#00696D",
                             "TRUE_12_40" = "#00696D",
                             "TRUE_48_0" = "#0B4151",
                             "TRUE_48_40" = "#0B4151"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits = c(0,100)) +
  labs(y="Rel. levels of\nNMD isoform [%]",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  # facet_wrap(~cell_line) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HUlong_SRSF2.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Fig3 HU-12 & 48h H1-2 & H2AC18 -----------------------------------------------------------

Fig3_HMD_HUlong_qPCR <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_analysis.xlsx", 
                                sheet = "Fig3_HMD_HUlong_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HUlong_qPCR %>% 
  ggplot(aes(x=IAA_time,
             y=log2fc,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "TRUE_12_0" = "#00696D",
                             "TRUE_12_40" = "#00696D",
                             "TRUE_48_0" = "#0B4151",
                             "TRUE_48_40" = "#0B4151"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="log2FC",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HUlong_qPCR.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

##
# Raw cT values -----------------------------------------------------------
##

## Fig1 KD GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_KD_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                           sheet = "Fig1_KD_qPCR") %>% 
  clean_names() %>% 
  mutate(knockdown_time = paste0(knockdown,"_",time)) %>% 
  mutate(knockdown_time = fct_inorder(knockdown_time))

Fig1_KD_qPCR_raw %>% 
  ggplot(aes(x=knockdown_time,
             y=c_t_value,
             fill=knockdown)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("Luc" = "#B7B7BA",
                             "SMG6+7" = "#CBAC7C",
                             "UPF1" = "#66B2B2",
                             "siPOOL" = "#006666")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="",
       fill="",
       title="Figure 1 | RT-PCR | SRSF2") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_KD_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Fig1 AID GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_AID_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                            sheet = "Fig1_AID_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_qPCR_raw %>% 
  ggplot(aes(x=IAA_time,
             y=c_t_value,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0" = "#F9F9F9",
                             "FALSE_48" = "#636363",
                             "TRUE_0" = "#65BB8C",
                             "TRUE_2" = "#43A787",
                             "TRUE_4" = "#229380",
                             "TRUE_8" = "#007E78",
                             "TRUE_12" = "#00696D",
                             "TRUE_24" = "#005560",
                             "TRUE_48" = "#0B4151"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(45, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Fig1 AID recovery GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_AID_recovery_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                                     sheet = "Fig1_AID_recovery_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa,"_",time,"_",recovery)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_AID_recovery_qPCR_raw %>% 
  ggplot(aes(x=IAA_time,
             y=c_t_value,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("TRUE_0_0" = "#65BB8C",
                             "TRUE_24_0" = "#682714",
                             "TRUE_24_2" = "#953C00",
                             "TRUE_24_4" = "#BD5500",
                             "TRUE_24_8" = "#E17100",
                             "TRUE_24_12" = "#EF9300",
                             "TRUE_24_24" = "#F5B300",
                             "TRUE_24_48" = "#F8CF69"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1",
       subtitle="Outlier TRUE_24_2 replicate #3 -> removed from further analysis") +
  facet_wrap(~cell_line+target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_AID_recovery_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Fig1 FKBP GAS5 & ZFAS1 -----------------------------------------------------------

Fig1_FKBP_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                             sheet = "Fig1_FKBP_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(d_tag,"_",time,"_",cell_line)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig1_FKBP_qPCR_raw %>% 
  ggplot(aes(x=IAA_time,
             y=c_t_value,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_HCT116" = "#F9F9F9",
                             "TRUE_0_HCT116" = "#C8C3E2",
                             "TRUE_24_HCT116" = "#6B56A7",
                             "FALSE_0_HEK293" = "#F9F9F9",
                             "TRUE_0_HEK293" = "#DD8492",
                             "TRUE_24_HEK293" = "#963B5A"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig1_FKBP_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Fig3 HU-6h H1-2 & H2AC18 -----------------------------------------------------------

Fig3_HMD_HU6_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                                sheet = "Fig3_HMD_HU6_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HU6_qPCR_raw %>% 
  ggplot(aes(x=IAA_time,
             y=c_t_value,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_20" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "FALSE_0_60" = "#65BB8C",
                             "TRUE_6_0" = "#11897C",
                             "TRUE_6_20" = "#11897C",
                             "TRUE_6_40" = "#11897C",
                             "TRUE_6_60" = "#11897C"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(40, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HU6_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Fig3 HU-12 & 48h H1-2 & H2AC18 -----------------------------------------------------------

Fig3_HMD_HUlong_qPCR_raw <- readxl::read_excel("Resources/PCR/UPF1_NMDRHT_PCR_rawData.xlsx", 
                                   sheet = "Fig3_HMD_HUlong_qPCR") %>% 
  clean_names() %>% 
  mutate(IAA_time = paste0(iaa_time,"_",time)) %>% 
  mutate(IAA_time = fct_inorder(IAA_time))

Fig3_HMD_HUlong_qPCR_raw %>% 
  ggplot(aes(x=IAA_time,
             y=c_t_value,
             fill=IAA_time)) +
  stat_summary(fun=mean, geom="bar", show.legend=FALSE,
               # fill="black",
               width=0.25,
               alpha=0.75)  +
  stat_summary(fun=mean, color="black",geom="crossbar", show.legend=FALSE,
               width=0.8,
               linewidth=0.1329,
               lineend = "round")  +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              stackratio=1.5, dotsize=2.2,
  #              stroke=0.2) +
  geom_point(shape=21,
             size=1.1965,
             show.legend=FALSE,
             position = position_dodge2(width = 0.8, preserve = "total"),
             stroke=0.2,
             color="black") +
  theme_bw(base_rect_size = 0.1) +
  theme(strip.background = element_rect(fill="white",
                                        linewidth=0.2), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5,
                                    color="gray20"),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_fill_manual(values=c("FALSE_0_0" = "#65BB8C",
                             "FALSE_0_40" = "#65BB8C",
                             "TRUE_12_0" = "#00696D",
                             "TRUE_12_40" = "#00696D",
                             "TRUE_48_0" = "#0B4151",
                             "TRUE_48_40" = "#0B4151"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_y_continuous(breaks=c(0,25,50,75,100)) +
  labs(y="cT value",
       x="IAA_time",
       fill="",
       title="Figure 1 | qPCR | GAS5 & ZFAS1") +
  facet_wrap(~target) +
  force_panelsizes(cols = unit(30, "mm"),
                   rows = unit(13, "mm"))

ggsave(filename = "Plots/PCR/Fig3_HMD_HUlong_qPCR_raw_cT.pdf", 
       width = 15,
       height = 10,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")
