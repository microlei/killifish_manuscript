# Code for plotting the venn diagram and bar plots

library(tidyverse)
library(phyloseq)
library(speedyseq)
library(VennDiagram)
library(corncob)
library(vegan)
library(cowplot)
library(patchwork)
library(magrittr)
library(here)

source(here("code/helpful_functions.R"))
source(here("code/plot_venn.R"))
source(here("code/plot_bar_f2_wild.R"))
source(here("code/plot_vibrio_pie.R"))

p_venn <- ggdraw() + draw_image(here("figures/Venn_gut_color.png"))
p_bar_full <- readRDS(here("data/processed/p_bar_full.rds"))
p_vibrio <- ggdraw() + draw_image(here("figures/plot_vibrio_pie.png"))

p_combined <- plot_grid(p_venn, p_vibrio, labels=c("A", "B")) %>% plot_grid(p_bar_full, nrow=2, labels=c("", "C"))
ggsave(p_combined, filename = here("figures/figure_2.pdf"), width=6, height=5, unit="in")
