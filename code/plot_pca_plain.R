# Code to plot the PCA of all samples and then subplots of wild vs captive
library(zCompositions)
library(tidyverse)
library(ggrepel)
library(speedyseq)
library(microbiome)
library(vegan)
library(here)

source(here("code/helpful_functions.R"))

set.seed(100)

# load data
ps <- readRDS(here("data/processed/ps.rds"))
metadata <- sample_data(ps) %>% as("data.frame")
taxonomy <- ps %>% tax_table %>% as_tibble %>% rename(ASV=".otu")

# Use zComposition to impute zeros before performing a center log ratio transform
f <- zCompositions::cmultRepl(t(otu_table(ps)), method="CZM", label=0, output="p-counts") %>% t()
ps.zcomp <- ps
otu_table(ps.zcomp) <- otu_table(f, taxa_are_rows = TRUE)
ps.zcomp <- transform(ps.zcomp, "clr")

# RDA with no other options is the same as CA/PCA
pca <- rda(t(otu_table(ps.zcomp)))
# Join the sample scores with the metadata
pca.df <- scores(pca, display="sites") %>% data.frame() %>% rownames_to_column(var="sample") %>% left_join(metadata %>% dplyr::select(sample, site, wild_or_F2, fishType))
# Join the species (ASV) loadings with the taxonomy
pca.species <- scores(pca, display="species") %>% data.frame() %>% rownames_to_column(var="ASV") %>% left_join(taxonomy)
# get summary
pca.summary <- summary(pca)

# normal PCA with no arrows
pca.plot <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=fishType), size=4, position="jitter") +
  scale_color_manual(name="", values=cols.fishType) +
  labs(x=str_c("PC1 [", round(pca.summary$cont$importance[2,1]*100,2),"%]"),
       y=str_c("PC2 [", round(pca.summary$cont$importance[2,2]*100,2),"%]"))

#ggsave(filename=here("figures/plot_pca_plain.png"), width = 220, height =180, unit="mm" )

# Alternate with ggforce() ellipses
library(ggforce)
pca.plot <- pca.plot +
  geom_mark_ellipse(aes(fill=wild_or_F2, label=wild_or_F2), alpha=0, label.fontsize=10) +
  guides(color = guide_legend(ncol=1), fill="none") +
  theme(legend.position = "right", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm"), legend.text=element_text(size=10)) +
  coord_cartesian()

wild_subplot_pca <- ps.zcomp %>% subset_samples(wild_or_F2=="wild") %>% otu_table() %>% t() %>% rda()
wild_summary <- summary(wild_subplot_pca)
wild_subplot <- scores(wild_subplot_pca, display="sites") %>% data.frame() %>% rownames_to_column(var="sample") %>% left_join(metadata %>% dplyr::select(sample, site, wild_or_F2, fishType)) %>% 
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=fishType), size=4, position="jitter") +
  scale_color_manual(name="", values=cols.fishType) + theme(legend.position = "none") + 
  labs(x=str_c("PC1 [", round(wild_summary$cont$importance[2,1]*100,2),"%]"),
       y=str_c("PC2 [", round(wild_summary$cont$importance[2,2]*100,2),"%]"),
       title="Wild fish")

F2_subplot_pca <- ps.zcomp %>% subset_samples(wild_or_F2=="F2") %>% otu_table() %>% t() %>% rda()
F2_summary <- summary(F2_subplot_pca)
F2_subplot <- scores(F2_subplot_pca, display="sites") %>% data.frame() %>% rownames_to_column(var="sample") %>% left_join(metadata %>% dplyr::select(sample, site, wild_or_F2, fishType)) %>% 
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=fishType), size=4, position="jitter") +
  scale_color_manual(name="", values=cols.fishType) + theme(legend.position = "none") + 
  labs(x=str_c("PC1 [", round(F2_summary$cont$importance[2,1]*100,2),"%]"),
       y=str_c("PC2 [", round(F2_summary$cont$importance[2,2]*100,2),"%]"),
       title = "F2 fish")
pca_full <- plot_grid(wild_subplot, F2_subplot, ncol=1) %>% plot_grid(pca.plot, rel_widths = c(.5,1))
# Manually saved/expoted 900/400 px image
# ggsave(filename=here("figures/plot_pca_plain_full.png"), width = 6.5, height=3.25, units = "in")

