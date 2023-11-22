# generate permanova and dispersion tests results table
library(tidyverse)
library(phyloseq)
library(speedyseq)
library(microbiome)
library(vegan)
library(broom)
library(here)

source(here("code/helpful_functions.R"))

set.seed(100)

ps <- readRDS(here("data/processed/ps.rds"))
ps.rel <- readRDS(here("data/processed/psrel.rds"))
f <- zCompositions::cmultRepl(t(otu_table(ps)), method="CZM", label=0, output="p-counts") %>% t()
ps.zcomp <- ps
otu_table(ps.zcomp) <- otu_table(f, taxa_are_rows = TRUE)
# NOTE SUBSET ONLY GUT SAMPLES HERE
ps.zcomp <- transform(ps.zcomp, "clr") %>% subset_samples(sampleType=="gut") %>% prune_zeros()
dist.aitch <- distance(ps.zcomp, "euclidean")
dat <- ps.zcomp %>% sample_data() %>% as("data.frame")

# PERMANOVA
perm <- how(nperm=999)
perm_fishType <- adonis2(dist.aitch ~ dat$fishType, data = dat) # PERMANOVA by fish type
perm_wild_F2 <- adonis2(dist.aitch ~ dat$wild_or_F2, data=dat) # PERMANOVA by wild or F2
dist_wild <- ps.zcomp %>% subset_samples(wild_or_F2 == "wild") %>% prune_zeros() %>% distance("euclidean")
dat_wild <- dat %>% filter(wild_or_F2=="wild")
dist_F2 <- ps.zcomp %>% subset_samples(wild_or_F2 == "F2") %>% prune_zeros() %>% distance("euclidean")
dat_F2 <- dat %>% filter(wild_or_F2 == "F2")
perm_wild <- adonis2(dist_wild ~ dat_wild$fishType, data=dat_wild) # PERMANOVA SC WILD VS NBH WILD
perm_F2 <- adonis2(dist_F2 ~ dat_F2$fishType, data=dat_F2) # PERMANOVA SC F2 VS NBH F2
dist_env <- subset_samples(ps.zcomp, !is.na(weight.g) & !is.na(length.cm) & !is.na(sex)) %>% distance("euclidean") # just making sure that fish physical traits are included
perm_physiology <- adonis2(dist_env ~ dat$fishType + dat$weight.g + dat$length.cm + dat$sex, data = dat, by="margin")

# Dispersion 
disper_fishType <-  betadisper(d=dist.aitch, dat$fishType)
permutest(disper_fishType)
TukeyHSD(disper_fishType) # signif NBH wild vs all and SC wild vs NBH F2. not signif F2 fish, SC wild vs SC F2. 


# Tidy permanova into a table
perm_results <- bind_rows(tidy(perm_fishType)[1, 2:6], tidy(perm_wild_F2)[1,2:6], tidy(perm_wild)[1,2:6], tidy(perm_F2)[1, 2:6]) %>% 
  bind_cols(tibble("Comparison" = c("All fish types", "Wild vs F2", "SC wild vs NBH wild", "SC F2 vs NBH F2")),.)
colnames(perm_results) <- c("Comparison", "Df", "Sum Sq", "R2", "Pseudo-F", "P")

perm_supp <- tidy(perm_physiology)[2:6] %>% bind_cols(tibble("Term" = c("Fish type", "Weight", "Length", "Sex", "Residual", "Total")),.)
colnames(perm_supp) <- c("Term", "Df", "Sum Sq", "R2", "Pseudo-F", "P")
write_csv(perm_results, here("figures/perm_results.csv"))
write_csv(perm_supp, here("figures/perm_supp.csv"))
