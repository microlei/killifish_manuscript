# Differential abundance results
library(corncob)
library(phyloseq)
library(tidyverse)
library(here)

source(here("code/helpful_functions.R"))

ps <- readRDS(here("data/processed/ps.rds"))

# differential test comparing wild vs F2 fish
# to save compute time, prun low abundance taxa
ps.abund <- ps %>% subset_samples(sampleType=="gut") %>% prune_taxa(taxa_sums(.)>20, .)
da_wild_f2 <- differentialTest(formula = ~wild_or_F2,
                               phi.formula = ~wild_or_F2,
                               formula_null = ~1,
                               phi.formula_null = ~wild_or_F2,
                               test="Wald", boot=FALSE,
                               data=ps.abund,
                               fdr_cutoff = 0.05)

# save the test for later
saveRDS(da_wild_f2, file=here("data/processed/da_wild_f2.rds"))
# reading in the .rds
da_wild_f2 <-readRDS(file=here("data/processed/da_wild_f2.rds"))
da_wild_f2_clean <- cleanDA(da_wild_f2, "wild_or_F2")%>% filter(measure=="mu")

# differential test comparing F2 fish
ps.f2.abund <- ps %>% subset_samples(wild_or_F2=="F2" & sampleType=="gut") %>% prune_taxa(taxa_sums(.)>20, .)
da_site_f2 <- differentialTest(formula = ~site,
                               phi.formula = ~site,
                               formula_null = ~1,
                               phi.formula_null = ~site,
                               test="Wald", boot=FALSE,
                               data=ps.f2.abund,
                               fdr_cutoff = 0.05)

saveRDS(da_site_f2, file=here("data/processed/da_site_f2.rds"))
da_site_f2 <- readRDS(here("data/processed/da_site_f2.rds"))
da_site_f2_clean <- cleanDA(da_site_f2, "site")%>% filter(measure=="mu")


# differential test comparing wild fish
ps.wild.abund <- ps %>% subset_samples(wild_or_F2=="wild" & sampleType=="gut") %>% prune_taxa(taxa_sums(.)>20, .)
da_site_wild <- differentialTest(formula = ~site,
                                 phi.formula = ~site,
                                 formula_null = ~1,
                                 phi.formula_null = ~site,
                                 test="Wald", boot=FALSE,
                                 data=ps.wild.abund,
                                 fdr_cutoff = 0.05)

saveRDS(da_site_wild, file=here("data/processed/da_site_wild.rds"))
da_site_wild <- readRDS(file=here("data/processed/da_site_wild.rds"))
da_site_wild_clean <- cleanDA(da_site_wild, "site") %>% filter(measure=="mu")


# Abundance of significantly different ASVs between wild fish in the water samples
ps.rel <- readRDS(here("data/processed/psrel.rds"))
sig_water <- ps.rel %>% prune_taxa(taxa_names(ps.rel) %in% da_site_wild_clean$ASV,.) %>% 
  subset_samples(sampleType=="water") %>% 
  psmelt() %>% 
  rename(ASV="OTU") %>% 
  left_join(da_site_wild_clean %>% filter(measure=="mu") %>% select(Estimate:ASV), by="ASV")

# Make a table of whether the ASV was enriched in both SC fish and SC water or in vice versa
sig_water <- sig_water %>% select(ASV, site, Estimate, Taxonomy, Abundance, Sample) %>% 
  mutate(SC_enriched = Estimate >0) %>% 
  group_by(ASV, site) %>% 
  summarise(avg=mean(Abundance), SC_enriched = SC_enriched, Taxonomy=Taxonomy) %>% 
  distinct() %>% 
  pivot_wider(id_cols=c(ASV,site, SC_enriched, Taxonomy), names_from=site, values_from=avg) %>% 
  mutate(expected = ifelse(SC_enriched, ifelse(`Scorton Creek` > `New Bedford Harbor`, "yes", "no"), ifelse(`New Bedford Harbor` > `Scorton Creek`, "yes", "no")))
sig_water <- rename(sig_water, NBH_abund="New Bedford Harbor", SC_abund = "Scorton Creek", SC_fish_enriched="SC_enriched")

idx <- taxa_names(ps) %in% sig_water$ASV
idx <- taxa_names(ps)[idx]
sig_water <- sig_water[match(idx, sig_water$ASV),]
write_csv(sig_water, file=here("figures/sig_water.csv"))
