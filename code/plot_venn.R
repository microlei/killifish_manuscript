# code for plotting venn diagram
library(tidyverse)
library(phyloseq)
library(speedyseq)
library(VennDiagram)
library(here)

source(here("code/helpful_functions.R"))

ps <- readRDS(here("data/processed/ps.rds"))

taxa_list <- list(nb_w = subset_samples(ps, fishType=="New Bedford Harbor wild") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  nb_f2 = subset_samples(ps, fishType=="New Bedford Harbor F2") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  sc_w = subset_samples(ps, fishType=="Scorton Creek wild") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  sc_f2 = subset_samples(ps, fishType=="Scorton Creek F2") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names())

taxa_overlap <- calculate.overlap(taxa_list)
f2 <- calculate.overlap(list(taxa_list$sc_f2, taxa_list$nb_f2))
wild <- calculate.overlap(list(taxa_list$sc_w, taxa_list$nb_w))

overlap <- list(all = length(taxa_overlap$a6), 
                f2 = length(f2$a3),
                f2_p = 100*length(f2$a3)/length(unique(c(f2$a1, f2$a2))),
                wild = length(wild$a3),
                wild_p = 100*length(wild$a3)/length(unique(c(wild$a1, wild$a2))))

saveRDS(overlap, here("data/processed/overlap.rds"))
asvSeqs <- read.delim(here("data/processed/asvSeqs.txt"), "\t", header=TRUE, row.names=1, check.names=FALSE)
overlap_csv<-tax_table(ps) %>% as.data.frame() %>% rownames_to_column(var="asv") %>% filter(asv %in% taxa_overlap$a6) %>% left_join(asvSeqs)
write_csv(overlap_csv, here("figures/shared_taxa.csv"))

# the below has already been done
venn.diagram(taxa_list, 
             filename = here("figures/Venn_gut_color.png"), disable.logging = TRUE, imagetype="png", 
             category.names = c("NBH wild", "SC F2", "SC wild", "SC F2"), 
             fill=cols.fishType[1:4], alpha=0.5, cat.col=cols.fishType[1:4], cat.cex=1.5, margin=.05, print.mode="percent")
