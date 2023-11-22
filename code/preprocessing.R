# Script to clean data as described in the R notebook 1_Raw_read_stats.Rmd

# libraries
library(tidyverse)
library(readr)
library(phyloseq)
library(decontam)
library(here)

# helpful functions
source("code/helpful_functions.R")

# read in the files
taxonomy <- read.delim(here("output/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("output/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("data/KLF_metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
asvSeqs <- read.delim(here("output/ASVseqs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
tree <- read_tree(here("output/tree.nwk"))
track <- read.delim(here("output/track.tsv"), "\t", header=TRUE, row.names=1, check.names=FALSE)
track <- mutate(track, prop.retained = nochim/reads.in) 

seqtab <- readRDS(here("output/seqtab_nochimeras.rds"))
rownames(seqtab) <- rownames(track)

# make phyloseq
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))

# add the dada2 tracked reads
ps <- merge_phyloseq(ps, sample_data(track))

# get the proportion of chloroplast reads
ps.chloroplast <- ps %>% transform_sample_counts(function(x) x/sum(x)) %>% subset_taxa(Order=="Chloroplast") %>% psmelt() %>% select(OTU, sample, Abundance)
ps.chloroplast <- ps.chloroplast %>% group_by(sample) %>% summarise(chl_abundance=sum(Abundance)) %>% column_to_rownames("sample")
ps <- merge_phyloseq(ps, sample_data(ps.chloroplast))

# remove mock samples (mock samples were analyzed in the R notebook)
ps <- subset_samples(ps, sampleType!="mock")

# use decontam to remove contaminants
sample_data(ps)$is.neg <- grepl("control", sample_data(ps)$sampleType)
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
contaminant_asvs <- rownames(contamdf.prev[which(contamdf.prev$contaminant),])
ps <- prune_taxa(!taxa_names(ps) %in% contaminant_asvs, ps)

# remove controls
ps <- subset_samples(ps, sampleType %in% c("gut", "water")) %>% prune_taxa(taxa_sums(.)>0,.)

# remove chloroplasts and Eukaryotes
ps <- subset_taxa(ps, Kingdom!="Eukaryota") %>% subset_taxa(Order!="Chloroplast" | is.na(Order))
ps <- add_taxonomy_column(ps)

# transform to relative abundances
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))

# make melted phyloseq object
ps.melted <- psmelt(ps)
ps.rel.melted <- psmelt(ps.rel)
ps.melted <- ps.melted %>% rename(abund_count = "Abundance") %>% 
  left_join(ps.rel.melted %>% rename(abund_rel="Abundance") %>% select(1:3), 
            by=c("OTU", "Sample"))

# export cleaned files
new.asvs <- otu_table(ps)
new.taxonomy <- tax_table(ps)
new.metadata <- as(sample_data(ps), "data.frame")
new.asvSeqs <- asvSeqs[which(asvSeqs$asv %in% taxa_names(ps)),]

# write the cleaned objects to file
saveRDS(ps.melted, file=here("data/processed/psmelted.rds"))
saveRDS(ps, file=here("data/processed/ps.rds"))
saveRDS(ps.rel, file=here("data/processed/psrel.rds"))
write.table(new.asvs, file=here("data/processed/ASVs.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.taxonomy, file=here("data/processed/taxonomy.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.metadata, file=here("data/processed/metadata.csv"), sep=",", col.names = NA, quote=FALSE)
write.table(new.asvSeqs, file=here("data/processed/asvSeqs.txt"), sep="\t", col.names = NA, quote = FALSE)
