# code for making networks of the F2 taxa
# no pruning of less abundant taxa

library(tidyverse)
library(phyloseq)
library(SpiecEasi)

set.seed(1)
ps <- readRDS("../data/processed/ps.rds")

ps_gut_abund_NBH <- ps %>% subset_samples(fishType=="New Bedford Harbor F2") %>% prune_taxa(taxa_sums(.)>0,.)

hpc_params <- list(rep.num = 50, seed = 99999, ncores = 10)

print("Running spieceasi for NBH wild\n")
res_NBH <- spiec.easi(ps_gut_abund_NBH, verbose=TRUE, sel.criterion="bstars",
                  method = 'glasso',
                  lambda.min.ratio=1e-2,
                  nlambda=20,
                  pulsar.params=hpc_params)
print("Running spieceasi for SC wild\n")

ps_gut_abund_SC <- ps %>% subset_samples(fishType=="Scorton Creek F2") %>% prune_taxa(taxa_sums(.)>0,.)

res_SC <- spiec.easi(ps_gut_abund_SC, verbose=TRUE, sel.criterion = "bstars",
		method = 'glasso',
		lambda.min.ratio=1e-2,
		nlambda=50,
		pulsar.params=hpc_params)

saveRDS(res_NBH, file="../data/processed/spieceasi_F2_NBH.rds")
saveRDS(res_SC, file="../data/processed/spieceasi_F2_SC.rds")
