# code for processing spieceasi results into the covariance matrix, the inverse covariance matrix, and the network fit
library(tidyverse)
library(SpiecEasi)
library(here)

se_NBH <- readRDS(here("data/processed/spieceasi_wild_NBH.rds"))
se_SC <- readRDS(here("data/processed/spieceasi_wild_SC.rds"))
se_NBH_F2 <- readRDS(here("data/processed/spieceasi_F2_NBH.rds"))
se_SC_F2 <- readRDS(here("data/processed/spieceasi_F2_SC.rds"))

se_NBH_cov <- getOptCov(se_NBH)
se_SC_cov <- getOptCov(se_SC)
se_NBH_F2_cov <- getOptCov(se_NBH_F2)
se_SC_F2_cov <- getOptCov(se_SC_F2)

se_NBH_icov <- getOptiCov(se_NBH)
se_SC_icov <- getOptiCov(se_SC)
se_NBH_F2_icov <- getOptiCov(se_NBH_F2)
se_SC_F2_icov <- getOptiCov(se_SC_F2)

se_NBH_refit <- getRefit(se_NBH)
se_SC_refit <- getRefit(se_SC)
se_NBH_F2_refit <- getRefit(se_NBH_F2)
se_SC_F2_refit <- getRefit(se_SC_F2)

save(se_NBH_cov, se_SC_cov, se_NBH_F2_cov, se_SC_F2_cov,
  se_NBH_icov, se_SC_icov, se_NBH_F2_icov, se_SC_F2_icov,
  se_NBH_refit, se_SC_refit, se_NBH_F2_refit, se_NBH_F2_refit,
  file=here("data/processed/se_results.RData"))
