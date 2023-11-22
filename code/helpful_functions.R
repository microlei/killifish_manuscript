# Plotting helpers
library(tidyverse)
library(magrittr)

library(phyloseq)
library(speedyseq)
library(vegan)

library(cowplot)
library(patchwork)
library(RColorBrewer)


factor.fishType <- c("New Bedford Harbor wild", "New Bedford Harbor F2", "Scorton Creek wild", "Scorton Creek F2", "New Bedford Harbor water", "Scorton Creek water")

cols.fishType <- c("New Bedford Harbor wild"="#779498",
                   "New Bedford Harbor F2" = "#849e67",
                   "Scorton Creek wild"= "#914e76",
                   "Scorton Creek F2"= "#c25f3b",
                   "New Bedford Harbor water" = "#8c58c6",
                   "Scorton Creek water" = "#a5c846")
cols.wild_or_F2 <- c("wild"="brown1",
                     "F2"="cadetblue")
cols.site <- c("New Bedford Harbor"="coral3",
               "Scorton Creek"="turquoise4")
labels.fishType <- c("New Bedford Harbor wild"="Tolerant wild",
                     "New Bedford Harbor F2" = "Tolerant F2",
                     "Scorton Creek wild" = "Sensitive wild",
                     "Scorton Creek F2" = "Sensitive F2")
labels.short <- c("New Bedford Harbor wild"="NBH wild",
                  "New Bedford Harbor F2" = "NBH F2",
                  "Scorton Creek wild" = "SC wild",
                  "Scorton Creek F2" = "SC F2")
theme_set(theme_light()+ theme(legend.text = element_text(size=8),
                               plot.title=element_text(size=9, hjust=0.5, face = "bold")))

# function for cleaning corncob DA differentialTest output
# retrieves only the mu (abundance) and phi (variances) scores and not the intercepts for each significant model

cleanDA<- function(da, varname){
  # make empty tibble
  df <- tibble()
  # count the number of coefficients
  n <- da$significant_models[[1]]$coefficients %>% nrow()
  # select all the mu and phi values, dropping the intercept rows
  for(i in 1:length(da$significant_models)){
    df <- rbind(df, da$significant_models[[i]]$coefficients[seq(2,n,2),] %>% as.data.frame() %>% rownames_to_column(var="variable"))
  }
  # remove the variable name from the covariate column
  df$variable <- da$significant_models[[1]]$coefficients[c(2,4),] %>% as.data.frame() %>% rownames_to_column(var="variable") %>% use_series("variable") %>% gsub(x=., pattern=varname, replace="")
  # get the taxonomy of the significant taxa, repeated twice since (once for mu and once for phi)
  tax <- da$data %>% tax_table() %>% as("matrix") %>% data.frame() %>% rownames_to_column(var="ASV") %>%
    filter(ASV %in% da$significant_taxa) %>% dplyr::slice(rep(1:n(), each = n/2))
  # bind taxonomy columns to corncob measurements. Split the "variable" column into measurement (mu or phi) and the covariate name
  df <- bind_cols(df, tax) %>% as_tibble() %>% dplyr::rename(StdE = "Std. Error", t= "t value", p = "Pr(>|t|)") %>% 
    separate(col=variable, c("measure", varname), sep="[.]")
  return(df)
}

# small shortcut function that I add whenever I subset samples. Removes taxa with zero abundance
prune_zeros <- function(ps){
  prune_taxa(taxa_sums(ps)>0, ps)
}

# function makes a null matrix from an otu matrix by reshuffling observations for each ASV.
# first is the otu table and then second parameter is the number of matrices to make
make_null <- function(otu_table,x){
  randomize <- function(otu_table){
    nm <- sapply(1:nrow(otu_table), 
                 function(x){sample(as.numeric(otu_table[x,]))}) %>% t()
    rownames(nm) <- rownames(otu_table)
    colnames(nm) <- colnames(otu_table)
    return(nm)
  }
  replicate(x, randomize(otu_table))
}

# borrowed these functions from https://rdrr.io/github/mworkentine/mattsUtils/src/R/microbiome_helpers.R
#' Add taxonomy label
#'
#' add a column to the taxonomy table of a phyloseq object that lists the
#' lowest rank taxonomy assigned to that OTU along with a prefix indicating
#' the taxonomic rank.
#'
#' Example: g:Pseudomonas
#'
#' @param physeq a valid phyloseq object that contains a taxonomy table
#' @param num_species the number of species to retain if more than one are identified
#' @return a phyloseq object with an additional column on the taxonomy
#'         table called "Taxonomy"
#' @export

split_species = function(string, n = 2) {
  splits = str_split(string, "/", n + 1)
  res = map_if(splits, ~length(.x) > 2, ~.x[1:n]) %>%
    map_chr(str_c, collapse = "/")
  return(res)
}

add_taxonomy_column = function(physeq, num_species = 2) {
  tax_df = as.data.frame(tax_table(physeq)@.Data) %>%
    rownames_to_column("OTU") %>%
    mutate(Species = split_species(Species, n = num_species)) %>%
    mutate(Taxonomy =
             case_when(
               is.na(Class)  ~ str_c(OTU, "_p:", Phylum),
               is.na(Order)  ~ str_c(OTU, "_c:", Class),
               is.na(Family)  ~ str_c(OTU, "_o:", Order),
               is.na(Genus)   ~ str_c(OTU, "_f:", Family),
               is.na(Species) ~ str_c(OTU, "_g:", Genus),
               TRUE ~ str_c(OTU, "_", Genus, " ", Species)
             )
    )
  
  tax = as.matrix(tax_df[, -1])
  rownames(tax) = tax_df$OTU
  tax_table(physeq) = tax_table(tax)
  
  return(physeq)
}
