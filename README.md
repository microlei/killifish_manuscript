Environmental and Population influences on Atlantic Killifish Gut
Microbiomes
================
Lei Ma
Last compiled on 22 November, 2023

## Abstract

The Atlantic killifish, Fundulus heteroclitus, is an abundant estuarine
fish broadly distributed along the eastern coast of the U.S. which has
repeatedly evolved tolerance to otherwise lethal levels of aromatic
hydrocarbon exposure. This tolerance is linked to reduced activation of
the aryl hydrocarbon receptor (AHR) signaling pathway. In other animals,
the AHR has been shown to influence the gastrointestinal-associated
microbial community, or gut microbiome, particularly when activated by
the model toxic pollutant 3,3’,4,4’,5-pentachlorobiphenyl (PCB-126) and
other dioxin-like compounds. In order to understand host population and
PCB-126 exposure effects on killifish gut microbiota, we sampled two
populations of wild fish, one from from a PCB-exposed environment (New
Bedford Harbor) and the other from a pristine location (Scorton Creek),
as well as laboratory reared F2 generation fish originating from each of
these populations. We examined the bacteria and archaea associated with
the gut of these fish using amplicon sequencing of small subunit
ribosomal RNA genes and found that fish living in the PCB polluted site
had high microbial alpha and beta diversity and an altered network
structure compared to fish from the pristine site. These differences
were not present in laboratory reared F2 fish which originated from the
PCB polluted environment, though differences were also present between
wild and lab-based F2 fish from the pristine site. Microbial
compositional differences existed between the wild and lab reared fish,
with the wild dominated by Vibrionaceae and the lab reared by
Enterococceae. These results suggest that killifish habitat and/or
environmental conditions has a stronger influence on the killifish gut
microbiome compared to population or hereditary-based influences.
Atlantic killifish are important eco-evolutionary model organism and
this work reveal their importance for exploring
host-environmental-microbiome dynamics.

## Directory description

### Data

Contains the metadata and data files that have been generated during
analysis phase. Ignores large files (like .rds files) used in some
scripts and notebooks due to github storage limits.

### Logs

Logs of the Snakemake dada2 data processing of the raw sequence reads

### Output

Output of the Snakemake dada2 data processing, including the ASV table
and some intermediate products. Does not include the quality profiles of
the samples.

### Scripts

R scripts used in the Snakemake workflow

### Figures

Contains tables and figures from the manuscript. Includes intermediate
products used to generate composite figures.

### Code

Contains code used in the analysis phase (after dada2 processing) to
generate figures and tables. Code used to generate figures are labeled
“plot_figure\_\[x\].R”. Others are intermediate or pre-processing
scripts.

## Dependencies

For the dada2 Snakemake pipeline

- See my [MiSeq](https://github.com/microlei/apprill-miseq) processing
  workflow for dependencies

For the data analysis and manuscript

- See session info below

## Session Info

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] NetSwan_0.1          VennDiagram_1.6.20   futile.logger_1.4.3 
    ##  [4] SpiecEasi_1.1.2      microbiome_1.15.3    corncob_0.2.0       
    ##  [7] graphlayouts_0.8.0   tidygraph_1.2.0      igraph_1.2.10       
    ## [10] speedyseq_0.5.3.9018 skimr_2.1.3          decontam_1.12.0     
    ## [13] vegan_2.5-7          lattice_0.20-44      permute_0.9-5       
    ## [16] phyloseq_1.36.0      ggrepel_0.9.1        RColorBrewer_1.1-2  
    ## [19] patchwork_1.1.1      cowplot_1.1.1        broom_0.7.9         
    ## [22] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.8         
    ## [25] purrr_0.3.4          readr_2.0.1          tidyr_1.1.4         
    ## [28] tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1     
    ## [31] here_1.0.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] VGAM_1.1-5             Rtsne_0.15             colorspace_2.0-2      
    ##  [4] ellipsis_0.3.2         rprojroot_2.0.2        XVector_0.32.0        
    ##  [7] base64enc_0.1-3        fs_1.5.2               rstudioapi_0.13       
    ## [10] fansi_1.0.2            lubridate_1.7.10       xml2_1.3.3            
    ## [13] codetools_0.2-18       splines_4.1.1          knitr_1.37            
    ## [16] ade4_1.7-18            jsonlite_1.7.3         cluster_2.1.2         
    ## [19] dbplyr_2.1.1           compiler_4.1.1         httr_1.4.2            
    ## [22] backports_1.2.1        assertthat_0.2.1       Matrix_1.3-4          
    ## [25] fastmap_1.1.0          cli_3.2.0              formatR_1.11          
    ## [28] htmltools_0.5.2        tools_4.1.1            gtable_0.3.0          
    ## [31] glue_1.6.1             GenomeInfoDbData_1.2.6 reshape2_1.4.4        
    ## [34] Rcpp_1.0.7             Biobase_2.52.0         cellranger_1.1.0      
    ## [37] vctrs_0.3.8            Biostrings_2.60.2      rhdf5filters_1.4.0    
    ## [40] multtest_2.48.0        ape_5.5                nlme_3.1-153          
    ## [43] iterators_1.0.13       xfun_0.29              rvest_1.0.1           
    ## [46] lifecycle_1.0.1        zlibbioc_1.38.0        MASS_7.3-54           
    ## [49] scales_1.1.1           hms_1.1.0              parallel_4.1.1        
    ## [52] biomformat_1.20.0      huge_1.3.5             rhdf5_2.36.0          
    ## [55] lambda.r_1.2.4         yaml_2.3.4             stringi_1.7.6         
    ## [58] S4Vectors_0.30.2       foreach_1.5.1          BiocGenerics_0.38.0   
    ## [61] shape_1.4.6            repr_1.1.3             GenomeInfoDb_1.28.4   
    ## [64] rlang_1.0.1            pkgconfig_2.0.3        bitops_1.0-7          
    ## [67] evaluate_0.15          Rhdf5lib_1.14.2        tidyselect_1.1.1      
    ## [70] plyr_1.8.6             magrittr_2.0.2         R6_2.5.1              
    ## [73] IRanges_2.26.0         generics_0.1.2         DBI_1.1.1             
    ## [76] pillar_1.7.0           haven_2.4.3            withr_2.4.3           
    ## [79] mgcv_1.8-36            survival_3.2-13        RCurl_1.98-1.5        
    ## [82] pulsar_0.3.7           modelr_0.1.8           crayon_1.5.0          
    ## [85] futile.options_1.0.1   utf8_1.2.2             tzdb_0.1.2            
    ## [88] rmarkdown_2.11         readxl_1.3.1           data.table_1.14.2     
    ## [91] reprex_2.0.1           digest_0.6.29          glmnet_4.1-3          
    ## [94] stats4_4.1.1           munsell_0.5.0
