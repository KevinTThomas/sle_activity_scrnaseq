Preprocessing 01: Ambient RNA removal
================
Kevin Thomas
14 January, 2022

-   [Using SoupX package to remove counts likely associated with ambient
    RNA in each
    run](#using-soupx-package-to-remove-counts-likely-associated-with-ambient-rna-in-each-run)

## Using SoupX package to remove counts likely associated with ambient RNA in each run

``` r
# variable set-up
file_directory <- "../demo_data/droplets/"
contam_gene_sets = list(
  HG = c(
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBE1",
    "HBG1",
    "HBG2",
    "HBM",
    "HBQ1",
    "HBZ"
  )
) %>% lapply(function(x) paste0("GRCh38_", x))
```

``` r
# Search file_directory (directory containing cellRanger output files) for raw and filtered barcode matrices
dataDirs <- list.files(
  path = file_directory,
  pattern = "raw_feature_bc_matrix|filtered_feature_bc_matrix|analysis",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
) %>%
  str_remove("/raw_feature_bc_matrix|/filtered_feature_bc_matrix|/analysis") %>%
  unique()
```

``` r
# Read in gene expression barcode data
sc_list <- future_map(
  .x = 1:length(dataDirs),
  .f = function (i) {
    message(str_glue("Loading {i} of {length(dataDirs)}"))
    sc <- load10X(
      dataDir = dataDirs[[i]],
      channelName = basename(dataDirs[[i]]),
      includeFeatures = c("Gene Expression")
    )
    if (mean(rowSums(sc$toc)) > 10 && mean(colSums(sc$toc)) > 100) {
      message(glue("Keeping sample {basename(dataDirs[[i]])}"))
      return(sc)
    }else{
      message(glue("Discarding sample {basename(dataDirs[[i]])}; insufficient quality"))
      sc <- NULL
      return(sc)
    }
  }
) %>% purrr::discard(.p = is.null)
```

    ## Loading 1 of 1

    ## Loading raw count data

    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.

    ## Loading cell-only count data

    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.

    ## Loading extra analysis data where available

    ## Keeping sample s1a

``` r
# Calculate the contamination using the human hemoglobin genes
contamCalcList <- future_map(
  .x = seq_along(sc_list),
  .f = function(i) {
    message(str_glue("Estimating non-expressing cells for {i} of {length(sc_list)}: {sc_list[[i]]$channelName}"))
    if (any(estimateNonExpressingCells(sc_list[[i]], contam_gene_sets))) {
      message(str_glue("Using clusters without expressing cells."))
      message(str_glue("Estimating contamination for {i} of {length(sc_list)}: {sc_list[[i]]$channelName}"))
      calculateContaminationFraction(
        sc = sc_list[[i]],
        nonExpressedGeneList = contam_gene_sets,
        useToEst = estimateNonExpressingCells(sc_list[[i]], contam_gene_sets)
      )
    }else{
      message(str_glue("No clusters found without expressing cells. Using individual cells instead."))
      message(str_glue("Estimating contamination for {i} of {length(sc_list)}: {sc_list[[i]]$channelName}"))
      calculateContaminationFraction(
        sc = sc_list[[i]],
        nonExpressedGeneList = contam_gene_sets,
        useToEst = estimateNonExpressingCells(sc_list[[i]], contam_gene_sets, clusters = FALSE)
      )
    }
  }
) %>% setNames(map_chr(sc_list, `$`, "channelName"))
```

    ## Estimating non-expressing cells for 1 of 1: s1a

    ## Using clusters without expressing cells.

    ## Estimating contamination for 1 of 1: s1a

    ## Estimated global contamination fraction of 8.61%

``` r
# Adjust counts using calculated contamination fractions
adjustCountList <- map(
  .x = contamCalcList,
  .f = adjustCounts,
  verbose = 2,
  roundToInt = TRUE
) %>% setNames(
  nm = map_chr(.x = sc_list, .f = `$`, "channelName")
  )
```

    ## Warning in sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w], : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you

    ## Expanding counts from 11 clusters to 15735 cells.

    ## Expanding cluster 1

    ## Expanding cluster 10

    ## Expanding cluster 11

    ## Expanding cluster 2

    ## Expanding cluster 3

    ## Expanding cluster 4

    ## Expanding cluster 5

    ## Expanding cluster 6

    ## Expanding cluster 7

    ## Expanding cluster 8

    ## Expanding cluster 9

    ## Rounding to integers.

``` r
# Save data for next step
saveRDS(object = adjustCountList, file = "../demo_data/adjustCountList.RDS")
```

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 4.1.0 (2021-05-18)
    ##  os       Ubuntu 20.04.2 LTS          
    ##  system   x86_64, linux-gnu           
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       Etc/UTC                     
    ##  date     2022-01-14                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  package              * version    date       lib source                                        
    ##  abind                  1.4-5      2016-07-21 [1] RSPM (R 4.1.0)                                
    ##  AnnotationDbi          1.54.1     2021-06-08 [1] RSPM (R 4.1.0)                                
    ##  assertthat             0.2.1      2019-03-21 [1] RSPM (R 4.1.0)                                
    ##  backports              1.2.1      2020-12-09 [1] RSPM (R 4.1.0)                                
    ##  base64url              1.4        2018-05-14 [1] RSPM (R 4.1.0)                                
    ##  bayesm                 3.1-4      2019-10-15 [1] RSPM (R 4.1.0)                                
    ##  beachmat               2.8.0      2021-05-19 [1] Bioconductor                                  
    ##  beeswarm               0.4.0      2021-06-01 [1] RSPM (R 4.1.0)                                
    ##  Biobase              * 2.52.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  BiocGenerics         * 0.38.0     2021-08-13 [1] bioc_git2r (@1db849a)                         
    ##  BiocManager          * 1.30.16    2021-06-15 [1] CRAN (R 4.1.0)                                
    ##  BiocNeighbors        * 1.10.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  BiocParallel         * 1.26.1     2021-07-04 [1] RSPM (R 4.1.0)                                
    ##  BiocSingular           1.8.1      2021-06-08 [1] Bioconductor                                  
    ##  Biostrings             2.60.2     2021-08-05 [1] Bioconductor                                  
    ##  bit                    4.0.4      2020-08-04 [1] RSPM (R 4.1.0)                                
    ##  bit64                  4.0.5      2020-08-30 [1] RSPM (R 4.1.0)                                
    ##  bitops                 1.0-7      2021-04-24 [1] RSPM (R 4.1.0)                                
    ##  blob                   1.2.2      2021-07-23 [1] RSPM (R 4.1.0)                                
    ##  broom                  0.7.9      2021-07-27 [1] RSPM (R 4.1.0)                                
    ##  cachem                 1.0.5      2021-05-15 [1] RSPM (R 4.1.0)                                
    ##  callr                  3.7.0      2021-04-20 [1] RSPM (R 4.1.0)                                
    ##  car                    3.0-11     2021-06-27 [1] RSPM (R 4.1.0)                                
    ##  carData                3.0-4      2020-05-22 [1] RSPM (R 4.1.0)                                
    ##  cellranger             1.1.0      2016-07-27 [1] RSPM (R 4.1.0)                                
    ##  cli                    3.0.1      2021-07-17 [1] RSPM (R 4.1.0)                                
    ##  cluster                2.1.2      2021-04-17 [2] CRAN (R 4.1.0)                                
    ##  codetools              0.2-18     2020-11-04 [2] CRAN (R 4.1.0)                                
    ##  colorspace             2.0-2      2021-06-24 [1] RSPM (R 4.1.0)                                
    ##  compositions           2.0-2      2021-07-14 [1] RSPM (R 4.1.0)                                
    ##  cowplot              * 1.1.1      2020-12-30 [1] RSPM (R 4.1.0)                                
    ##  crayon                 1.4.1      2021-02-08 [1] RSPM (R 4.1.0)                                
    ##  curl                   4.3.2      2021-06-23 [1] RSPM (R 4.1.0)                                
    ##  data.table           * 1.14.0     2021-02-21 [1] RSPM (R 4.1.0)                                
    ##  DBI                    1.1.1      2021-01-15 [1] RSPM (R 4.1.0)                                
    ##  dbplyr                 2.1.1      2021-04-06 [1] RSPM (R 4.1.0)                                
    ##  DelayedArray           0.18.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  DelayedMatrixStats     1.14.2     2021-08-08 [1] Bioconductor                                  
    ##  deldir                 0.2-10     2021-02-16 [1] RSPM (R 4.1.0)                                
    ##  deMULTIplex          * 1.0.2      2021-08-13 [1] Github (chris-mcginnis-ucsf/MULTI-seq@6e2a142)
    ##  DEoptimR               1.0-9      2021-05-24 [1] RSPM (R 4.1.0)                                
    ##  desc                   1.3.0      2021-03-05 [1] RSPM (R 4.1.0)                                
    ##  devtools               2.4.2      2021-06-07 [1] RSPM (R 4.1.0)                                
    ##  digest                 0.6.27     2020-10-24 [1] RSPM (R 4.1.0)                                
    ##  dplyr                * 1.0.7      2021-06-18 [1] RSPM (R 4.1.0)                                
    ##  drake                * 7.13.2     2021-04-22 [1] RSPM (R 4.1.0)                                
    ##  dsb                  * 0.2.0      2021-08-13 [1] Github (niaid/dsb@768691f)                    
    ##  ellipsis               0.3.2      2021-04-29 [1] RSPM (R 4.1.0)                                
    ##  evaluate               0.14       2019-05-28 [1] RSPM (R 4.1.0)                                
    ##  fansi                  0.5.0      2021-05-25 [1] RSPM (R 4.1.0)                                
    ##  farver                 2.1.0      2021-02-28 [1] RSPM (R 4.1.0)                                
    ##  fastmap                1.1.0      2021-01-25 [1] RSPM (R 4.1.0)                                
    ##  filelock               1.0.2      2018-10-05 [1] RSPM (R 4.1.0)                                
    ##  fitdistrplus           1.1-5      2021-05-28 [1] RSPM (R 4.1.0)                                
    ##  forcats              * 0.5.1      2021-01-27 [1] RSPM (R 4.1.0)                                
    ##  foreign                0.8-81     2020-12-22 [2] CRAN (R 4.1.0)                                
    ##  formattable          * 0.2.1      2021-01-07 [1] RSPM (R 4.1.0)                                
    ##  fs                     1.5.0      2020-07-31 [1] RSPM (R 4.1.0)                                
    ##  furrr                * 0.2.3      2021-06-25 [1] RSPM (R 4.1.0)                                
    ##  future               * 1.21.0     2020-12-10 [1] RSPM (R 4.1.0)                                
    ##  future.apply           1.7.0      2021-01-04 [1] RSPM (R 4.1.0)                                
    ##  future.callr         * 0.6.1      2021-05-04 [1] RSPM (R 4.1.0)                                
    ##  generics               0.1.0      2020-10-31 [1] RSPM (R 4.1.0)                                
    ##  GenomeInfoDb         * 1.28.1     2021-07-01 [1] RSPM (R 4.1.0)                                
    ##  GenomeInfoDbData       1.2.6      2021-08-04 [1] RSPM (R 4.1.0)                                
    ##  GenomicRanges        * 1.44.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  ggbeeswarm             0.6.0      2017-08-07 [1] RSPM (R 4.1.0)                                
    ##  ggforce              * 0.3.3      2021-03-05 [1] RSPM (R 4.1.0)                                
    ##  ggplot2              * 3.3.5      2021-06-25 [1] RSPM (R 4.1.0)                                
    ##  ggpubr                 0.4.0      2020-06-27 [1] RSPM (R 4.1.0)                                
    ##  ggrepel              * 0.9.1      2021-01-15 [1] RSPM (R 4.1.0)                                
    ##  ggridges               0.5.3      2021-01-08 [1] RSPM (R 4.1.0)                                
    ##  ggsignif               0.6.2      2021-06-14 [1] RSPM (R 4.1.0)                                
    ##  globals                0.14.0     2020-11-22 [1] RSPM (R 4.1.0)                                
    ##  glue                 * 1.4.2      2020-08-27 [1] RSPM (R 4.1.0)                                
    ##  goftest                1.2-2      2019-12-02 [1] RSPM (R 4.1.0)                                
    ##  gridExtra              2.3        2017-09-09 [1] RSPM (R 4.1.0)                                
    ##  gtable                 0.3.0      2019-03-25 [1] RSPM (R 4.1.0)                                
    ##  gtools                 3.9.2      2021-06-06 [1] RSPM (R 4.1.0)                                
    ##  harmony              * 0.1.0      2021-08-13 [1] Github (immunogenomics/harmony@c93de54)       
    ##  haven                  2.4.3      2021-08-04 [1] RSPM (R 4.1.0)                                
    ##  HGNChelper             0.8.1      2019-10-24 [1] RSPM (R 4.1.0)                                
    ##  hms                    1.1.0      2021-05-17 [1] RSPM (R 4.1.0)                                
    ##  htmltools              0.5.1.1    2021-01-22 [1] RSPM (R 4.1.0)                                
    ##  htmlwidgets            1.5.3      2020-12-10 [1] RSPM (R 4.1.0)                                
    ##  httpuv                 1.6.1      2021-05-07 [1] RSPM (R 4.1.0)                                
    ##  httr                   1.4.2      2020-07-20 [1] RSPM (R 4.1.0)                                
    ##  ica                    1.0-2      2018-05-24 [1] RSPM (R 4.1.0)                                
    ##  igraph                 1.2.6      2020-10-06 [1] RSPM (R 4.1.0)                                
    ##  IRanges              * 2.26.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  irlba                  2.3.3      2019-02-05 [1] RSPM (R 4.1.0)                                
    ##  janitor              * 2.1.0      2021-01-05 [1] RSPM (R 4.1.0)                                
    ##  jsonlite               1.7.2      2020-12-09 [1] RSPM (R 4.1.0)                                
    ##  kableExtra           * 1.3.4      2021-02-20 [1] RSPM (R 4.1.0)                                
    ##  KEGGREST               1.32.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  KernSmooth             2.23-20    2021-05-03 [2] CRAN (R 4.1.0)                                
    ##  knitr                  1.33       2021-04-24 [1] RSPM (R 4.1.0)                                
    ##  later                  1.2.0      2021-04-23 [1] RSPM (R 4.1.0)                                
    ##  lattice                0.20-44    2021-05-02 [2] CRAN (R 4.1.0)                                
    ##  lazyeval               0.2.2      2019-03-15 [1] RSPM (R 4.1.0)                                
    ##  leiden                 0.3.9      2021-07-27 [1] RSPM (R 4.1.0)                                
    ##  lifecycle              1.0.0      2021-02-15 [1] RSPM (R 4.1.0)                                
    ##  limma                  3.48.3     2021-08-10 [1] Bioconductor                                  
    ##  listenv                0.8.0      2019-12-05 [1] RSPM (R 4.1.0)                                
    ##  lmtest                 0.9-38     2020-09-09 [1] RSPM (R 4.1.0)                                
    ##  lubridate            * 1.7.10     2021-02-26 [1] RSPM (R 4.1.0)                                
    ##  magrittr             * 2.0.1      2020-11-17 [1] RSPM (R 4.1.0)                                
    ##  MASS                   7.3-54     2021-05-03 [2] CRAN (R 4.1.0)                                
    ##  Matrix               * 1.3-4      2021-06-01 [2] RSPM (R 4.1.0)                                
    ##  MatrixGenerics       * 1.4.2      2021-08-08 [1] Bioconductor                                  
    ##  matrixStats          * 0.60.0     2021-07-26 [1] RSPM (R 4.1.0)                                
    ##  mclust                 5.4.7      2020-11-20 [1] RSPM (R 4.1.0)                                
    ##  memoise                2.0.0      2021-01-26 [1] RSPM (R 4.1.0)                                
    ##  mgcv                   1.8-36     2021-06-01 [2] RSPM (R 4.1.0)                                
    ##  mime                   0.11       2021-06-23 [1] RSPM (R 4.1.0)                                
    ##  miniUI                 0.1.1.1    2018-05-18 [1] RSPM (R 4.1.0)                                
    ##  modelr                 0.1.8      2020-05-19 [1] RSPM (R 4.1.0)                                
    ##  munsell                0.5.0      2018-06-12 [1] RSPM (R 4.1.0)                                
    ##  nlme                   3.1-152    2021-02-04 [2] CRAN (R 4.1.0)                                
    ##  openxlsx               4.2.4      2021-06-16 [1] RSPM (R 4.1.0)                                
    ##  org.Hs.eg.db           3.13.0     2021-08-13 [1] Bioconductor                                  
    ##  paletteer            * 1.4.0      2021-07-20 [1] RSPM (R 4.1.0)                                
    ##  parallelly             1.27.0     2021-07-19 [1] RSPM (R 4.1.0)                                
    ##  patchwork              1.1.1      2020-12-17 [1] RSPM (R 4.1.0)                                
    ##  pbapply                1.4-3      2020-08-18 [1] RSPM (R 4.1.0)                                
    ##  pillar                 1.6.2      2021-07-29 [1] RSPM (R 4.1.0)                                
    ##  pkgbuild               1.2.0      2020-12-15 [1] RSPM (R 4.1.0)                                
    ##  pkgconfig              2.0.3      2019-09-22 [1] RSPM (R 4.1.0)                                
    ##  pkgload                1.2.1      2021-04-06 [1] RSPM (R 4.1.0)                                
    ##  plotly                 4.9.4.1    2021-06-18 [1] RSPM (R 4.1.0)                                
    ##  plyr                   1.8.6      2020-03-03 [1] RSPM (R 4.1.0)                                
    ##  png                    0.1-7      2013-12-03 [1] RSPM (R 4.1.0)                                
    ##  polyclip               1.10-0     2019-03-14 [1] RSPM (R 4.1.0)                                
    ##  presto               * 1.0.0      2021-08-13 [1] Github (immunogenomics/presto@052085d)        
    ##  prettyunits            1.1.1      2020-01-24 [1] RSPM (R 4.1.0)                                
    ##  processx               3.5.2      2021-04-30 [1] RSPM (R 4.1.0)                                
    ##  progress               1.2.2      2019-05-16 [1] RSPM (R 4.1.0)                                
    ##  promises               1.2.0.1    2021-02-11 [1] RSPM (R 4.1.0)                                
    ##  ps                     1.6.0      2021-02-28 [1] RSPM (R 4.1.0)                                
    ##  purrr                * 0.3.4      2020-04-17 [1] RSPM (R 4.1.0)                                
    ##  R6                     2.5.0      2020-10-28 [1] RSPM (R 4.1.0)                                
    ##  RANN                   2.6.1      2019-01-08 [1] RSPM (R 4.1.0)                                
    ##  RColorBrewer         * 1.1-2      2014-12-07 [1] RSPM (R 4.1.0)                                
    ##  Rcpp                 * 1.0.7      2021-07-07 [1] RSPM (R 4.1.0)                                
    ##  RcppAnnoy              0.0.19     2021-07-30 [1] RSPM (R 4.1.0)                                
    ##  RCurl                  1.98-1.3   2021-03-16 [1] RSPM (R 4.1.0)                                
    ##  readr                * 2.0.0      2021-07-20 [1] RSPM (R 4.1.0)                                
    ##  readxl                 1.3.1      2019-03-13 [1] RSPM (R 4.1.0)                                
    ##  ReductionWrappers    * 2.5.3      2021-08-13 [1] Github (milescsmith/ReductionWrappers@cb0f3a8)
    ##  rematch2               2.1.2      2020-05-01 [1] RSPM (R 4.1.0)                                
    ##  remotes                2.4.0      2021-06-02 [1] RSPM (R 4.1.0)                                
    ##  reprex                 2.0.1      2021-08-05 [1] RSPM (R 4.1.0)                                
    ##  reshape2               1.4.4      2020-04-09 [1] RSPM (R 4.1.0)                                
    ##  reticulate           * 1.20       2021-05-03 [1] RSPM (R 4.1.0)                                
    ##  rio                    0.5.27     2021-06-21 [1] RSPM (R 4.1.0)                                
    ##  rlang                * 0.4.11     2021-04-30 [1] RSPM (R 4.1.0)                                
    ##  rmarkdown              2.10       2021-08-06 [1] RSPM (R 4.1.0)                                
    ##  robustbase             0.93-8     2021-06-02 [1] RSPM (R 4.1.0)                                
    ##  ROCR                   1.0-11     2020-05-02 [1] RSPM (R 4.1.0)                                
    ##  rpart                  4.1-15     2019-04-12 [2] CRAN (R 4.1.0)                                
    ##  rprojroot              2.0.2      2020-11-15 [1] RSPM (R 4.1.0)                                
    ##  RSQLite                2.2.7      2021-04-22 [1] RSPM (R 4.1.0)                                
    ##  rstatix                0.7.0      2021-02-13 [1] RSPM (R 4.1.0)                                
    ##  rstudioapi             0.13       2020-11-12 [1] RSPM (R 4.1.0)                                
    ##  rsvd                 * 1.0.5      2021-04-16 [1] RSPM (R 4.1.0)                                
    ##  Rtsne                  0.15       2018-11-10 [1] RSPM (R 4.1.0)                                
    ##  rvest                  1.0.1      2021-07-26 [1] RSPM (R 4.1.0)                                
    ##  s2a                  * 0.3.2      2021-08-13 [1] Github (milescsmith/s2a@0c729f4)              
    ##  S4Vectors            * 0.30.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  ScaledMatrix           1.0.0      2021-05-19 [1] Bioconductor                                  
    ##  scales                 1.1.1      2020-05-11 [1] RSPM (R 4.1.0)                                
    ##  scater               * 1.20.1     2021-06-15 [1] Bioconductor                                  
    ##  scattermore            0.7        2020-11-24 [1] RSPM (R 4.1.0)                                
    ##  scClustViz           * 1.3.8      2021-08-13 [1] Github (BaderLab/scClustViz@32064ce)          
    ##  sctransform            0.3.2      2020-12-16 [1] RSPM (R 4.1.0)                                
    ##  scuttle              * 1.2.1      2021-08-05 [1] Bioconductor                                  
    ##  sessioninfo            1.1.1      2018-11-05 [1] RSPM (R 4.1.0)                                
    ##  Seurat               * 4.0.3      2021-08-13 [1] Github (satijalab/Seurat@9b38929)             
    ##  SeuratBubblePlot     * 0.4.2      2021-08-13 [1] Github (milescsmith/SeuratBubblePlot@e0f2abf) 
    ##  SeuratObject         * 4.0.2      2021-06-09 [1] RSPM (R 4.1.0)                                
    ##  shiny                * 1.6.0      2021-01-25 [1] RSPM (R 4.1.0)                                
    ##  SingleCellExperiment * 1.14.1     2021-08-13 [1] bioc_git2r (@5357eff)                         
    ##  SingleR              * 1.6.1      2021-05-20 [1] Bioconductor                                  
    ##  snakecase              0.11.0     2019-05-25 [1] RSPM (R 4.1.0)                                
    ##  SoupX                * 1.5.2      2021-05-17 [1] RSPM (R 4.1.0)                                
    ##  sparseMatrixStats      1.4.2      2021-08-08 [1] Bioconductor                                  
    ##  spatstat.core          2.3-0      2021-07-16 [1] RSPM (R 4.1.0)                                
    ##  spatstat.data          2.1-0      2021-03-21 [1] RSPM (R 4.1.0)                                
    ##  spatstat.geom          2.2-2      2021-07-12 [1] RSPM (R 4.1.0)                                
    ##  spatstat.sparse        2.0-0      2021-03-16 [1] RSPM (R 4.1.0)                                
    ##  spatstat.utils         2.2-0      2021-06-14 [1] RSPM (R 4.1.0)                                
    ##  storr                  1.2.5      2020-12-01 [1] RSPM (R 4.1.0)                                
    ##  stringi                1.7.3      2021-07-16 [1] RSPM (R 4.1.0)                                
    ##  stringr              * 1.4.0      2019-02-10 [1] RSPM (R 4.1.0)                                
    ##  SummarizedExperiment * 1.22.0     2021-08-13 [1] bioc_git2r (@7d1110e)                         
    ##  survival               3.2-11     2021-04-26 [2] CRAN (R 4.1.0)                                
    ##  svglite                2.0.0      2021-02-20 [1] RSPM (R 4.1.0)                                
    ##  systemfonts            1.0.2      2021-05-11 [1] RSPM (R 4.1.0)                                
    ##  tensor                 1.5        2012-05-05 [1] RSPM (R 4.1.0)                                
    ##  tensorA                0.36.2     2020-11-19 [1] RSPM (R 4.1.0)                                
    ##  testthat               3.0.4      2021-07-01 [1] RSPM (R 4.1.0)                                
    ##  tibble               * 3.1.3      2021-07-23 [1] RSPM (R 4.1.0)                                
    ##  tidyr                * 1.1.3      2021-03-03 [1] RSPM (R 4.1.0)                                
    ##  tidyselect             1.1.1      2021-04-30 [1] RSPM (R 4.1.0)                                
    ##  tidyverse            * 1.3.1      2021-04-15 [1] RSPM (R 4.1.0)                                
    ##  tweenr                 1.0.2      2021-03-23 [1] RSPM (R 4.1.0)                                
    ##  txtq                   0.2.4      2021-03-27 [1] RSPM (R 4.1.0)                                
    ##  tzdb                   0.1.2      2021-07-20 [1] RSPM (R 4.1.0)                                
    ##  usethis                2.0.1      2021-02-10 [1] RSPM (R 4.1.0)                                
    ##  utf8                   1.2.2      2021-07-24 [1] RSPM (R 4.1.0)                                
    ##  uwot                   0.1.10     2020-12-15 [1] RSPM (R 4.1.0)                                
    ##  vctrs                  0.3.8      2021-04-29 [1] RSPM (R 4.1.0)                                
    ##  ViolinEnsemble       * 0.0.0.9300 2021-08-13 [1] Github (milescsmith/ViolinEnsemble@f383f08)   
    ##  vipor                  0.4.5      2017-03-22 [1] RSPM (R 4.1.0)                                
    ##  viridis              * 0.6.1      2021-05-11 [1] RSPM (R 4.1.0)                                
    ##  viridisLite          * 0.4.0      2021-04-13 [1] RSPM (R 4.1.0)                                
    ##  webshot                0.5.2      2019-11-22 [1] RSPM (R 4.1.0)                                
    ##  withr                  2.4.2      2021-04-18 [1] RSPM (R 4.1.0)                                
    ##  xfun                   0.25       2021-08-06 [1] RSPM (R 4.1.0)                                
    ##  xml2                 * 1.3.2      2020-04-23 [1] RSPM (R 4.1.0)                                
    ##  xtable                 1.8-4      2019-04-21 [1] RSPM (R 4.1.0)                                
    ##  XVector                0.32.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  yaml                   2.2.1      2020-02-01 [1] RSPM (R 4.1.0)                                
    ##  zip                    2.2.0      2021-05-31 [1] RSPM (R 4.1.0)                                
    ##  zlibbioc               1.38.0     2021-05-19 [1] RSPM (R 4.1.0)                                
    ##  zoo                    1.8-9      2021-03-09 [1] RSPM (R 4.1.0)                                
    ## 
    ## [1] /usr/local/lib/R/site-library
    ## [2] /usr/local/lib/R/library
