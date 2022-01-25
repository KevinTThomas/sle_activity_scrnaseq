Preprocessing 01: Ambient RNA removal
================
Kevin Thomas
20 January, 2022

-   [Using SoupX package to remove counts likely associated with ambient
    RNA in each
    run](#using-soupx-package-to-remove-counts-likely-associated-with-ambient-rna-in-each-run)

## Using SoupX package to remove counts likely associated with ambient RNA in each run

``` r
# variable set-up
file_directory <- here("demo_data", "droplets")
contam_gene_sets <- list(
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
) |> lapply(function(x) paste0("GRCh38_", x))
```

``` r
# Search file_directory (directory containing cellRanger output files) for raw and filtered barcode matrices
dataDirs <- list.files(
  path = file_directory,
  pattern = "raw_feature_bc_matrix|filtered_feature_bc_matrix|analysis",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
) |>
  str_remove("/raw_feature_bc_matrix|/filtered_feature_bc_matrix|/analysis") |>
  unique()
```

``` r
# Read in gene expression barcode data
sc_list <- furrr::future_imap(
  .x = dataDirs,
  .f = function(i, j) {
    # message(stringr::str_glue("Loading {j} of {length(dataDirs)}"))
    message(i)
    sc <- load10X(
      dataDir = i,
      channelName = basename(i),
      includeFeatures = c("Gene Expression")
    )
    if (mean(rowSums(sc$toc)) > 10 && mean(colSums(sc$toc)) > 100) {
      message(str_glue("Keeping sample {basename(i)}"))
      return(sc)
    } else {
      message(str_glue("Discarding sample {basename(i)}; insufficient quality"))
      sc <- NULL
      return(sc)
    }
  }
) |>
  discard(.p = is.null)
```

    ## /home/rstudio/workspace/sle_activity_scrnaseq/demo_data/droplets/s1a

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
    } else {
      message(str_glue("No clusters found without expressing cells. Using individual cells instead."))
      message(str_glue("Estimating contamination for {i} of {length(sc_list)}: {sc_list[[i]]$channelName}"))
      calculateContaminationFraction(
        sc = sc_list[[i]],
        nonExpressedGeneList = contam_gene_sets,
        useToEst = estimateNonExpressingCells(sc_list[[i]], contam_gene_sets, clusters = FALSE)
      )
    }
  }
) |> set_names(map_chr(sc_list, use_series, "channelName"))
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
) |> set_names(
  nm = map_chr(.x = sc_list, .f = use_series, "channelName")
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
saveRDS(object = adjustCountList, file = here("demo_data", "adjustCountList.RDS"))
```

``` r
session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.1.2 (2021-11-01)
    ##  os       Ubuntu 20.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Etc/UTC
    ##  date     2022-01-20
    ##  pandoc   2.14.0.3 @ /usr/lib/rstudio-server/bin/pandoc/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  package         * version date (UTC) lib source
    ##  abind             1.4-5   2016-07-21 [1] RSPM (R 4.1.0)
    ##  BiocManager     * 1.30.16 2021-06-15 [1] RSPM (R 4.1.0)
    ##  BiocParallel    * 1.28.3  2021-12-09 [1] Bioconductor
    ##  cachem            1.0.6   2021-08-19 [1] RSPM (R 4.1.0)
    ##  callr             3.7.0   2021-04-20 [1] RSPM (R 4.1.0)
    ##  cli               3.1.0   2021-10-27 [1] RSPM (R 4.1.0)
    ##  cluster           2.1.2   2021-04-17 [2] CRAN (R 4.1.2)
    ##  codetools         0.2-18  2020-11-04 [2] CRAN (R 4.1.2)
    ##  colorspace        2.0-2   2021-06-24 [1] RSPM (R 4.1.0)
    ##  cowplot           1.1.1   2020-12-30 [1] RSPM (R 4.1.0)
    ##  crayon            1.4.2   2021-10-29 [1] RSPM (R 4.1.0)
    ##  data.table        1.14.2  2021-09-27 [1] RSPM (R 4.1.0)
    ##  deldir            1.0-6   2021-10-23 [1] RSPM (R 4.1.0)
    ##  desc              1.4.0   2021-09-28 [1] RSPM (R 4.1.0)
    ##  devtools        * 2.4.3   2021-11-30 [1] RSPM (R 4.1.0)
    ##  digest            0.6.29  2021-12-01 [1] RSPM (R 4.1.0)
    ##  dplyr           * 1.0.7   2021-06-18 [1] RSPM (R 4.1.0)
    ##  ellipsis          0.3.2   2021-04-29 [1] RSPM (R 4.1.0)
    ##  evaluate          0.14    2019-05-28 [1] RSPM (R 4.1.0)
    ##  fansi             1.0.2   2022-01-14 [1] RSPM (R 4.1.0)
    ##  fastmap           1.1.0   2021-01-25 [1] RSPM (R 4.1.0)
    ##  fitdistrplus      1.1-6   2021-09-28 [1] RSPM (R 4.1.0)
    ##  fs                1.5.2   2021-12-08 [1] RSPM (R 4.1.0)
    ##  furrr           * 0.2.3   2021-06-25 [1] RSPM (R 4.1.0)
    ##  future          * 1.23.0  2021-10-31 [1] RSPM (R 4.1.0)
    ##  future.apply      1.8.1   2021-08-10 [1] RSPM (R 4.1.0)
    ##  generics          0.1.1   2021-10-25 [1] RSPM (R 4.1.0)
    ##  ggplot2           3.3.5   2021-06-25 [1] RSPM (R 4.1.0)
    ##  ggrepel           0.9.1   2021-01-15 [1] RSPM (R 4.1.0)
    ##  ggridges          0.5.3   2021-01-08 [1] RSPM (R 4.1.0)
    ##  globals           0.14.0  2020-11-22 [1] RSPM (R 4.1.0)
    ##  glue              1.6.0   2021-12-17 [1] RSPM (R 4.1.0)
    ##  goftest           1.2-3   2021-10-07 [1] RSPM (R 4.1.0)
    ##  gridExtra         2.3     2017-09-09 [1] RSPM (R 4.1.0)
    ##  gtable            0.3.0   2019-03-25 [1] RSPM (R 4.1.0)
    ##  here            * 1.0.1   2020-12-13 [1] RSPM (R 4.1.0)
    ##  htmltools         0.5.2   2021-08-25 [1] RSPM (R 4.1.0)
    ##  htmlwidgets       1.5.4   2021-09-08 [1] RSPM (R 4.1.0)
    ##  httpuv            1.6.5   2022-01-05 [1] RSPM (R 4.1.0)
    ##  httr              1.4.2   2020-07-20 [1] RSPM (R 4.1.0)
    ##  ica               1.0-2   2018-05-24 [1] RSPM (R 4.1.0)
    ##  igraph            1.2.11  2022-01-04 [1] RSPM (R 4.1.0)
    ##  irlba             2.3.5   2021-12-06 [1] RSPM (R 4.1.0)
    ##  jsonlite          1.7.3   2022-01-17 [1] RSPM (R 4.1.0)
    ##  KernSmooth        2.23-20 2021-05-03 [2] CRAN (R 4.1.2)
    ##  knitr             1.37    2021-12-16 [1] RSPM (R 4.1.0)
    ##  later             1.3.0   2021-08-18 [1] RSPM (R 4.1.0)
    ##  lattice           0.20-45 2021-09-22 [2] CRAN (R 4.1.2)
    ##  lazyeval          0.2.2   2019-03-15 [1] RSPM (R 4.1.0)
    ##  leiden            0.3.9   2021-07-27 [1] RSPM (R 4.1.0)
    ##  lifecycle         1.0.1   2021-09-24 [1] RSPM (R 4.1.0)
    ##  listenv           0.8.0   2019-12-05 [1] RSPM (R 4.1.0)
    ##  lmtest            0.9-39  2021-11-07 [1] RSPM (R 4.1.0)
    ##  magrittr        * 2.0.1   2020-11-17 [1] RSPM (R 4.1.0)
    ##  MASS              7.3-54  2021-05-03 [2] CRAN (R 4.1.2)
    ##  Matrix          * 1.3-4   2021-06-01 [2] CRAN (R 4.1.2)
    ##  matrixStats       0.61.0  2021-09-17 [1] RSPM (R 4.1.0)
    ##  memoise           2.0.1   2021-11-26 [1] RSPM (R 4.1.0)
    ##  mgcv              1.8-38  2021-10-06 [2] CRAN (R 4.1.2)
    ##  mime              0.12    2021-09-28 [1] RSPM (R 4.1.0)
    ##  miniUI            0.1.1.1 2018-05-18 [1] RSPM (R 4.1.0)
    ##  munsell           0.5.0   2018-06-12 [1] RSPM (R 4.1.0)
    ##  nlme              3.1-153 2021-09-07 [2] CRAN (R 4.1.2)
    ##  parallelly        1.30.0  2021-12-17 [1] RSPM (R 4.1.0)
    ##  patchwork         1.1.1   2020-12-17 [1] RSPM (R 4.1.0)
    ##  pbapply           1.5-0   2021-09-16 [1] RSPM (R 4.1.0)
    ##  pillar            1.6.4   2021-10-18 [1] RSPM (R 4.1.0)
    ##  pkgbuild          1.3.1   2021-12-20 [1] RSPM (R 4.1.0)
    ##  pkgconfig         2.0.3   2019-09-22 [1] RSPM (R 4.1.0)
    ##  pkgload           1.2.4   2021-11-30 [1] RSPM (R 4.1.0)
    ##  plotly            4.10.0  2021-10-09 [1] RSPM (R 4.1.0)
    ##  plyr              1.8.6   2020-03-03 [1] RSPM (R 4.1.0)
    ##  png               0.1-7   2013-12-03 [1] RSPM (R 4.1.0)
    ##  polyclip          1.10-0  2019-03-14 [1] RSPM (R 4.1.0)
    ##  prettyunits       1.1.1   2020-01-24 [1] RSPM (R 4.1.0)
    ##  processx          3.5.2   2021-04-30 [1] RSPM (R 4.1.0)
    ##  promises          1.2.0.1 2021-02-11 [1] RSPM (R 4.1.0)
    ##  ps                1.6.0   2021-02-28 [1] RSPM (R 4.1.0)
    ##  purrr           * 0.3.4   2020-04-17 [1] RSPM (R 4.1.0)
    ##  R6                2.5.1   2021-08-19 [1] RSPM (R 4.1.0)
    ##  RANN              2.6.1   2019-01-08 [1] RSPM (R 4.1.0)
    ##  RColorBrewer      1.1-2   2014-12-07 [1] RSPM (R 4.1.0)
    ##  Rcpp              1.0.8   2022-01-13 [1] RSPM (R 4.1.0)
    ##  RcppAnnoy         0.0.19  2021-07-30 [1] RSPM (R 4.1.0)
    ##  remotes           2.4.2   2021-11-30 [1] RSPM (R 4.1.0)
    ##  reshape2          1.4.4   2020-04-09 [1] RSPM (R 4.1.0)
    ##  reticulate        1.23    2022-01-14 [1] RSPM (R 4.1.0)
    ##  rlang           * 0.4.12  2021-10-18 [1] RSPM (R 4.1.0)
    ##  rmarkdown         2.11    2021-09-14 [1] RSPM (R 4.1.0)
    ##  ROCR              1.0-11  2020-05-02 [1] RSPM (R 4.1.0)
    ##  rpart             4.1-15  2019-04-12 [2] CRAN (R 4.1.2)
    ##  rprojroot         2.0.2   2020-11-15 [1] RSPM (R 4.1.0)
    ##  rstudioapi        0.13    2020-11-12 [1] RSPM (R 4.1.0)
    ##  Rtsne             0.15    2018-11-10 [1] RSPM (R 4.1.0)
    ##  scales            1.1.1   2020-05-11 [1] RSPM (R 4.1.0)
    ##  scattermore       0.7     2020-11-24 [1] RSPM (R 4.1.0)
    ##  sctransform       0.3.3   2022-01-13 [1] RSPM (R 4.1.0)
    ##  sessioninfo       1.2.2   2021-12-06 [1] RSPM (R 4.1.0)
    ##  Seurat            4.1.0   2022-01-14 [1] RSPM (R 4.1.0)
    ##  SeuratObject      4.0.4   2021-11-23 [1] RSPM (R 4.1.0)
    ##  shiny             1.7.1   2021-10-02 [1] RSPM (R 4.1.0)
    ##  SoupX           * 1.5.2   2021-05-17 [1] RSPM (R 4.1.0)
    ##  spatstat.core     2.3-2   2021-11-26 [1] RSPM (R 4.1.0)
    ##  spatstat.data     2.1-2   2021-12-17 [1] RSPM (R 4.1.0)
    ##  spatstat.geom     2.3-1   2021-12-10 [1] RSPM (R 4.1.0)
    ##  spatstat.sparse   2.1-0   2021-12-17 [1] RSPM (R 4.1.0)
    ##  spatstat.utils    2.3-0   2021-12-12 [1] RSPM (R 4.1.0)
    ##  stringi           1.7.6   2021-11-29 [1] RSPM (R 4.1.0)
    ##  stringr         * 1.4.0   2019-02-10 [1] RSPM (R 4.1.0)
    ##  survival          3.2-13  2021-08-24 [2] CRAN (R 4.1.2)
    ##  tensor            1.5     2012-05-05 [1] RSPM (R 4.1.0)
    ##  testthat          3.1.1   2021-12-03 [1] RSPM (R 4.1.0)
    ##  tibble            3.1.6   2021-11-07 [1] RSPM (R 4.1.0)
    ##  tidyr             1.1.4   2021-09-27 [1] RSPM (R 4.1.0)
    ##  tidyselect        1.1.1   2021-04-30 [1] RSPM (R 4.1.0)
    ##  usethis         * 2.1.5   2021-12-09 [1] RSPM (R 4.1.0)
    ##  utf8              1.2.2   2021-07-24 [1] RSPM (R 4.1.0)
    ##  uwot              0.1.11  2021-12-02 [1] RSPM (R 4.1.0)
    ##  vctrs             0.3.8   2021-04-29 [1] RSPM (R 4.1.0)
    ##  viridisLite       0.4.0   2021-04-13 [1] RSPM (R 4.1.0)
    ##  withr             2.4.3   2021-11-30 [1] RSPM (R 4.1.0)
    ##  xfun              0.29    2021-12-14 [1] RSPM (R 4.1.0)
    ##  xtable            1.8-4   2019-04-21 [1] RSPM (R 4.1.0)
    ##  yaml              2.2.1   2020-02-01 [1] RSPM (R 4.1.0)
    ##  zoo               1.8-9   2021-03-09 [1] RSPM (R 4.1.0)
    ## 
    ##  [1] /usr/local/lib/R/site-library
    ##  [2] /usr/local/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
