Preprocessing 02: Multimodal object construction
================
Kevin Thomas
20 January, 2022

-   [Construct multimodal object with filtered counts from previous
    step](#construct-multimodal-object-with-filtered-counts-from-previous-step)

## Construct multimodal object with filtered counts from previous step

``` r
# Load data from last step and samplesheet for object metadata
# Here, gene expression matrices are simply taken as the adjusted counts output of the SoupX analysis
gem_matrices <- readRDS(file = here("demo_data", "adjustCountList.RDS"))
samplesheet <- read_csv(file = here("demo_data", "da_samplesheet_final.csv")) |>
  clean_names()
```

    ## Rows: 46 Columns: 6

    ## ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (5): Subject_id, ancestry, classification, run, Hashtag
    ## dbl (1): age

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# We still need to find cite counts, if they exist, but only the ones we've carried through to this point
filtered_matrices <- map(
  .x = here(dataDirs, "filtered_feature_bc_matrix"),
  .f = Read10X
)
```

    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.

``` r
# Build adt_matrices if filtered_matrices contain "Antibody Capture" data and the sample is included in gem_matrices
adt_matrices <- if_else(
  condition = map2_lgl(
    .x = filtered_matrices,
    .y = basename(dataDirs),
    .f = function(x, y) {
      "Antibody Capture" %in% names(x) && y %in% names(gem_matrices)
    }
  ),
  true = filtered_matrices %>% map("Antibody Capture"),
  false = NULL
) |>
  discard(.p = is.null) |>
  set_names(names(gem_matrices))

# Free memory
rm(filtered_matrices)
gc()
```

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  3795776 202.8    6866628 366.8  4865037 259.9
    ## Vcells 21131171 161.3   84446138 644.3 89831369 685.4

``` r
# Separate mouse from human cells in gene expression
gem_classification <- calcGEMClassification(gem_matrices, primary_species_prefix, secondary_species_prefix)
```

    ## Joining, by = c("id", "run", "cells")

``` r
# DSB-normalize each ADT matrix, but only if CITE-seq data is sufficient
dsb_adt_matrices <- if (!is.null(adt_matrices)) {
  future_map(
    .x = seq_along(adt_matrices),
    .f = function(i) {
      message(str_glue("Evaluating CITE counts for sample {names(adt_matrices)[[i]]} of {length(adt_matrices)}."))
      if (mean(rowSums(adt_matrices[[i]])) > 50 && mean(colSums(adt_matrices[[i]])) > 10) {
        message(str_glue("Keeping sample {names(adt_matrices)[[i]]}."))

        cells_use <-
          gem_classification |>
          filter(run == names(x = adt_matrices)[[i]]) |>
          filter(call == str_remove_all(string = secondary_species_prefix, pattern = "_")) |>
          pull(cells)

        adt <- DSBNormalizeProtein(
          cell_protein_matrix = adt_matrices[[i]],
          empty_drop_matrix =
            adt_matrices[[i]][, cells_use],
          use.isotype.control = FALSE,
          isotype.control.name.vec = NULL
        )
        adt
      } else {
        message(str_glue("Not normalizing sample {names(adt_matrices)[[i]]}. CITE data insufficient."))
        adt <- adt_matrices[[i]]
        adt
      }
    }
  ) |>
    setNames(names(adt_matrices)) |>
    discard(.p = is.null)
} else {
  NULL
}
```

    ## [1] "potential isotype controls detected: "
    ## [1] "isotype-control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s1a of 1.

    ## Keeping sample s1a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE not recommended if isotype controls are available.
    ##  If data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE
    ##  and set `isotype.control.name.vec` to a vector of isotype control rownames from cell_protein_matrix

``` r
# Construct each run as a multimodal seurat object, stored in a list
## DO NOT FUTURE_MAP THIS STEP (parallelization is already happening in the `constructGEMs` function)
gems <- map(
  .x = seq_along(gem_matrices[names(dsb_adt_matrices)]),
  .f = function(i) {
    message(str_glue("Constructing seurat object for run {names(gem_matrices[names(dsb_adt_matrices)])[[i]]}"))
    constructGEMs(
      sample = names(gem_matrices[names(dsb_adt_matrices)])[[i]],
      gem_matrix = gem_matrices[names(dsb_adt_matrices)][[i]],
      adt_matrix = adt_matrices[names(dsb_adt_matrices)][[i]],
      norm_adt_matrix = dsb_adt_matrices[[i]],
      hashtags = if (multisample == "hash") {
        levels(factor(samplesheet$hashtag))
      } else {
        FALSE
      }
    )
  }
) |>
  set_names(names(gem_matrices[names(dsb_adt_matrices)])) |>
  discard(.p = is.null)
```

    ## Constructing seurat object for run s1a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 13180

    ## Cutoff for HT-1: 16 at quantile = 0.999

    ## Cutoff for HT-2: 25 at quantile = 0.9

    ## Cutoff for HT-3: 50 at quantile = 0.85

    ## Cutoff for HT-4: 202 at quantile = 0.99

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Loading required package: KernSmooth

    ## KernSmooth 2.23 loaded
    ## Copyright M. P. Wand 1997-2009

    ## Loading required package: reshape2

    ## 
    ## Attaching package: 'reshape2'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     dcast, melt

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

    ## Round 1: 755 Negative cells removed

    ## Round 2: 1 Negative cells removed

    ## Round 3: Done! All cells classified.

    ## HT-1 cells: 2085

    ## HT-2 cells: 3441

    ## HT-3 cells: 3343

    ## HT-4 cells: 4266

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

``` r
# Add in the metadata
gems_meta <- future_map(
  .x = seq_along(gems),
  .f = function(i) {
    x <- gems[[i]]
    assignMetaData(
      seurat_obj = x,
      samplesheet = samplesheet,
      multiplex = multisample,
      barcode_df = NULL
    )
  }
) |>
  set_names(names(gems)) |>
  discard(.p = is.null)
```

    ## Registered S3 method overwritten by 'cli':
    ##   method     from         
    ##   print.boxx spatstat.geom

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(x)` instead of `x` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

``` r
# Free memory
rm(gems)
gc()
```

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  4163010 222.4    6866628 366.8  6866628 366.8
    ## Vcells 52597927 401.3   84446138 644.3 89831369 685.4

``` r
# Merge data together
pre_qc <- if (length(gems_meta) > 1) {
  merge(
    gems_meta[[1]],
    map(
      2:length(gems_meta),
      function(x) gems_meta[[x]]
    )
  )
} else {
  gems_meta[[1]]
}

# Free memory
rm(gems_meta)
gc()
```

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  4162990 222.4    6866628 366.8  6866628 366.8
    ## Vcells 52597910 401.3   84446138 644.3 89831369 685.4

``` r
# Save data for next step
saveRDS(object = pre_qc, file = here("demo_data", "pre_qc.RDS"))
saveRDS(object = gem_classification, file = here("demo_data", "gem_classification.RDS"))
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
    ##  BiocGenerics      0.40.0  2021-10-26 [1] Bioconductor
    ##  BiocManager     * 1.30.16 2021-06-15 [1] RSPM (R 4.1.0)
    ##  BiocNeighbors   * 1.12.0  2021-10-26 [1] Bioconductor
    ##  BiocParallel    * 1.28.3  2021-12-09 [1] Bioconductor
    ##  bit               4.0.4   2020-08-04 [1] RSPM (R 4.1.0)
    ##  bit64             4.0.5   2020-08-30 [1] RSPM (R 4.1.0)
    ##  cachem            1.0.6   2021-08-19 [1] RSPM (R 4.1.0)
    ##  callr             3.7.0   2021-04-20 [1] RSPM (R 4.1.0)
    ##  cli               3.1.0   2021-10-27 [1] RSPM (R 4.1.0)
    ##  cluster         * 2.1.2   2021-04-17 [2] CRAN (R 4.1.2)
    ##  codetools         0.2-18  2020-11-04 [2] CRAN (R 4.1.2)
    ##  colorspace        2.0-2   2021-06-24 [1] RSPM (R 4.1.0)
    ##  cowplot         * 1.1.1   2020-12-30 [1] RSPM (R 4.1.0)
    ##  crayon            1.4.2   2021-10-29 [1] RSPM (R 4.1.0)
    ##  data.table      * 1.14.2  2021-09-27 [1] RSPM (R 4.1.0)
    ##  deldir            1.0-6   2021-10-23 [1] RSPM (R 4.1.0)
    ##  deMULTIplex     * 1.0.2   2022-01-20 [1] Github (chris-mcginnis-ucsf/MULTI-seq@233a0c0)
    ##  desc              1.4.0   2021-09-28 [1] RSPM (R 4.1.0)
    ##  devtools        * 2.4.3   2021-11-30 [1] RSPM (R 4.1.0)
    ##  digest            0.6.29  2021-12-01 [1] RSPM (R 4.1.0)
    ##  dplyr           * 1.0.7   2021-06-18 [1] RSPM (R 4.1.0)
    ##  dsb             * 0.3.0   2022-01-05 [1] RSPM (R 4.1.2)
    ##  ellipsis          0.3.2   2021-04-29 [1] RSPM (R 4.1.0)
    ##  evaluate          0.14    2019-05-28 [1] RSPM (R 4.1.0)
    ##  fansi             1.0.2   2022-01-14 [1] RSPM (R 4.1.0)
    ##  farver            2.1.0   2021-02-28 [1] RSPM (R 4.1.0)
    ##  fastmap           1.1.0   2021-01-25 [1] RSPM (R 4.1.0)
    ##  fitdistrplus    * 1.1-6   2021-09-28 [1] RSPM (R 4.1.0)
    ##  fs                1.5.2   2021-12-08 [1] RSPM (R 4.1.0)
    ##  furrr           * 0.2.3   2021-06-25 [1] RSPM (R 4.1.0)
    ##  future          * 1.23.0  2021-10-31 [1] RSPM (R 4.1.0)
    ##  future.apply    * 1.8.1   2021-08-10 [1] RSPM (R 4.1.0)
    ##  generics          0.1.1   2021-10-25 [1] RSPM (R 4.1.0)
    ##  ggplot2         * 3.3.5   2021-06-25 [1] RSPM (R 4.1.0)
    ##  ggrepel           0.9.1   2021-01-15 [1] RSPM (R 4.1.0)
    ##  ggridges          0.5.3   2021-01-08 [1] RSPM (R 4.1.0)
    ##  globals           0.14.0  2020-11-22 [1] RSPM (R 4.1.0)
    ##  glue              1.6.0   2021-12-17 [1] RSPM (R 4.1.0)
    ##  goftest           1.2-3   2021-10-07 [1] RSPM (R 4.1.0)
    ##  gridExtra         2.3     2017-09-09 [1] RSPM (R 4.1.0)
    ##  gtable            0.3.0   2019-03-25 [1] RSPM (R 4.1.0)
    ##  here            * 1.0.1   2020-12-13 [1] RSPM (R 4.1.0)
    ##  hms               1.1.1   2021-09-26 [1] RSPM (R 4.1.0)
    ##  htmltools         0.5.2   2021-08-25 [1] RSPM (R 4.1.0)
    ##  htmlwidgets       1.5.4   2021-09-08 [1] RSPM (R 4.1.0)
    ##  httpuv            1.6.5   2022-01-05 [1] RSPM (R 4.1.0)
    ##  httr              1.4.2   2020-07-20 [1] RSPM (R 4.1.0)
    ##  ica               1.0-2   2018-05-24 [1] RSPM (R 4.1.0)
    ##  igraph            1.2.11  2022-01-04 [1] RSPM (R 4.1.0)
    ##  irlba             2.3.5   2021-12-06 [1] RSPM (R 4.1.0)
    ##  janitor         * 2.1.0   2021-01-05 [1] RSPM (R 4.1.0)
    ##  jsonlite          1.7.3   2022-01-17 [1] RSPM (R 4.1.0)
    ##  KernSmooth      * 2.23-20 2021-05-03 [2] CRAN (R 4.1.2)
    ##  knitr             1.37    2021-12-16 [1] RSPM (R 4.1.0)
    ##  labeling          0.4.2   2020-10-20 [1] RSPM (R 4.1.0)
    ##  later             1.3.0   2021-08-18 [1] RSPM (R 4.1.0)
    ##  lattice           0.20-45 2021-09-22 [2] CRAN (R 4.1.2)
    ##  lazyeval          0.2.2   2019-03-15 [1] RSPM (R 4.1.0)
    ##  leiden            0.3.9   2021-07-27 [1] RSPM (R 4.1.0)
    ##  lifecycle         1.0.1   2021-09-24 [1] RSPM (R 4.1.0)
    ##  limma             3.50.0  2021-10-26 [1] Bioconductor
    ##  listenv           0.8.0   2019-12-05 [1] RSPM (R 4.1.0)
    ##  lmtest            0.9-39  2021-11-07 [1] RSPM (R 4.1.0)
    ##  lubridate         1.8.0   2021-10-07 [1] RSPM (R 4.1.0)
    ##  magrittr        * 2.0.1   2020-11-17 [1] RSPM (R 4.1.0)
    ##  MASS            * 7.3-54  2021-05-03 [2] CRAN (R 4.1.2)
    ##  Matrix          * 1.3-4   2021-06-01 [2] CRAN (R 4.1.2)
    ##  matrixStats       0.61.0  2021-09-17 [1] RSPM (R 4.1.0)
    ##  mclust            5.4.9   2021-12-17 [1] RSPM (R 4.1.0)
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
    ##  readr           * 2.1.1   2021-11-30 [1] RSPM (R 4.1.0)
    ##  remotes           2.4.2   2021-11-30 [1] RSPM (R 4.1.0)
    ##  reshape2        * 1.4.4   2020-04-09 [1] RSPM (R 4.1.0)
    ##  reticulate        1.23    2022-01-14 [1] RSPM (R 4.1.0)
    ##  rlang           * 0.4.12  2021-10-18 [1] RSPM (R 4.1.0)
    ##  rmarkdown         2.11    2021-09-14 [1] RSPM (R 4.1.0)
    ##  ROCR              1.0-11  2020-05-02 [1] RSPM (R 4.1.0)
    ##  rpart             4.1-15  2019-04-12 [2] CRAN (R 4.1.2)
    ##  rprojroot         2.0.2   2020-11-15 [1] RSPM (R 4.1.0)
    ##  rstudioapi        0.13    2020-11-12 [1] RSPM (R 4.1.0)
    ##  Rtsne             0.15    2018-11-10 [1] RSPM (R 4.1.0)
    ##  S4Vectors         0.32.3  2021-11-21 [1] Bioconductor
    ##  scales            1.1.1   2020-05-11 [1] RSPM (R 4.1.0)
    ##  scattermore       0.7     2020-11-24 [1] RSPM (R 4.1.0)
    ##  sctransform       0.3.3   2022-01-13 [1] RSPM (R 4.1.0)
    ##  sessioninfo       1.2.2   2021-12-06 [1] RSPM (R 4.1.0)
    ##  Seurat          * 4.1.0   2022-01-14 [1] RSPM (R 4.1.0)
    ##  SeuratObject    * 4.0.4   2021-11-23 [1] RSPM (R 4.1.0)
    ##  shiny             1.7.1   2021-10-02 [1] RSPM (R 4.1.0)
    ##  snakecase         0.11.0  2019-05-25 [1] RSPM (R 4.1.0)
    ##  sp              * 1.4-6   2021-11-14 [1] RSPM (R 4.1.0)
    ##  spatstat.core     2.3-2   2021-11-26 [1] RSPM (R 4.1.0)
    ##  spatstat.data     2.1-2   2021-12-17 [1] RSPM (R 4.1.0)
    ##  spatstat.geom     2.3-1   2021-12-10 [1] RSPM (R 4.1.0)
    ##  spatstat.sparse   2.1-0   2021-12-17 [1] RSPM (R 4.1.0)
    ##  spatstat.utils    2.3-0   2021-12-12 [1] RSPM (R 4.1.0)
    ##  stringi           1.7.6   2021-11-29 [1] RSPM (R 4.1.0)
    ##  stringr         * 1.4.0   2019-02-10 [1] RSPM (R 4.1.0)
    ##  survival        * 3.2-13  2021-08-24 [2] CRAN (R 4.1.2)
    ##  tensor            1.5     2012-05-05 [1] RSPM (R 4.1.0)
    ##  testthat          3.1.1   2021-12-03 [1] RSPM (R 4.1.0)
    ##  tibble          * 3.1.6   2021-11-07 [1] RSPM (R 4.1.0)
    ##  tidyr           * 1.1.4   2021-09-27 [1] RSPM (R 4.1.0)
    ##  tidyselect        1.1.1   2021-04-30 [1] RSPM (R 4.1.0)
    ##  tzdb              0.2.0   2021-10-27 [1] RSPM (R 4.1.0)
    ##  usethis         * 2.1.5   2021-12-09 [1] RSPM (R 4.1.0)
    ##  utf8              1.2.2   2021-07-24 [1] RSPM (R 4.1.0)
    ##  uwot              0.1.11  2021-12-02 [1] RSPM (R 4.1.0)
    ##  vctrs             0.3.8   2021-04-29 [1] RSPM (R 4.1.0)
    ##  viridisLite       0.4.0   2021-04-13 [1] RSPM (R 4.1.0)
    ##  vroom             1.5.7   2021-11-30 [1] RSPM (R 4.1.0)
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
