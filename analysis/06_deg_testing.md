Analysis 06: Differential Expression Testing
================
Kevin Thomas
25 April, 2022

-   [Differential expression testing on each
    cluster](#differential-expression-testing-on-each-cluster)

## Differential expression testing on each cluster

Normally, this loop would be run on all clusters in the dataset. Across
all 25 clusters in the dataset, this analysis takes just under 3.5 hours
on a Google virtual machine with 32 cores and 128 Gb of RAM. For the
sake of demonstration, the following code runs only on the first cluster
of Naive B cells. On the same system, analysis of this first cluster
takes approximately 10 minutes.

``` r
# Convert to sce object
sce <- as.SingleCellExperiment(sc_total, assay = "RNA")

# Perform DE testing on each cluster
for (ct in levels(sce$clusters_annotated)[1]) {
## Remove this subset indicator "[1]".....^ to run DE analysis on all clusters
  message("Testing cell type: ", ct)
  # Set up object
  sc <- sce[,sce$clusters_annotated == ct]
  sc <- logNormCounts(sc)
  sca <- SceToSingleCellAssay(sc, check_sanity = TRUE)
  
  # Scaled cellular detection rate as a covariate
  cdr2 <- colSums(assay(sca) > 0)
  colData(sca)$ngeneson <- scale(cdr2)
  
  # Set up variables for comparisons
  colData(sca)$ancestry <- factor(colData(sca)$ancestry, levels = c("EA", "AA"))
  colData(sca)$classification <- factor(colData(sca)$classification, levels = c("Control", "SLE INACT", "SLE ACT"))
  colData(sca)$run <- factor(colData(sca)$run)
  colData(sca)$subject_id <- factor(colData(sca)$subject_id)
  colData(sca)$subject_id <- as.numeric(colData(sca)$subject_id)
  
  # Get names for file headers
  prefix <- sprintf("%02d", which(levels(sce$clusters_annotated) == ct))
  clustername <- str_replace_all(
    string = str_replace_all(ct, "[^[:alnum:]]", "_"),
    pattern = "_+",
    replacement = "_"
  )
  
  # Grid of all groups to be compared
  groups <- expand.grid(
    ancestry = str_sort(unique(sca$ancestry)),
    classification = levels(sca$classification)
  )
  # List of group indices for each comparison
  comparisons <- list(
    "INACTvsControl"     = list(c(3,4), c(1,2)),
    "ACTvsControl"       = list(c(5,6), c(1,2)),
    "ACTvsINACT"         = list(c(5,6), c(3,4)),
    "AA.CTRLvsEA.CTRL"   = list(c(1), c(2)),
    "AA.INACTvsEA.INACT" = list(c(3), c(4)),
    "AA.ACTvsEA.ACT"     = list(c(5), c(6)),
    "AA.INACTvsAA.CTRL"  = list(c(3), c(1)),
    "AA.ACTvsAA.CTRL"    = list(c(5), c(1)),
    "AA.ACTvsAA.INACT"   = list(c(5), c(3)),
    "EA.INACTvsEA.CTRL"  = list(c(4), c(2)),
    "EA.ACTvsEA.CTRL"    = list(c(6), c(2)),
    "EA.ACTvsEA.INACT"   = list(c(6), c(4))
  )
  # Do DE Testing for all comparisons
  res_list <- imap(
    .x = comparisons,
    .f = function(idxs, comp) {
      # idxs = comparisons[[4]]
      # comp = names(comparisons)[4]
      message(paste0("Comparison: ", comp))
      # Get the ancestry and classification groups necessary for the comparison
      anc <- unique(groups[c(idxs[[1]], idxs[[2]]), "ancestry"])
      class <- unique(groups[c(idxs[[1]], idxs[[2]]), "classification"])
      
      # Subset data to cells of interest
      sc_test <- sca[, sca$ancestry %in% anc & sca$classification %in% class]
      # Define group1 and group2
      sc_test$group <- ifelse(
        test = sc_test$ancestry %in% groups[idxs[[1]], "ancestry"] &
          sc_test$classification %in% groups[idxs[[1]], "classification"],
        yes = "02",
        no = "01"
      )
      
      # Filter on genes that are expressed in at least 5% of cells (or at least 10 cells, whichever is higher) in either group
      group1_thresh <- names(
        which(
          rowSums(assay(sc_test[,sc_test$group == "01"]) > 0) >= max(0.05*ncol(assay(sc_test[,sc_test$group == "01"])), 10)
        )
      )
      group2_thresh <- names(
        which(
          rowSums(assay(sc_test[,sc_test$group == "02"]) > 0) >= max(0.05*ncol(assay(sc_test[,sc_test$group == "02"])), 10)
        )
      )
      genes_thresh <- union(group1_thresh, group2_thresh)
      message(paste0("     Testing ", length(genes_thresh), " genes."))
      
      # Set up formula
      formula <- as.formula(paste0("~ group + ngeneson + run + (1 | subject_id)"))
      ## If multiple ancestries are in a comparison group, add ancestry as a blocking variable
      if (
        length(unique(sc_test[,sc_test$group == "01"]$ancestry)) > 1 ||
        length(unique(sc_test[,sc_test$group == "02"]$ancestry)) > 1
      ) {
        formula <- update(formula, ~ . + ancestry)
      }
      message(paste0("     Final formula: ", paste(as.character(formula), collapse = " ")))
      
      t_start <- Sys.time()
      # Fit zero-inflated model
      zlmCond <- zlm(
        formula = formula,
        sca = sc_test[genes_thresh,],
        method = "glmer",
        ebayes = FALSE,
        strictConvergence = FALSE,
        force = TRUE,
        fitArgsD = list(nAGQ=0),
        parallel = TRUE
      )
      
      # Tested contrast stats
      summaryCond <- summary(zlmCond, doLRT = "group02", parallel = TRUE, fitArgsD = list(nAGQ=0))
      t_end <- Sys.time()
      
      # Save necessary fit objects
      if (!dir.exists(paths = here(paste0("data/MAST_out/", prefix, "_", clustername, "_fit/")))) {
        dir.create(
          path = here(paste0("data/MAST_out/", prefix, "_", clustername, "_fit/")),
          recursive = TRUE
        )
      }
      saveRDS(
        object = list(zlmCond = zlmCond, summaryCond = summaryCond),
        file = here(
          paste0("data/MAST_out/", prefix, "_", clustername, "_fit/"),
          paste0(prefix, "_", clustername, "_", str_remove_all(comp, "[^[:alnum:]]"), ".RDS")
        )
      )
      
      # Get the output
      summaryDt <- summaryCond$datatable
      de.out <- merge(
        summaryDt[summaryDt$component == "H" & contrast == "group02", c(1, 3, 4)],
        summaryDt[summaryDt$component == "logFC" & contrast == "group02", c(1, 3, 7, 5, 6, 8)],
        by = c("primerid", "contrast")
      ) |>
        as_tibble() |>
        dplyr::rename(feature = primerid, p_val = `Pr(>Chisq)`, logFC = coef, logFC_hi = ci.hi, logFC_lo = ci.lo) |>
        select(-contrast) |>
        mutate(p_val_adj = p.adjust(p_val, method = "BH"))
      
      message(paste0("Test completed in ", round(as.numeric(t_end-t_start, units = "mins"), 3), " minutes"))
      return(de.out)
    }
  )
  # Write needed data to .csv for input to IPA
  df <- res_list |>
    map(select, feature, logFC, p_val, p_val_adj) |>
    map(dplyr::rename, P = p_val, FDR = p_val_adj) |>
    imap(
      .f = function(x, nm) {
        set_colnames(x, str_remove_all(c("feature", paste0(colnames(x)[-1], "_", prefix, "_", nm)), "^logFC_"))
      }
    ) |>
    purrr::reduce(full_join, by = "feature")
  write.csv(
    x = df,
    file = here(
      "data/MAST_out/",
      paste0(prefix, "_", clustername, "_IPA.csv")
    ),
    row.names = FALSE
  )
}
```

    ## Testing cell type: Naive B cells

    ## `fData` has no primerid.  I'll make something up.

    ## `cData` has no wellKey.  I'll make something up.

    ## Assuming data assay in position 2, with name logcounts is log-transformed.

    ## Comparison: INACTvsControl

    ##      Testing 674 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id) + ancestry

    ## Loading required namespace: lme4

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.749 minutes

    ## Comparison: ACTvsControl

    ##      Testing 774 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id) + ancestry

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.598 minutes

    ## Comparison: ACTvsINACT

    ##      Testing 798 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id) + ancestry

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.598 minutes

    ## Comparison: AA.CTRLvsEA.CTRL

    ##      Testing 626 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.363 minutes

    ## Comparison: AA.INACTvsEA.INACT

    ##      Testing 856 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs5b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.451 minutes

    ## Comparison: AA.ACTvsEA.ACT

    ##      Testing 1097 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients

    ## 
    ## Done!

    ## Test completed in 0.391 minutes

    ## Comparison: AA.INACTvsAA.CTRL

    ##      Testing 473 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs1b, runs2a, runs2b, runs5a, runs5b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.263 minutes

    ## Comparison: AA.ACTvsAA.CTRL

    ##      Testing 466 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs1b, runs2a, runs2b, runs5b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.234 minutes

    ## Comparison: AA.ACTvsAA.INACT

    ##      Testing 500 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs1b, runs2a, runs2b, runs5b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.253 minutes

    ## Comparison: EA.INACTvsEA.CTRL

    ##      Testing 869 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs3a, runs3b, runs4a, runs4b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.469 minutes

    ## Comparison: EA.ACTvsEA.CTRL

    ##      Testing 1083 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs3a, runs3b, runs4a, runs4b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.48 minutes

    ## Comparison: EA.ACTvsEA.INACT

    ##      Testing 1130 genes.

    ##      Final formula: ~ group + ngeneson + run + (1 | subject_id)

    ## Warning in .nextMethod(object = object, value = value): Coefficients runs3a, runs3b, runs4a, runs4b are never estimible and will be dropped.

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ## Test completed in 0.476 minutes

``` r
session_info()
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
    ##  date     2022-04-25                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  package              * version  date       lib source                                
    ##  abind                  1.4-5    2016-07-21 [1] RSPM (R 4.1.0)                        
    ##  assertthat             0.2.1    2019-03-21 [1] RSPM (R 4.1.0)                        
    ##  beachmat               2.8.0    2021-05-19 [1] Bioconductor                          
    ##  Biobase              * 2.52.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  BiocGenerics         * 0.38.0   2021-08-13 [1] bioc_git2r (@1db849a)                 
    ##  BiocManager          * 1.30.16  2021-06-15 [1] CRAN (R 4.1.0)                        
    ##  BiocNeighbors        * 1.10.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  BiocParallel         * 1.26.1   2021-07-04 [1] RSPM (R 4.1.0)                        
    ##  bitops                 1.0-7    2021-04-24 [1] RSPM (R 4.1.0)                        
    ##  boot                   1.3-28   2021-05-03 [2] CRAN (R 4.1.0)                        
    ##  cachem                 1.0.5    2021-05-15 [1] RSPM (R 4.1.0)                        
    ##  callr                  3.7.0    2021-04-20 [1] RSPM (R 4.1.0)                        
    ##  cli                    3.0.1    2021-07-17 [1] RSPM (R 4.1.0)                        
    ##  cluster                2.1.2    2021-04-17 [2] CRAN (R 4.1.0)                        
    ##  codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)                        
    ##  colorspace             2.0-2    2021-06-24 [1] RSPM (R 4.1.0)                        
    ##  cowplot              * 1.1.1    2020-12-30 [1] RSPM (R 4.1.0)                        
    ##  crayon                 1.4.1    2021-02-08 [1] RSPM (R 4.1.0)                        
    ##  data.table           * 1.14.0   2021-02-21 [1] RSPM (R 4.1.0)                        
    ##  DBI                    1.1.1    2021-01-15 [1] RSPM (R 4.1.0)                        
    ##  DelayedArray           0.18.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  DelayedMatrixStats     1.14.2   2021-08-08 [1] Bioconductor                          
    ##  deldir                 1.0-6    2021-10-23 [1] RSPM (R 4.1.0)                        
    ##  desc                   1.3.0    2021-03-05 [1] RSPM (R 4.1.0)                        
    ##  devtools             * 2.4.2    2021-06-07 [1] RSPM (R 4.1.0)                        
    ##  digest                 0.6.27   2020-10-24 [1] RSPM (R 4.1.0)                        
    ##  dplyr                * 1.0.7    2021-06-18 [1] RSPM (R 4.1.0)                        
    ##  ellipsis               0.3.2    2021-04-29 [1] RSPM (R 4.1.0)                        
    ##  evaluate               0.14     2019-05-28 [1] RSPM (R 4.1.0)                        
    ##  fansi                  0.5.0    2021-05-25 [1] RSPM (R 4.1.0)                        
    ##  fastmap                1.1.0    2021-01-25 [1] RSPM (R 4.1.0)                        
    ##  fitdistrplus           1.1-5    2021-05-28 [1] RSPM (R 4.1.0)                        
    ##  fs                     1.5.0    2020-07-31 [1] RSPM (R 4.1.0)                        
    ##  future               * 1.21.0   2020-12-10 [1] RSPM (R 4.1.0)                        
    ##  future.apply           1.7.0    2021-01-04 [1] RSPM (R 4.1.0)                        
    ##  generics               0.1.0    2020-10-31 [1] RSPM (R 4.1.0)                        
    ##  GenomeInfoDb         * 1.28.1   2021-07-01 [1] RSPM (R 4.1.0)                        
    ##  GenomeInfoDbData       1.2.6    2021-08-04 [1] RSPM (R 4.1.0)                        
    ##  GenomicRanges        * 1.44.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  ggplot2              * 3.3.5    2021-06-25 [1] RSPM (R 4.1.0)                        
    ##  ggrepel              * 0.9.1    2021-01-15 [1] RSPM (R 4.1.0)                        
    ##  ggridges               0.5.3    2021-01-08 [1] RSPM (R 4.1.0)                        
    ##  globals                0.14.0   2020-11-22 [1] RSPM (R 4.1.0)                        
    ##  glue                   1.4.2    2020-08-27 [1] RSPM (R 4.1.0)                        
    ##  goftest                1.2-2    2019-12-02 [1] RSPM (R 4.1.0)                        
    ##  gridExtra              2.3      2017-09-09 [1] RSPM (R 4.1.0)                        
    ##  gtable                 0.3.0    2019-03-25 [1] RSPM (R 4.1.0)                        
    ##  here                 * 1.0.1    2020-12-13 [1] RSPM (R 4.1.0)                        
    ##  hms                    1.1.0    2021-05-17 [1] RSPM (R 4.1.0)                        
    ##  htmltools              0.5.1.1  2021-01-22 [1] RSPM (R 4.1.0)                        
    ##  htmlwidgets            1.5.3    2020-12-10 [1] RSPM (R 4.1.0)                        
    ##  httpuv                 1.6.1    2021-05-07 [1] RSPM (R 4.1.0)                        
    ##  httr                   1.4.2    2020-07-20 [1] RSPM (R 4.1.0)                        
    ##  ica                    1.0-2    2018-05-24 [1] RSPM (R 4.1.0)                        
    ##  igraph                 1.2.6    2020-10-06 [1] RSPM (R 4.1.0)                        
    ##  IRanges              * 2.26.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  irlba                  2.3.3    2019-02-05 [1] RSPM (R 4.1.0)                        
    ##  jsonlite               1.7.2    2020-12-09 [1] RSPM (R 4.1.0)                        
    ##  KernSmooth             2.23-20  2021-05-03 [2] CRAN (R 4.1.0)                        
    ##  knitr                  1.33     2021-04-24 [1] RSPM (R 4.1.0)                        
    ##  later                  1.2.0    2021-04-23 [1] RSPM (R 4.1.0)                        
    ##  lattice                0.20-44  2021-05-02 [2] CRAN (R 4.1.0)                        
    ##  lazyeval               0.2.2    2019-03-15 [1] RSPM (R 4.1.0)                        
    ##  leiden                 0.3.9    2021-07-27 [1] RSPM (R 4.1.0)                        
    ##  lifecycle              1.0.0    2021-02-15 [1] RSPM (R 4.1.0)                        
    ##  listenv                0.8.0    2019-12-05 [1] RSPM (R 4.1.0)                        
    ##  lme4                   1.1-27.1 2021-06-22 [1] RSPM (R 4.1.0)                        
    ##  lmtest                 0.9-38   2020-09-09 [1] RSPM (R 4.1.0)                        
    ##  magrittr             * 2.0.1    2020-11-17 [1] RSPM (R 4.1.0)                        
    ##  MASS                   7.3-54   2021-05-03 [2] CRAN (R 4.1.0)                        
    ##  MAST                 * 1.18.0   2021-05-19 [1] Bioconductor                          
    ##  Matrix                 1.3-4    2021-06-01 [2] RSPM (R 4.1.0)                        
    ##  MatrixGenerics       * 1.4.2    2021-08-08 [1] Bioconductor                          
    ##  matrixStats          * 0.60.0   2021-07-26 [1] RSPM (R 4.1.0)                        
    ##  memoise                2.0.0    2021-01-26 [1] RSPM (R 4.1.0)                        
    ##  mgcv                   1.8-36   2021-06-01 [2] RSPM (R 4.1.0)                        
    ##  mime                   0.11     2021-06-23 [1] RSPM (R 4.1.0)                        
    ##  miniUI                 0.1.1.1  2018-05-18 [1] RSPM (R 4.1.0)                        
    ##  minqa                  1.2.4    2014-10-09 [1] RSPM (R 4.1.0)                        
    ##  munsell                0.5.0    2018-06-12 [1] RSPM (R 4.1.0)                        
    ##  nlme                   3.1-152  2021-02-04 [2] CRAN (R 4.1.0)                        
    ##  nloptr                 1.2.2.2  2020-07-02 [1] RSPM (R 4.1.0)                        
    ##  parallelly             1.27.0   2021-07-19 [1] RSPM (R 4.1.0)                        
    ##  patchwork              1.1.1    2020-12-17 [1] RSPM (R 4.1.0)                        
    ##  pbapply                1.4-3    2020-08-18 [1] RSPM (R 4.1.0)                        
    ##  pillar                 1.6.2    2021-07-29 [1] RSPM (R 4.1.0)                        
    ##  pkgbuild               1.2.0    2020-12-15 [1] RSPM (R 4.1.0)                        
    ##  pkgconfig              2.0.3    2019-09-22 [1] RSPM (R 4.1.0)                        
    ##  pkgload                1.2.1    2021-04-06 [1] RSPM (R 4.1.0)                        
    ##  plotly                 4.9.4.1  2021-06-18 [1] RSPM (R 4.1.0)                        
    ##  plyr                   1.8.6    2020-03-03 [1] RSPM (R 4.1.0)                        
    ##  png                    0.1-7    2013-12-03 [1] RSPM (R 4.1.0)                        
    ##  polyclip               1.10-0   2019-03-14 [1] RSPM (R 4.1.0)                        
    ##  presto               * 1.0.0    2021-08-13 [1] Github (immunogenomics/presto@052085d)
    ##  prettyunits            1.1.1    2020-01-24 [1] RSPM (R 4.1.0)                        
    ##  processx               3.5.2    2021-04-30 [1] RSPM (R 4.1.0)                        
    ##  progress               1.2.2    2019-05-16 [1] RSPM (R 4.1.0)                        
    ##  promises               1.2.0.1  2021-02-11 [1] RSPM (R 4.1.0)                        
    ##  ps                     1.6.0    2021-02-28 [1] RSPM (R 4.1.0)                        
    ##  purrr                * 0.3.4    2020-04-17 [1] RSPM (R 4.1.0)                        
    ##  R6                     2.5.0    2020-10-28 [1] RSPM (R 4.1.0)                        
    ##  RANN                   2.6.1    2019-01-08 [1] RSPM (R 4.1.0)                        
    ##  RColorBrewer           1.1-2    2014-12-07 [1] RSPM (R 4.1.0)                        
    ##  Rcpp                 * 1.0.7    2021-07-07 [1] RSPM (R 4.1.0)                        
    ##  RcppAnnoy              0.0.19   2021-07-30 [1] RSPM (R 4.1.0)                        
    ##  RCurl                  1.98-1.3 2021-03-16 [1] RSPM (R 4.1.0)                        
    ##  remotes                2.4.0    2021-06-02 [1] RSPM (R 4.1.0)                        
    ##  reshape2               1.4.4    2020-04-09 [1] RSPM (R 4.1.0)                        
    ##  reticulate             1.20     2021-05-03 [1] RSPM (R 4.1.0)                        
    ##  rlang                  0.4.11   2021-04-30 [1] RSPM (R 4.1.0)                        
    ##  rmarkdown              2.10     2021-08-06 [1] RSPM (R 4.1.0)                        
    ##  ROCR                   1.0-11   2020-05-02 [1] RSPM (R 4.1.0)                        
    ##  rpart                  4.1-15   2019-04-12 [2] CRAN (R 4.1.0)                        
    ##  rprojroot              2.0.2    2020-11-15 [1] RSPM (R 4.1.0)                        
    ##  rstudioapi             0.13     2020-11-12 [1] RSPM (R 4.1.0)                        
    ##  Rtsne                  0.15     2018-11-10 [1] RSPM (R 4.1.0)                        
    ##  S4Vectors            * 0.30.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  scales                 1.1.1    2020-05-11 [1] RSPM (R 4.1.0)                        
    ##  scattermore            0.7      2020-11-24 [1] RSPM (R 4.1.0)                        
    ##  sctransform            0.3.2    2020-12-16 [1] RSPM (R 4.1.0)                        
    ##  scuttle              * 1.2.1    2021-08-05 [1] Bioconductor                          
    ##  sessioninfo            1.1.1    2018-11-05 [1] RSPM (R 4.1.0)                        
    ##  Seurat               * 4.0.3    2021-08-13 [1] Github (satijalab/Seurat@9b38929)     
    ##  SeuratObject         * 4.0.2    2021-06-09 [1] RSPM (R 4.1.0)                        
    ##  shiny                  1.6.0    2021-01-25 [1] RSPM (R 4.1.0)                        
    ##  SingleCellExperiment * 1.14.1   2021-08-13 [1] bioc_git2r (@5357eff)                 
    ##  sparseMatrixStats      1.4.2    2021-08-08 [1] Bioconductor                          
    ##  spatstat.core          2.3-0    2021-07-16 [1] RSPM (R 4.1.0)                        
    ##  spatstat.data          2.1-0    2021-03-21 [1] RSPM (R 4.1.0)                        
    ##  spatstat.geom          2.4-0    2022-03-29 [1] RSPM (R 4.1.0)                        
    ##  spatstat.sparse        2.0-0    2021-03-16 [1] RSPM (R 4.1.0)                        
    ##  spatstat.utils         2.2-0    2021-06-14 [1] RSPM (R 4.1.0)                        
    ##  stringi                1.7.3    2021-07-16 [1] RSPM (R 4.1.0)                        
    ##  stringr              * 1.4.0    2019-02-10 [1] RSPM (R 4.1.0)                        
    ##  SummarizedExperiment * 1.22.0   2021-08-13 [1] bioc_git2r (@7d1110e)                 
    ##  survival               3.2-11   2021-04-26 [2] CRAN (R 4.1.0)                        
    ##  tensor                 1.5      2012-05-05 [1] RSPM (R 4.1.0)                        
    ##  testthat               3.0.4    2021-07-01 [1] RSPM (R 4.1.0)                        
    ##  tibble               * 3.1.3    2021-07-23 [1] RSPM (R 4.1.0)                        
    ##  tidyr                  1.1.3    2021-03-03 [1] RSPM (R 4.1.0)                        
    ##  tidyselect             1.1.1    2021-04-30 [1] RSPM (R 4.1.0)                        
    ##  usethis              * 2.0.1    2021-02-10 [1] RSPM (R 4.1.0)                        
    ##  utf8                   1.2.2    2021-07-24 [1] RSPM (R 4.1.0)                        
    ##  uwot                   0.1.10   2020-12-15 [1] RSPM (R 4.1.0)                        
    ##  vctrs                  0.3.8    2021-04-29 [1] RSPM (R 4.1.0)                        
    ##  viridisLite            0.4.0    2021-04-13 [1] RSPM (R 4.1.0)                        
    ##  withr                  2.4.2    2021-04-18 [1] RSPM (R 4.1.0)                        
    ##  xfun                   0.25     2021-08-06 [1] RSPM (R 4.1.0)                        
    ##  xtable                 1.8-4    2019-04-21 [1] RSPM (R 4.1.0)                        
    ##  XVector                0.32.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  yaml                   2.2.1    2020-02-01 [1] RSPM (R 4.1.0)                        
    ##  zlibbioc               1.38.0   2021-05-19 [1] RSPM (R 4.1.0)                        
    ##  zoo                    1.8-9    2021-03-09 [1] RSPM (R 4.1.0)                        
    ## 
    ## [1] /usr/local/lib/R/site-library
    ## [2] /usr/local/lib/R/library
