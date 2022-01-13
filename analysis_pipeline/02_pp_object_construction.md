Preprocessing 02: Multimodal object construction
================
Kevin Thomas
17 November, 2021

-   [Construct multimodal object with filtered counts from previous
    step](#construct-multimodal-object-with-filtered-counts-from-previous-step)

## Construct multimodal object with filtered counts from previous step

``` r
# Load data from last step and samplesheet for object metadata
# Here, gene expression matrices are simply taken as the adjusted counts output of the SoupX analysis
gem_matrices <- readRDS(file = "../data/adjustCountList.RDS")
samplesheet <- read_csv(file = "../data/da_samplesheet_final.csv") %>%
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
  .x = paste0(dataDirs, "/filtered_feature_bc_matrix"),
  .f = Read10X
)
```

    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.
    ## 10X data contains more than one type and is being returned as a list containing matrices of each type.

``` r
# Build adt_matrices if filtered_matrices contain "Antibody Capture" data and the sample is included in gem_matrices
adt_matrices <- ifelse(
  test = mapply(
    FUN = function(x,y){"Antibody Capture" %in% names(x) && y %in% names(gem_matrices)},
    x = filtered_matrices,
    y = basename(dataDirs),
    SIMPLIFY = FALSE
  ),
  yes = filtered_matrices %>% map("Antibody Capture"),
  no = NULL
) %>%
  purrr::discard(.p = is.null) %>%
  setNames(names(gem_matrices))

# Free memory
rm(filtered_matrices)
gc()
```

    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   9063148  484.1   15075078  805.1  10055418  537.1
    ## Vcells 196009794 1495.5  568406240 4336.6 565644408 4315.6

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
        adt <- DSBNormalizeProtein(
          cell_protein_matrix = adt_matrices[[i]],
          empty_drop_matrix = adt_matrices[[i]][ , gem_classification %>%
                                                   dplyr::filter(run == names(adt_matrices)[[i]]) %>%
                                                   dplyr::filter(call == str_remove_all(secondary_species_prefix, "_")) %>%
                                                   pull(cells)], 
          use.isotype.control = FALSE,
          isotype.control.name.vec = NULL
        )
        return(adt)
      }else{
        message(str_glue("Not normalizing sample {names(adt_matrices)[[i]]}. CITE data insufficient."))
        adt <- adt_matrices[[i]]
        return(adt)
      }
    }
  ) %>%
    setNames(names(adt_matrices)) %>%
    purrr::discard(.p = is.null)
}else{
  NULL
}
```

    ## [1] "potential isotype controls detected"
    ## [1] "isotype-control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s1a of 10.

    ## Keeping sample s1a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype-control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s1b of 10.

    ## Keeping sample s1b.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype-control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s2a of 10.

    ## Keeping sample s2a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype-control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s2b of 10.

    ## Keeping sample s2b.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s3a of 10.

    ## Keeping sample s3a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s3b of 10.

    ## Keeping sample s3b.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s4a of 10.

    ## Keeping sample s4a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s4b of 10.

    ## Keeping sample s4b.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s5a of 10.

    ## Keeping sample s5a.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

    ## [1] "potential isotype controls detected"
    ## [1] "isotype_control"
    ## [1] "correcting ambient protein background noise"
    ## [1] "calculating dsb technical component for each cell to remove cell to cell techncial noise"

    ## Evaluating CITE counts for sample s5b of 10.

    ## Keeping sample s5b.

    ## Warning in DSBNormalizeProtein(cell_protein_matrix = adt_matrices[[i]], : denoise.counts = TRUE with use.isotype.control = FALSE is not recommended if isotype controls are available. If the raw data include isotype controls, set `denoise.counts` = TRUE `use.isotype.control` = TRUE and set `isotype.control.name.vec` to a vector of isotype control protien names from cell_protein_matrix

``` r
# Construct each run as a multimodal seurat object, stored in a list
## DO NOT FUTURE_MAP THIS STEP (parallelization is already happening in the `constructGEMs` function)
gems <- map(
  .x = seq_along(gem_matrices[names(dsb_adt_matrices)]),
  .f = function(i) {
    message(glue("Constructing seurat object for run {names(gem_matrices[names(dsb_adt_matrices)])[[i]]}"))
    x <- constructGEMs(
      sample = names(gem_matrices[names(dsb_adt_matrices)])[[i]],
      gem_matrix = gem_matrices[names(dsb_adt_matrices)][[i]],
      adt_matrix = adt_matrices[names(dsb_adt_matrices)][[i]],
      norm_adt_matrix = dsb_adt_matrices[[i]],
      hashtags = if (multisample == "hash") {levels(factor(samplesheet$hashtag))}else{FALSE}
    )
    return(x)
  }
) %>%
  `names<-`(names(gem_matrices[names(dsb_adt_matrices)])) %>%
  purrr::discard(.p = is.null)
```

    ## Constructing seurat object for run s1a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.2961512
    ## 
    ## $start.arg$mu
    ## [1] 1.419929
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.07631566
    ## 
    ## $start.arg$mu
    ## [1] 9.259357
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.08655052
    ## 
    ## $start.arg$mu
    ## [1] 24.82159
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.3397566
    ## 
    ## $start.arg$mu
    ## [1] 37.63203
    ## 
    ## 
    ## $fix.arg
    ## NULL

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

    ## Constructing seurat object for run s1b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1260973
    ## 
    ## $start.arg$mu
    ## [1] 6.968729
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.8028103
    ## 
    ## $start.arg$mu
    ## [1] 1.849123
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.7735155
    ## 
    ## $start.arg$mu
    ## [1] 35.46082
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.08945614
    ## 
    ## $start.arg$mu
    ## [1] 131.1311
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 12270

    ## Cutoff for HT-1: 54 at quantile = 0.995

    ## Cutoff for HT-2: 12 at quantile = 0.999

    ## Cutoff for HT-3: 82 at quantile = 0.9

    ## Cutoff for HT-4: 214 at quantile = 0.8

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 991 Negative cells removed

    ## Round 2: 51 Negative cells removed

    ## Round 3: 7 Negative cells removed

    ## Round 4: Done! All cells classified.

    ## HT-1 cells: 2951

    ## HT-2 cells: 420

    ## HT-3 cells: 3100

    ## HT-4 cells: 5745

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s2a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.4142544
    ## 
    ## $start.arg$mu
    ## [1] 9.271155
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.126265
    ## 
    ## $start.arg$mu
    ## [1] 39.74274
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.088517
    ## 
    ## $start.arg$mu
    ## [1] 73.79957
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5278577
    ## 
    ## $start.arg$mu
    ## [1] 39.27332
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11691

    ## Cutoff for HT-1: 68 at quantile = 0.999

    ## Cutoff for HT-2: 131 at quantile = 0.995

    ## Cutoff for HT-3: 161 at quantile = 0.95

    ## Cutoff for HT-4: 239 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 1: 3578 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 2: 835 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 3: 709 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 4: 894 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 5: 978 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 6: 863 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 7: 142 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 8: 8 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 9: Done! All cells classified.

    ## HT-1 cells: 3317

    ## HT-2 cells: 1563

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."

    ## HT-3 cells rescued: 585

    ## HT-4 cells: 4369

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s2b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.7073275
    ## 
    ## $start.arg$mu
    ## [1] 4.873098
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.077578
    ## 
    ## $start.arg$mu
    ## [1] 23.30568
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 2.402003
    ## 
    ## $start.arg$mu
    ## [1] 71.99864
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.746445
    ## 
    ## $start.arg$mu
    ## [1] 36.10997
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11139

    ## Cutoff for HT-1: 25 at quantile = 0.995

    ## Cutoff for HT-2: 79 at quantile = 0.99

    ## Cutoff for HT-3: 134 at quantile = 0.9

    ## Cutoff for HT-4: 192 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 2264 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 2: 2356 Negative cells removed

    ## Round 3: 1132 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 4: 657 Negative cells removed

    ## Round 5: 55 Negative cells removed

    ## Round 6: Done! All cells classified.

    ## HT-1 cells: 2659

    ## HT-2 cells: 3647

    ## HT-3 cells: 1590

    ## HT-4 cells: 673

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Constructing seurat object for run s3a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1106102
    ## 
    ## $start.arg$mu
    ## [1] 0.04272134
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1281504
    ## 
    ## $start.arg$mu
    ## [1] 0.3605543
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1492889
    ## 
    ## $start.arg$mu
    ## [1] 0.6325736
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1453724
    ## 
    ## $start.arg$mu
    ## [1] 0.4694766
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1360989
    ## 
    ## $start.arg$mu
    ## [1] 1.234796
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 7584

    ## Cutoff for HT-1: 2 at quantile = 0.999

    ## Cutoff for HT-2: 10 at quantile = 0.999

    ## Cutoff for HT-3: 15 at quantile = 0.999

    ## Cutoff for HT-4: 11 at quantile = 0.999

    ## Cutoff for HT-5: 22 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## Round 1: 14812 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 2: 17065 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 3: 7463 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 4: 16376 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 5: 4431 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 6: 3917 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 7: 3848 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 8: 9959 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 9: 2272 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 10: 2536 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 11: 3303 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 12: 2407 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 13: 1641 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 14: 1598 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 15: 697 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## Round 16: 1176 Negative cells removed

    ## Round 17: Stopping, as no threshold extrema could be found.

    ## HT-1 cells: 998

    ## HT-2 cells: 620

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## HT-3 cells rescued: 3268

    ## HT-4 cells: 740

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## HT-5 cells rescued: 5009

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning in max(.): no non-missing arguments to max; returning -Inf

    ## Warning: Removed 46 rows containing missing values (geom_point).

    ## Constructing seurat object for run s3b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 1.941319
    ## 
    ## $start.arg$mu
    ## [1] 30.60114
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5576581
    ## 
    ## $start.arg$mu
    ## [1] 3.490444
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.021632
    ## 
    ## $start.arg$mu
    ## [1] 13.25639
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5876622
    ## 
    ## $start.arg$mu
    ## [1] 13.28543
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.624614
    ## 
    ## $start.arg$mu
    ## [1] 28.8136
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 15258

    ## Cutoff for HT-1: 48 at quantile = 0.8

    ## Cutoff for HT-2: 21 at quantile = 0.995

    ## Cutoff for HT-3: 76 at quantile = 0.999

    ## Cutoff for HT-4: 91 at quantile = 0.999

    ## Cutoff for HT-5: 190 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 2562 Negative cells removed

    ## Round 2: 206 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 3178 Negative cells removed

    ## Round 4: 517 Negative cells removed

    ## Round 5: 159 Negative cells removed

    ## Round 6: Done! All cells classified.

    ## HT-1 cells: 3290

    ## HT-2 cells: 2090

    ## HT-3 cells: 2334

    ## HT-4 cells: 1868

    ## HT-5 cells: 2957

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s4a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 12.51834
    ## 
    ## $start.arg$mu
    ## [1] 2.65494
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 24.02663
    ## 
    ## $start.arg$mu
    ## [1] 5.294374
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 5.772354
    ## 
    ## $start.arg$mu
    ## [1] 4.549455
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.188791
    ## 
    ## $start.arg$mu
    ## [1] 17.75699
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.795221
    ## 
    ## $start.arg$mu
    ## [1] 95.14064
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 27263

    ## Cutoff for HT-1: 4 at quantile = 0.8

    ## Cutoff for HT-2: 7 at quantile = 0.8

    ## Cutoff for HT-3: 14 at quantile = 0.999

    ## Cutoff for HT-4: 52 at quantile = 0.999

    ## Cutoff for HT-5: 245 at quantile = 0.99

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 23961 Negative cells removed

    ## Round 2: 5210 Negative cells removed

    ## Round 3: Done! All cells classified.

    ## HT-1 cells: 3418

    ## HT-2 cells: 5661

    ## HT-3 cells: 1176

    ## HT-4 cells: 4040

    ## HT-5 cells: 8264

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: did not converge in 10 iterations

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning in max(.): no non-missing arguments to max; returning -Inf

    ## Constructing seurat object for run s4b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 3.018579
    ## 
    ## $start.arg$mu
    ## [1] 11.0267
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 2.984052
    ## 
    ## $start.arg$mu
    ## [1] 13.99061
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.529136
    ## 
    ## $start.arg$mu
    ## [1] 9.3195
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.839596
    ## 
    ## $start.arg$mu
    ## [1] 25.19017
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.686838
    ## 
    ## $start.arg$mu
    ## [1] 120.8595
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 21673

    ## Cutoff for HT-1: 16 at quantile = 0.8

    ## Cutoff for HT-2: 20 at quantile = 0.8

    ## Cutoff for HT-3: 29 at quantile = 0.999

    ## Cutoff for HT-4: 71 at quantile = 0.995

    ## Cutoff for HT-5: 330 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-1..."

    ## Round 1: 18251 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 2: 5395 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 3896 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 4: 4295 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 5: 3325 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 6: 2799 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 7: 2392 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 8: 2247 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 9: 1841 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 10: 1920 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 11: 1426 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 12: 758 Negative cells removed

    ## Round 13: 27 Negative cells removed

    ## Round 14: 2 Negative cells removed

    ## Round 15: 1 Negative cells removed

    ## Round 16: Done! All cells classified.

    ## HT-1 cells: 1047

    ## HT-2 cells: 863

    ## HT-3 cells: 705

    ## HT-4 cells: 827

    ## HT-5 cells: 870

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s5a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 4.044071
    ## 
    ## $start.arg$mu
    ## [1] 96.01433
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.194757
    ## 
    ## $start.arg$mu
    ## [1] 5.450637
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.553368
    ## 
    ## $start.arg$mu
    ## [1] 9.690797
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.9665492
    ## 
    ## $start.arg$mu
    ## [1] 34.00716
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.18407
    ## 
    ## $start.arg$mu
    ## [1] 80.78738
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 14900

    ## Cutoff for HT-1: 138 at quantile = 0.8

    ## Cutoff for HT-2: 30 at quantile = 0.999

    ## Cutoff for HT-3: 42 at quantile = 0.999

    ## Cutoff for HT-4: 176 at quantile = 0.999

    ## Cutoff for HT-5: 434 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-1..."

    ## Round 1: 7411 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 2: 7518 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 120 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 4: Done! All cells classified.

    ## [1] "No threshold found for HT-1..."

    ## HT-1 cells rescued: 8461

    ## HT-2 cells: 1889

    ## HT-3 cells: 2046

    ## HT-4 cells: 2372

    ## HT-5 cells: 3437

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Removed 25 rows containing missing values (geom_point).

    ## Constructing seurat object for run s5b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.880844
    ## 
    ## $start.arg$mu
    ## [1] 3.671127
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.8538241
    ## 
    ## $start.arg$mu
    ## [1] 4.350469
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.326974
    ## 
    ## $start.arg$mu
    ## [1] 10.37911
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.6993285
    ## 
    ## $start.arg$mu
    ## [1] 12.89977
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.779457
    ## 
    ## $start.arg$mu
    ## [1] 198.8209
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11225

    ## Cutoff for HT-1: 23 at quantile = 0.999

    ## Cutoff for HT-2: 28 at quantile = 0.999

    ## Cutoff for HT-3: 39 at quantile = 0.99

    ## Cutoff for HT-4: 75 at quantile = 0.999

    ## Cutoff for HT-5: 1036 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 4444 Negative cells removed

    ## Round 2: 98 Negative cells removed

    ## Round 3: 25 Negative cells removed

    ## Round 4: 2 Negative cells removed

    ## Round 5: Done! All cells classified.

    ## HT-1 cells: 2070

    ## HT-2 cells: 2124

    ## HT-3 cells: 2270

    ## HT-4 cells: 2077

    ## HT-5 cells: 2610

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

``` r
gems <- map(
  .x = seq_along(gem_matrices[names(dsb_adt_matrices)]),
  .f = function(i) {
    message(glue("Constructing seurat object for run {names(gem_matrices[names(dsb_adt_matrices)])[[i]]}"))
    x <- constructGEMs(
      sample = names(gem_matrices[names(dsb_adt_matrices)])[[i]],
      gem_matrix = gem_matrices[names(dsb_adt_matrices)][[i]],
      adt_matrix = adt_matrices[names(dsb_adt_matrices)][[i]],
      norm_adt_matrix = dsb_adt_matrices[[i]],
      hashtags = if (multisample == "hash") {levels(factor(samplesheet$hashtag))}else{FALSE}
    )
    return(x)
  }
) %>%
  `names<-`(names(gem_matrices[names(dsb_adt_matrices)])) %>%
  purrr::discard(.p = is.null)
```

    ## Constructing seurat object for run s1a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.2961512
    ## 
    ## $start.arg$mu
    ## [1] 1.419929
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.07631566
    ## 
    ## $start.arg$mu
    ## [1] 9.259357
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.08655052
    ## 
    ## $start.arg$mu
    ## [1] 24.82159
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.3397566
    ## 
    ## $start.arg$mu
    ## [1] 37.63203
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 13180

    ## Cutoff for HT-1: 16 at quantile = 0.999

    ## Cutoff for HT-2: 25 at quantile = 0.9

    ## Cutoff for HT-3: 50 at quantile = 0.85

    ## Cutoff for HT-4: 202 at quantile = 0.99

    ## Demuxing hashtags with deMULTIplex to rescue singlets

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

    ## Constructing seurat object for run s1b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1260973
    ## 
    ## $start.arg$mu
    ## [1] 6.968729
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.8028103
    ## 
    ## $start.arg$mu
    ## [1] 1.849123
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.7735155
    ## 
    ## $start.arg$mu
    ## [1] 35.46082
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.08945614
    ## 
    ## $start.arg$mu
    ## [1] 131.1311
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 12270

    ## Cutoff for HT-1: 54 at quantile = 0.995

    ## Cutoff for HT-2: 12 at quantile = 0.999

    ## Cutoff for HT-3: 82 at quantile = 0.9

    ## Cutoff for HT-4: 214 at quantile = 0.8

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 991 Negative cells removed

    ## Round 2: 51 Negative cells removed

    ## Round 3: 7 Negative cells removed

    ## Round 4: Done! All cells classified.

    ## HT-1 cells: 2951

    ## HT-2 cells: 420

    ## HT-3 cells: 3100

    ## HT-4 cells: 5745

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s2a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.4142544
    ## 
    ## $start.arg$mu
    ## [1] 9.271155
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.126265
    ## 
    ## $start.arg$mu
    ## [1] 39.74274
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.088517
    ## 
    ## $start.arg$mu
    ## [1] 73.79957
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5278577
    ## 
    ## $start.arg$mu
    ## [1] 39.27332
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11691

    ## Cutoff for HT-1: 68 at quantile = 0.999

    ## Cutoff for HT-2: 131 at quantile = 0.995

    ## Cutoff for HT-3: 161 at quantile = 0.95

    ## Cutoff for HT-4: 239 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 1: 3578 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 2: 835 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 3: 709 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 4: 894 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 5: 978 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."

    ## Round 6: 863 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 7: 142 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 8: 8 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 9: Done! All cells classified.

    ## HT-1 cells: 3317

    ## HT-2 cells: 1563

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."

    ## HT-3 cells rescued: 585

    ## HT-4 cells: 4369

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s2b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.7073275
    ## 
    ## $start.arg$mu
    ## [1] 4.873098
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.077578
    ## 
    ## $start.arg$mu
    ## [1] 23.30568
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 2.402003
    ## 
    ## $start.arg$mu
    ## [1] 71.99864
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.746445
    ## 
    ## $start.arg$mu
    ## [1] 36.10997
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11139

    ## Cutoff for HT-1: 25 at quantile = 0.995

    ## Cutoff for HT-2: 79 at quantile = 0.99

    ## Cutoff for HT-3: 134 at quantile = 0.9

    ## Cutoff for HT-4: 192 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 2264 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 2: 2356 Negative cells removed

    ## Round 3: 1132 Negative cells removed

    ## [1] "No threshold found for HT-3..."

    ## Round 4: 657 Negative cells removed

    ## Round 5: 55 Negative cells removed

    ## Round 6: Done! All cells classified.

    ## HT-1 cells: 2659

    ## HT-2 cells: 3647

    ## HT-3 cells: 1590

    ## HT-4 cells: 673

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Removed 18 rows containing missing values (geom_point).

    ## Constructing seurat object for run s3a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1106102
    ## 
    ## $start.arg$mu
    ## [1] 0.04272134
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1281504
    ## 
    ## $start.arg$mu
    ## [1] 0.3605543
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1492889
    ## 
    ## $start.arg$mu
    ## [1] 0.6325736
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1453724
    ## 
    ## $start.arg$mu
    ## [1] 0.4694766
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.1360989
    ## 
    ## $start.arg$mu
    ## [1] 1.234796
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 7584

    ## Cutoff for HT-1: 2 at quantile = 0.999

    ## Cutoff for HT-2: 10 at quantile = 0.999

    ## Cutoff for HT-3: 15 at quantile = 0.999

    ## Cutoff for HT-4: 11 at quantile = 0.999

    ## Cutoff for HT-5: 22 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## Round 1: 14812 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 2: 17065 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 3: 7463 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 4: 16376 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 5: 4431 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 6: 3917 Negative cells removed

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 7: 3848 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 8: 9959 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 9: 2272 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 10: 2536 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 11: 3303 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 12: 2407 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 13: 1641 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 14: 1598 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-5..."

    ## Round 15: 697 Negative cells removed

    ## [1] "No threshold found for HT-1..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## Round 16: 1176 Negative cells removed

    ## Round 17: Stopping, as no threshold extrema could be found.

    ## HT-1 cells: 998

    ## HT-2 cells: 620

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## HT-3 cells rescued: 3268

    ## HT-4 cells: 740

    ## [1] "No threshold found for HT-2..."
    ## [1] "No threshold found for HT-3..."
    ## [1] "No threshold found for HT-4..."
    ## [1] "No threshold found for HT-5..."

    ## HT-5 cells rescued: 5009

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning in max(.): no non-missing arguments to max; returning -Inf

    ## Warning: Removed 46 rows containing missing values (geom_point).

    ## Constructing seurat object for run s3b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 1.941319
    ## 
    ## $start.arg$mu
    ## [1] 30.60114
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5576581
    ## 
    ## $start.arg$mu
    ## [1] 3.490444
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.021632
    ## 
    ## $start.arg$mu
    ## [1] 13.25639
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.5876622
    ## 
    ## $start.arg$mu
    ## [1] 13.28543
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.624614
    ## 
    ## $start.arg$mu
    ## [1] 28.8136
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 15258

    ## Cutoff for HT-1: 48 at quantile = 0.8

    ## Cutoff for HT-2: 21 at quantile = 0.995

    ## Cutoff for HT-3: 76 at quantile = 0.999

    ## Cutoff for HT-4: 91 at quantile = 0.999

    ## Cutoff for HT-5: 190 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 2562 Negative cells removed

    ## Round 2: 206 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 3178 Negative cells removed

    ## Round 4: 517 Negative cells removed

    ## Round 5: 159 Negative cells removed

    ## Round 6: Done! All cells classified.

    ## HT-1 cells: 3290

    ## HT-2 cells: 2090

    ## HT-3 cells: 2334

    ## HT-4 cells: 1868

    ## HT-5 cells: 2957

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s4a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 12.51834
    ## 
    ## $start.arg$mu
    ## [1] 2.65494
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 24.02663
    ## 
    ## $start.arg$mu
    ## [1] 5.294374
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 5.772354
    ## 
    ## $start.arg$mu
    ## [1] 4.549455
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.188791
    ## 
    ## $start.arg$mu
    ## [1] 17.75699
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.795221
    ## 
    ## $start.arg$mu
    ## [1] 95.14064
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 27263

    ## Cutoff for HT-1: 4 at quantile = 0.8

    ## Cutoff for HT-2: 7 at quantile = 0.8

    ## Cutoff for HT-3: 14 at quantile = 0.999

    ## Cutoff for HT-4: 52 at quantile = 0.999

    ## Cutoff for HT-5: 245 at quantile = 0.99

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 23961 Negative cells removed

    ## Round 2: 5210 Negative cells removed

    ## Round 3: Done! All cells classified.

    ## HT-1 cells: 3418

    ## HT-2 cells: 5661

    ## HT-3 cells: 1176

    ## HT-4 cells: 4040

    ## HT-5 cells: 8264

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: did not converge in 10 iterations

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 1369150)

    ## Warning in max(.): no non-missing arguments to max; returning -Inf

    ## Constructing seurat object for run s4b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 3.018579
    ## 
    ## $start.arg$mu
    ## [1] 11.0267
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 2.984052
    ## 
    ## $start.arg$mu
    ## [1] 13.99061
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.529136
    ## 
    ## $start.arg$mu
    ## [1] 9.3195
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.839596
    ## 
    ## $start.arg$mu
    ## [1] 25.19017
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 3.686838
    ## 
    ## $start.arg$mu
    ## [1] 120.8595
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 21673

    ## Cutoff for HT-1: 16 at quantile = 0.8

    ## Cutoff for HT-2: 20 at quantile = 0.8

    ## Cutoff for HT-3: 29 at quantile = 0.999

    ## Cutoff for HT-4: 71 at quantile = 0.995

    ## Cutoff for HT-5: 330 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-1..."

    ## Round 1: 18251 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 2: 5395 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 3896 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 4: 4295 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 5: 3325 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 6: 2799 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 7: 2392 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 8: 2247 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 9: 1841 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 10: 1920 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 11: 1426 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 12: 758 Negative cells removed

    ## Round 13: 27 Negative cells removed

    ## Round 14: 2 Negative cells removed

    ## Round 15: 1 Negative cells removed

    ## Round 16: Done! All cells classified.

    ## HT-1 cells: 1047

    ## HT-2 cells: 863

    ## HT-3 cells: 705

    ## HT-4 cells: 827

    ## HT-5 cells: 870

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Constructing seurat object for run s5a

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 4.044071
    ## 
    ## $start.arg$mu
    ## [1] 96.01433
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.194757
    ## 
    ## $start.arg$mu
    ## [1] 5.450637
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.553368
    ## 
    ## $start.arg$mu
    ## [1] 9.690797
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.9665492
    ## 
    ## $start.arg$mu
    ## [1] 34.00716
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.18407
    ## 
    ## $start.arg$mu
    ## [1] 80.78738
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 14900

    ## Cutoff for HT-1: 138 at quantile = 0.8

    ## Cutoff for HT-2: 30 at quantile = 0.999

    ## Cutoff for HT-3: 42 at quantile = 0.999

    ## Cutoff for HT-4: 176 at quantile = 0.999

    ## Cutoff for HT-5: 434 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## [1] "No threshold found for HT-1..."

    ## Round 1: 7411 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 2: 7518 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 3: 120 Negative cells removed

    ## [1] "No threshold found for HT-1..."

    ## Round 4: Done! All cells classified.

    ## [1] "No threshold found for HT-1..."

    ## HT-1 cells rescued: 8461

    ## HT-2 cells: 1889

    ## HT-3 cells: 2046

    ## HT-4 cells: 2372

    ## HT-5 cells: 3437

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

    ## Warning: Removed 25 rows containing missing values (geom_point).

    ## Constructing seurat object for run s5b

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

    ## Demuxing hashtags with HTODemux to maximize singlets

    ## $start.arg
    ## $start.arg$size
    ## [1] 0.880844
    ## 
    ## $start.arg$mu
    ## [1] 3.671127
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.8538241
    ## 
    ## $start.arg$mu
    ## [1] 4.350469
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.326974
    ## 
    ## $start.arg$mu
    ## [1] 10.37911
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 0.6993285
    ## 
    ## $start.arg$mu
    ## [1] 12.89977
    ## 
    ## 
    ## $fix.arg
    ## NULL
    ## 
    ## $start.arg
    ## $start.arg$size
    ## [1] 1.779457
    ## 
    ## $start.arg$mu
    ## [1] 198.8209
    ## 
    ## 
    ## $fix.arg
    ## NULL

    ## Calculating cutoffs for maximum singlet recovery

    ## Maximum singlets = 11225

    ## Cutoff for HT-1: 23 at quantile = 0.999

    ## Cutoff for HT-2: 28 at quantile = 0.999

    ## Cutoff for HT-3: 39 at quantile = 0.99

    ## Cutoff for HT-4: 75 at quantile = 0.999

    ## Cutoff for HT-5: 1036 at quantile = 0.999

    ## Demuxing hashtags with deMULTIplex to rescue singlets

    ## Round 1: 4444 Negative cells removed

    ## Round 2: 98 Negative cells removed

    ## Round 3: 25 Negative cells removed

    ## Round 4: 2 Negative cells removed

    ## Round 5: Done! All cells classified.

    ## HT-1 cells: 2070

    ## HT-2 cells: 2124

    ## HT-3 cells: 2270

    ## HT-4 cells: 2077

    ## HT-5 cells: 2610

    ## [1] "Normalizing barode data..."
    ## [1] "Pre-allocating data structures..."
    ## [1] "Determining classifications for negative cells across all BCs, all q..."
    ## [1] "Computing classification stability..."
    ## [1] "Extracting rescued classifications..."

``` r
# Save data for next step
saveRDS(object = gems, file = "../data/gems.RDS")
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
    ##  date     2021-11-17                  
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
    ##  KernSmooth           * 2.23-20    2021-05-03 [2] CRAN (R 4.1.0)                                
    ##  knitr                  1.33       2021-04-24 [1] RSPM (R 4.1.0)                                
    ##  labeling               0.4.2      2020-10-20 [1] RSPM (R 4.1.0)                                
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
    ##  reshape2             * 1.4.4      2020-04-09 [1] RSPM (R 4.1.0)                                
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
    ##  sp                     1.4-5      2021-01-10 [1] RSPM (R 4.1.0)                                
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
    ##  vroom                  1.5.4      2021-08-05 [1] RSPM (R 4.1.0)                                
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