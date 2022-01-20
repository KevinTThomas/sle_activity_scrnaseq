`%nin%` <- purrr::compose(`!`, `%in%`)

#' @title PurgeXenoGenes
#'
#' @description Removes counts and gene loadings for a secondary species and strips the primary species
#' prefix from the remaining genes.  Useful when cleaning an object that had controls cells added for
#' multi-modal assays such as CITE-seq or REAP-seq.
#'
#' @param object Formerly mixed species Seurat object
#' @param primary_species_prefix Literal string or regex prefix denoting
#' genes from the species of interest
#' @param secondary_species_prefix Literal string or regex prefix
#' denoting genes from the species used as control
#'
#' @importFrom stringr str_detect str_remove
#' @importFrom Seurat CreateAssayObject
#' @importFrom SeuratObject LogSeuratCommand
#' @return Seurat object
#' @export
#'
#' @examples pbmc <- PurgeXenoGenes(object = pbmc, primary_species_prefix = "hg38-", secondary_species_prefix = "mm10-")
PurgeXenoGenes <-
    function(object,
    primary_species_prefix,
    secondary_species_prefix) {

        # clean assays
        for (i in names(object@assays)) {
            counts <- object[[i]]@counts[str_detect(
                string = rownames(object[[i]]@counts),
                pattern = secondary_species_prefix,
                negate = TRUE
            ), ]
            data <- object[[i]]@data[str_detect(
                string = rownames(object[[i]]@data),
                pattern = secondary_species_prefix,
                negate = TRUE
            ), ]
            meta.features <- object[[i]]@meta.features[str_detect(
                string = rownames(object[[i]]@meta.features),
                pattern = secondary_species_prefix,
                negate = TRUE
            ), ]
            rownames(counts) <- str_remove(string = rownames(counts), pattern = primary_species_prefix)
            rownames(data) <- str_remove(string = rownames(data), pattern = primary_species_prefix)
            rownames(meta.features) <- str_remove(string = rownames(meta.features), pattern = primary_species_prefix)

            replacement_assay <- CreateAssayObject(counts = counts)
            replacement_assay@data <- data
            replacement_assay@meta.features <- meta.features
            object[[i]] <- replacement_assay
        }

        # clean reductions
        for (j in names(object@reductions)) {
            feature.loadings <- object[[j]]@feature.loadings[str_detect(
                string = rownames(object[[j]]@feature.loadings),
                pattern = secondary_species_prefix,
                negate = TRUE
            ), ]
            rownames(feature.loadings) <-
                str_remove(string = rownames(feature.loadings), pattern = primary_species_prefix)
        }

        object <- SeuratObject::LogSeuratCommand(object = object, return.command = FALSE)
        object
    }

#' @title FindPCAElbow
#' @description Identify the upper bounds of
#' significant principal components
#'
#' @param object Seurat object
#' @param assay Assay to examine
#' @param slot Slot to use
#' @param reduction Dimensional reduction to look for (e.g. "pca" or "ica")
#'
#' @importFrom rsvd rpca
#' @importFrom SeuratObject GetAssayData
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#'
#' @return A copy of the object with the
#' largest significant principal component stored
#' as `max_pca_dim` in the object@misc slot
#' @export
#'
#' @examples
FindPCAElbow <-
    function(
      object,
      assay = "RNA",
      slot = "data",
      reduction = "pca",
      name = paste0("max_", reduction, "_dim"),
      n_pcs = 100,
      elbow_th = 0.025,
      perform_new_pca = FALSE
      ) {
        if (isTRUE(perform_new_pca) | !(reduction %in% names(object@reductions))) {
            expr_data <-
                GetAssayData(
                    object = object,
                    assay = assay,
                    slot = slot
                ) |>
                as.matrix()

            pca <-
                rpca(
                    A = t(expr_data),
                    k = n_pcs
                )

            pca_var_drop <- c(diff(pca$sdev^2), 0) / -pca$sdev^2
        } else {
            pca_var_drop <- c(diff(object[[reduction]]@stdev^2), 0) / -object[[reduction]]@stdev^2
        }
        # eval(parse(text = paste0("object@misc$", name, " <- rev(which(pca_sdev_drop > elbow_th))[1]")))

        pc_below_th <- which(pca_var_drop - elbow_th < 0)
        pc_th <- split(x = pc_below_th, f = cumsum(c(1, diff(pc_below_th) != 1))) %>%
            extract2(which(map_dbl(., length) > 1)[1]) |>
            extract(1) - 1
        eval(parse(text = paste0("object@misc$", name, " <- ", pc_th)))
        object
    }


#' @title constructGEMs
#' @description Constructs the Gene Expression Matrices (GEMs) with the given
#' parameters for multimodal data. Returns a final seurat object.
#' @param sample Character string to append onto each cell barcode.
#' @param gem_matrix Matrix containing gene expression data only. Assumes that
#' rows are genes and that columns are cells.
#' @param adt_matrix Matrix containing antibody capture data. If NULL, assumes
#' that no antibody capture data exists.
#' @param hashtags Character vector containing cell hashing tags (e.g. "HTO1")
#' found in the antibody capture data. If FALSE, assumes that no cell hashing
#' has taken place.
#'
#' @importFrom magrittr set_colnames
#' @importFrom Seurat NormalizeData
#' @importFrom SeuratObject AddMetaData CreateSeuratObject CreateAssayObject SetAssayData
#'
#' @return
#' @export
constructGEMs <- function(sample = NULL,
    gem_matrix,
    adt_matrix = NULL,
    norm_adt_matrix = NULL,
    hashtags = FALSE
    ) {
    # Add sample names to matrix
    if (!is.null(sample)) {
        rna_counts <-
            set_colnames(
                x = gem_matrix,
                value = paste(sample, colnames(gem_matrix), sep = "_")
            )
    } else {
        rna_counts <- gem_matrix
    }

    # If there's antibody data, make a separate matrix
    if (!is.null(adt_matrix)) {
        if (!is.null(sample)) {
            adt_counts <-
                set_colnames(
                    x = adt_matrix,
                    value = paste(sample, colnames(adt_matrix), sep = "_")
                )
            # If there's a normalized matrix, load that in, too
            if (!is.null(norm_adt_matrix)) {
                norm_adt_counts <- set_colnames(norm_adt_matrix, paste(sample, colnames(norm_adt_matrix), sep = "_"))
            }
        } else {
            adt_counts <- adt_matrix
            # If there's a normalized matrix, load that in, too
            if (!is.null(norm_adt_matrix)) {
                norm_adt_counts <- norm_adt_matrix
            }
        }

        # If there are hashtags in the antibody counts, make them a separate count matrix
        if (!is_false(hashtags) && any(hashtags %in% rownames(adt_counts))) {
            hto_counts <- adt_counts[hashtags[hashtags %in% rownames(adt_counts)], ]
            adt_counts <- adt_counts[setdiff(rownames(adt_counts), hashtags), ]
            if (!is.null(norm_adt_counts)) {
                norm_hto_counts <- norm_adt_counts[hashtags[hashtags %in% rownames(norm_adt_counts)], ]
                norm_adt_counts <- norm_adt_counts[setdiff(rownames(norm_adt_counts), hashtags), ]
            }
            # If the entire adt matrix is cell hashing, then the adt count matrix is now empty
            if (all(is.na(adt_counts))) {
                adt_counts <- NULL
                if (!is.null(norm_adt_counts)) {
                    norm_adt_counts <- NULL
                }
            }
        } else {
            hto_counts <- NULL
        }
    } else {
        adt_counts <- NULL
        hto_counts <- NULL
    }

    # Make the seurat object with the rna counts
    seurat_obj <- CreateSeuratObject(counts = rna_counts)

    # Make the adt assay object, if it exists
    if (!is.null(adt_counts)) {
        # If there are normalized counts, add them into the assay
        if (!is.null(norm_adt_counts)) {
            adt_obj <- CreateAssayObject(counts = adt_counts)
            adt_obj <- SetAssayData(object = adt_obj, slot = "data", new.data = norm_adt_counts)
        } else {
            adt_obj <- CreateAssayObject(counts = adt_counts)
        }
    } else {
        adt_obj <- NULL
    }

    # Make the hto assay object, if it exists
    if (!is.null(hto_counts)) {
        # If there are normalized counts, add them into the assay
        if (!is.null(norm_hto_counts)) {
            hto_obj <- CreateAssayObject(counts = hto_counts)
            hto_obj <- SetAssayData(object = hto_obj, slot = "data", new.data = norm_hto_counts)
        } else {
            hto_obj <- CreateAssayObject(counts = hto_counts)
        }
    } else {
        hto_obj <- NULL
    }

    # Check to see if the seurat object has >10 cells and >100 genes, and make the object
    if (ncol(seurat_obj) > 10 && nrow(seurat_obj) > 100) {
        # Add in multimodal data if it exists
        if (!is.null(adt_obj)) {
            seurat_obj[["CITE"]] <- adt_obj
        }
        if (!is.null(hto_obj)) {
            seurat_obj[["HTO"]] <- hto_obj
            seurat_obj <- NormalizeData(object = seurat_obj, assay = "HTO", normalization.method = "CLR")
            message("Demuxing hashtags with HTODemux to maximize singlets")
            seurat_obj <- HTODemux2(object = seurat_obj, assay = "HTO", positive.quantile.range = c(seq(0.80, 0.95, 0.05), 0.99, 0.995, 0.999))
            message("Demuxing hashtags with deMULTIplex to rescue singlets")
            seurat_obj <- AddMetaData(
                object = seurat_obj,
                metadata =
                    deMULTIplex(
                        counts = seurat_obj[["HTO"]]@counts,
                        plot.name = paste0("deMULTIplex_qc_plots/", sample, ".deMULTIplex.quantile.pdf")
                    )
            )
            seurat_singlets <- length(seurat_obj$hash.ID %nin% c("Doublet", "Negative"))
            deMULTIplex_singlets <- length(seurat_obj$deMULTIplex.calls.rescued %nin% c("Doublet", "Negative"))
            if (deMULTIplex_singlets >= seurat_singlets) {
                seurat_obj$final.HTO.ID <- seurat_obj$deMULTIplex.calls.rescued
            } else {
                seurat_obj$final.HTO.ID <- seurat_obj$hash.ID
            }
        }
        seurat_obj
    } else {
        NULL
    }
}

#' @title assignMetaData
#' @description Adds metadata to a seurat object from a samplesheet. Does this if the samples are singular or contain multiplex data. Returns a seurat object with metadata.
#' @param seurat_obj Seurat object for a sample to be annotated with metadata.
#' @param samplesheet Samplesheet containing per sample data to be added.
#' @param multiplex Character vector indicating what kind of multiplexing, if any, should be used. "none" assumes no multiplexing is occuring. "hash" looks for hashtags in the samplesheet and adds them according to exisiting metadata in the object for cell hashing. "demuxlet" indicates that samples have been demuxed in some other fashion (e.g. by genetics) and you are providing a dataframe indicating which cells are assigned to which identity.
#' @param barcode_df A dataframe containing barcode-sample assignments, such as the \code{[prefix].best} file generated by demuxlet.
#'
#' @importFrom dplyr filter select
#' @importFrom purrr reduce
#'
#' @return
#' @export
assignMetaData <- function(seurat_obj, samplesheet, multiplex = c("none", "hash", "demuxlet"), barcode_df = NULL) {
    # Determine the sample name, hooked at the beginning of the cell names
    sample <- gsub("_.*$", "", seurat_obj[["RNA"]]@counts@Dimnames[[2]][[1]])

    # If no multiplexing, assign by sample name
    if (multiplex == "none") {
        sample_meta <- samplesheet |>
            filter(sample_name == sample) |>
            select(-sample_name)
        for (x in names(sample_meta)) {
            seurat_obj[[x]] <- sample_meta[[x]]
        }
        seurat_obj@meta.data$run <- paste0(
            seurat_obj@meta.data[["run_batch"]],
            seurat_obj@meta.data[["sub_batch"]]
        )
        seurat_obj
    } else if (multiplex == "hash") { # If hashing, assign by hashtag
        seurat_obj@meta.data$run <- reduce(samplesheet$run[str_detect(sample, samplesheet$run)], intersect)
        sample_meta <-
            filter(
              .data = samplesheet,
              run == unique(seurat_obj@meta.data$run)
              ) |>
            select(-run)
        for (x in names(sample_meta)) {
            tag_meta <- character(length = ncol(seurat_obj))
            names(tag_meta) <- colnames(seurat_obj)
            for (y in sample_meta$hashtag) {
                tag_meta[seurat_obj$final.HTO.ID == y] <-
                    filter(
                      .data = sample_meta,
                      hashtag == y
                      ) |>
                    select(x)
            }
            seurat_obj[[x]] <- unlist(tag_meta)
        }
        seurat_obj
    }
}

#' @title getClosestFactor
#' @description Gets the closest factor pair for an integer
#' @param int A positive integer as a numeric value.
#'
#' @return
#' @export
getClosestFactors <- function(int) {
    if (length(int) != 1) {
        message("Argument must be length 1")
        stop()
    }
    if (int %% 1 != 0) {
        message("Argument must be a numeric integer")
        stop()
    }
    testNum <- floor(sqrt(int))
    while (int %% testNum != 0) {
        testNum <- testNum - 1
    }

    c(testNum, int / testNum)
}

#' @title calcGEMClassification
#' @param list A list of gene expression matrices to feed the function
#' @param primary_species_prefix
#' @param secondary_species_prefix
#'
#' @importFrom dplyr left_join bind_rows mutate relocate pull row_number
#' @importFrom stringr str_detect str_remove_all
#' @importFrom Matrix colSums
#' @importFrom tibble enframe
#' @importFrom purrr map
#' @importFrom sp point.in.polygon
#' @importFrom rlang set_names
#' @importFrom tidyr replace_na
#' @importFrom magrittr extract
#'
#' @return
#' @export
calcGEMClassification <- function(list, primary_species_prefix, secondary_species_prefix) {
    # create df of all cells with both species counts
    df1 <- map(
        list,
        function(x) {
            left_join(
                x = colSums(x[str_detect(string = rownames(x), primary_species_prefix, negate = FALSE), ]) |>
                    enframe(name = "cells", value = paste0(str_remove_all(primary_species_prefix, "_"), "_counts")),
                y = colSums(x[str_detect(string = rownames(x), secondary_species_prefix, negate = FALSE), ]) |>
                    enframe(name = "cells", value = paste0(str_remove_all(secondary_species_prefix, "_"), "_counts")),
                by = "cells"
            )
        }
    ) |>
        bind_rows(.id = "run") |>
        mutate(id = row_number()) |>
        relocate(id)
    # secondary polygons
    mouse_cutoffs <- list(
        s1a = data.frame(
            x = c(0, 0, 25, 25, 7, 7),
            y = c(400, 58100, 58100, 11000, 900, 400)
        ),
        s1b = data.frame(
            x = c(0, 0, 25, 25, 7, 7),
            y = c(400, 58100, 58100, 11000, 900, 400)
        ),
        s2a = data.frame(
            x = c(0, 0, 25, 25, 7, 7),
            y = c(400, 58100, 58100, 11000, 900, 400)
        ),
        s2b = data.frame(
            x = c(0, 0, 30, 30, 7, 7),
            y = c(400, 58100, 58100, 11000, 900, 400)
        ),
        s3a = data.frame(
            x = c(0, 0, 20, 20, 10, 10),
            y = c(50, 5500, 5500, 1200, 200, 50)
        ),
        s3b = data.frame(
            x = c(0, 0, 700, 700, 40, 40),
            y = c(120, 40000, 40000, 8800, 300, 120)
        ),
        s4a = data.frame(
            x = c(0, 0, 80, 80, 10, 4),
            y = c(300, 14000, 14000, 5000, 1000, 300)
        ),
        s4b = data.frame(
            x = c(0, 0, 100, 100, 15, 10),
            y = c(200, 20000, 20000, 5000, 400, 200)
        ),
        s5a = data.frame(
            x = c(0, 0, 100, 100, 20, 20),
            y = c(100, 15000, 15000, 3000, 300, 100)
        ),
        s5b = data.frame(
            x = c(0, 0, 100, 100, 20),
            y = c(300, 15000, 15000, 3000, 300)
        )
    ) |>
        bind_rows(.id = "run")
    # primary species bounding polygon coordinates
    human_cutoffs <- list(
        s1a = data.frame(
            x = c(100, 10000, 10000, 100),
            y = c(70, 70, 0, 0)
        ),
        s1b = data.frame(
            x = c(100, 10000, 11000, 100),
            y = c(70, 70, 0, 0)
        ),
        s2a = data.frame(
            x = c(100, 10000, 10000, 100),
            y = c(80, 80, 0, 0)
        ),
        s2b = data.frame(
            x = c(100, 10000, 10000, 100),
            y = c(80, 80, 0, 0)
        ),
        s3a = data.frame(
            x = c(15, 10000, 10000, 15),
            y = c(150, 150, 0, 0)
        ),
        s3b = data.frame(
            x = c(50, 10000, 10000, 50),
            y = c(150, 150, 0, 0)
        ),
        s4a = data.frame(
            x = c(20, 10000, 10000, 20),
            y = c(225, 225, 0, 0)
        ),
        s4b = data.frame(
            x = c(50, 100, 10000, 10000, 50),
            y = c(150, 300, 300, 0, 0)
        ),
        s5a = data.frame(
            x = c(25, 10000, 10000, 25),
            y = c(150, 150, 0, 0)
        ),
        s5b = data.frame(
            x = c(30, 10000, 10000, 30),
            y = c(150, 150, 0, 0)
        )
    ) |>
        bind_rows(.id = "run")
    # Generate cell dataframes
    mm_cells <- map(
        .x = unique(df1[["run"]]),
        .f = function(r) {
            mm_indices <- sp::point.in.polygon(
                point.x = filter(.data = df1, run == r) |>
                    pull(eval(parse(text = paste0("`", str_remove_all(primary_species_prefix, "_"), "_counts`")))),
                point.y = filter(.data = df1, run == r) |>
                    pull(eval(parse(text = paste0("`", str_remove_all(secondary_species_prefix, "_"), "_counts`")))),
                pol.x = filter(.data = mouse_cutoffs, run == r) |> pull(x),
                pol.y = filter(.data = mouse_cutoffs, run == r) |> pull(y)
            )

            mm_indices <- which(mm_indices != 0)
            barcodes <- filter(.data = df1, run == r) |>
            extract(mm_indices, c("id", "run", "cells"))

            barcodes
        }
    ) |>
        set_names(unique(df1[["run"]])) |>
        bind_rows(.id = "run") |>
        unique()
    hs_cells <- map(
        .x = unique(df1[["run"]]),
        .f = function(r) {
            hs_indices <- point.in.polygon(
                point.x =
                    filter(
                      .data = df1,
                      run == r
                    ) |>
                    pull(eval(parse(text = paste0("`", str_remove_all(primary_species_prefix, "_"), "_counts`")))),
                point.y =
                    filter(
                        .data = df1,
                        run == r
                        ) |>
                    pull(eval(parse(text = paste0("`", str_remove_all(secondary_species_prefix, "_"), "_counts`")))),
                pol.x =
                    filter(.data = human_cutoffs, run == r) |>
                    pull(x),
                pol.y =
                    filter(.data = human_cutoffs, run == r) |>
                    pull(y)
            )

            hs_indices <- which(hs_indices != 0)
            barcodes <- filter(.data = df1, run == r) |>
                extract(hs_indices, c("id", "run", "cells"))
            barcodes
        }
    ) |>
        set_names(unique(df1$run)) |>
        bind_rows(.id = "run")
    # Join with original df
    df1 <- left_join(
        x = df1,
        y = bind_rows(
            set_names(
                x = list(hs_cells, mm_cells),
                nm = c(str_remove_all(string = primary_species_prefix, pattern = "_"), str_remove_all(string = secondary_species_prefix, pattern = "_"))
            ),
            .id = "call"
        )
    ) |>
        mutate(call = replace_na(call, "Multiplet"))
    df1
}


#' Title
#'
#' @param object
#' @param assay
#' @param ig.hc
#' @param ig.lc
#' @param conserve.memory
#'
#' @importFrom Seurat SCTransform AddModuleScore CellCycleScoring
#' @importFrom magrittr extract2
#' @importFrom utils combn
#' @importFrom stringr str_remove
#' @importFrom dplyr bind_cols
#'
#' @return
#' @export
#'
#' @examples
ScoreIgDiffs <- function(object, assay = "RNA", ig.hc, ig.lc, conserve.memory = FALSE) {
    # Normalize the object
    message("Normalizing data...")
    # future::plan(strategy = future::multicore, workers = ceiling(0.25*parallel::detectCores()))
    norm_obj <- SCTransform(
        object = object,
        assay = assay,
        method = "glmGamPoi",
        conserve.memory = conserve.memory
    )
    # future::plan(strategy = future::multicore, workers = ceiling(parallel::detectCores()))

    # Score each heavy chain gene as a module
    message("Scoring heavy chain differences...")
    ig.hc.scores <- AddModuleScore(
        object = norm_obj,
        assay = "SCT",
        features = ig.hc[ig.hc %in% rownames(norm_obj[["SCT"]])],
        name = paste0(ig.hc[ig.hc %in% rownames(norm_obj[["SCT"]])], "."),
        nbin = 20
    ) %>%
        extract2(grep(pattern = "[.][0-9]+$", x = colnames(x = .[[]]), value = TRUE))

    # Calculate difference between all heavy chain module scores
    hc.diff <- combn(
        x = colnames(ig.hc.scores),
        m = 2,
        FUN = function(x) {
            ig.hc.scores[, x[1]] - ig.hc.scores[, x[2]]
        }
    )
    colnames(hc.diff) <- combn(
        x = colnames(ig.hc.scores),
        m = 2,
        FUN = function(x) {
            paste0(
              str_remove(
                string = x[1],
                pattern = "[.][0-9]+$"
                ),
              ".",
              str_remove(
                string = x[2],
                pattern = "[.][0-9]+$"
                ),
              ".diff"
              )
        }
    )
    rownames(hc.diff) <- rownames(ig.hc.scores)

    # Score each light chain gene as a module
    message("Scoring light chain differences...")
    ig.lc.scores <- AddModuleScore(
        object = norm_obj,
        assay = "SCT",
        features = ig.lc[ig.lc %in% rownames(norm_obj[["SCT"]])],
        name = paste0(ig.lc[ig.lc %in% rownames(norm_obj[["SCT"]])], "."),
        nbin = 20
    ) %>%
        extract2(grep(pattern = "[.][0-9]+$", x = colnames(x = .[[]]), value = TRUE))

    # Calculate difference between each light chain module scores
    lc.diff <- combn(
        x = colnames(ig.lc.scores),
        m = 2,
        FUN = function(x) {
            ig.lc.scores[, x[1]] - ig.lc.scores[, x[2]]
        }
    )
    colnames(lc.diff) <- combn(
        x = colnames(ig.lc.scores),
        m = 2,
        FUN = function(x) {
            paste0(
              str_remove(string = x[1], pattern = "[.][0-9]+$"), ".", str_remove(string = x[2], pattern = "[.][0-9]+$"), ".diff")
        }
    )
    rownames(lc.diff) <- rownames(ig.lc.scores)

    # Score cell cycle modules
    message("Scoring cell cycle differences...")
    cc.scores <- CellCycleScoring(
        object = norm_obj,
        s.features = cc.genes.updated.2019$s.genes,
        g2m.features = cc.genes.updated.2019$g2m.genes,
        set.ident = FALSE,
        nbin = 20
    ) |>
        extract2(c("S.Score", "G2M.Score", "Phase"))

    # Find difference in modules
    cc.scores$cc.diff <- cc.scores$S.Score - cc.scores$G2M.Score

    # Add to the metadata
    object@meta.data <- bind_cols(object@meta.data, as.data.frame(hc.diff), as.data.frame(lc.diff), as.data.frame(cc.scores[, c("cc.diff", "Phase")]))

    object
}


#' @title HTODemux2
#'
#' @param object
#' @param assay
#' @param positive.quantile.range
#' @param init
#' @param nstarts
#' @param kfunc
#' @param nsamples
#' @param seed
#' @param verbose
#'
#' @importFrom SeuratObject DefaultAssay GetAssayData Idents WhichCells
#' @importFrom Seurat AverageExpression
#' @importFrom rlang %||%
#' @importFrom stats quantile
#' @importFrom future.apply future_apply
#' @importFrom Matrix colSums
#' @importFrom fitdistrplus fitdist
#' @importFrom stats kmeans
#' @importFrom cluster clara
#'
#' @return
#' @export
#'
#' @examples
HTODemux2 <- function(object, assay = "HTO", positive.quantile.range = seq(0.8, 0.95, 0.05), init = NULL,
    nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42,
    verbose = TRUE) {
    if (!is.null(x = seed)) {
        set.seed(seed = seed)
    }
    assay <- assay %||% DefaultAssay(object = object)
    data <- GetAssayData(object = object, assay = assay)
    counts <- GetAssayData(object = object, assay = assay, slot = "counts")[
        ,
        colnames(x = object)
    ]
    counts <- as.matrix(x = counts)
    ncenters <- init %||% (nrow(x = data) + 1)
    switch(EXPR = kfunc,
        kmeans = {
            init.clusters <- kmeans(x = t(x = GetAssayData(
                object = object,
                assay = assay
            )), centers = ncenters, nstart = nstarts)
            Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
        },
        clara = {
            init.clusters <- clara(x = t(x = GetAssayData(
                object = object,
                assay = assay
            )), k = ncenters, samples = nsamples)
            Idents(
                object = object, cells = names(x = init.clusters$clustering),
                drop = TRUE
            ) <- init.clusters$clustering
        },
        stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'")
    )
    average.expression <- AverageExpression(
        object = object,
        assays = assay, verbose = FALSE
    )[[assay]]
    if (sum(average.expression == 0) > 0) {
        stop("Cells with zero counts exist as a cluster.")
    }
    discrete <- GetAssayData(object = object, assay = assay)
    # Set up
    cutoffs <- vector(mode = "list", length = nrow(x = data))
    names(cutoffs) <- rownames(x = data)
    # Calculate cutoff counts for each hashtag in the data
    for (iter in rownames(x = data)) {
        # Get the vector of count values
        values <- counts[iter, colnames(object)]
        # Take the values of cells in the cluster with the lowest average expression of the hashtag
        values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, ])]])]
        # Model the expression of the negative cluster cells with a negative binomial distribution
        fit <- suppressWarnings(expr = fitdist(data = values.use, distr = "nbinom"))
        # Get cutoffs for all positive quantiles for the hashtag
        cutoffs[[iter]] <- as.numeric(x = quantile(x = fit, probs = positive.quantile.range)$quantiles)
    }
    # List each combination of cutoffs as a dataframe
    cutoffs <- expand.grid(cutoffs)
    # Calculate the number of singlets for each of the cutoff values
    if (verbose) {
        message("Calculating cutoffs for maximum singlet recovery")
    }
    singlets <- future_apply(
        X = cutoffs,
        MARGIN = 1,
        FUN = function(co) {
            discrete[discrete > 0] <- 0
            for (iter in colnames(cutoffs)) {
                discrete[iter, names(x = which(x = counts[iter, colnames(object)] > co[[iter]]))] <- 1
            }
            npositive <- colSums(x = discrete)
            table(npositive)[["1"]]
        }
    )
    sample_quantiles <- rep(x = list(positive.quantile.range), times = nrow(x = data)) |>
        set_names(nm = rownames(x = data)) |>
        expand.grid()
    if (verbose) {
        message(paste0("Maximum singlets = ", max(singlets)))
        for (i in seq_along(cutoffs[which.max(singlets), ][1, ])) {
            message(paste0("Cutoff for ", names(cutoffs)[i], ": ", cutoffs[which.max(singlets), ][1, i], " at quantile = ", sample_quantiles[which.max(singlets), ][1, i]))
        }
    }
    # Calculate the number of positive hashtags for each cell at the cutoff values yielding the maximum number of singlets
    discrete[discrete > 0] <- 0
    for (iter in colnames(cutoffs[which.max(singlets), ][1, ])) {
        discrete[iter, names(x = which(x = counts[iter, colnames(object)] > cutoffs[which.max(singlets), ][1, iter]))] <- 1
    }
    npositive <- colSums(x = discrete)
    # Assign global classifications
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive > 1] <- "Doublet"
    donor.id <- rownames(x = data)
    hash.max <- apply(X = data, MARGIN = 2, FUN = max)
    hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
    hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
    hash.maxID <- as.character(x = donor.id[sapply(
        X = 1:ncol(x = data),
        FUN = function(x) {
            return(which(x = data[, x] == hash.max[x])[1])
        }
    )])
    hash.secondID <- as.character(x = donor.id[sapply(
        X = 1:ncol(x = data),
        FUN = function(x) {
            return(which(x = data[, x] == hash.second[x])[1])
        }
    )])
    hash.margin <- hash.max - hash.second
    doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
        paste(sort(x = c(hash.maxID[x], hash.secondID[x])),
            collapse = "_"
        )
    })
    classification <- classification.global
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global ==
        "Singlet")]
    classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global ==
        "Doublet")]
    classification.metadata <- data.frame(
        hash.maxID, hash.secondID,
        hash.margin, classification, classification.global
    )
    colnames(x = classification.metadata) <- paste(assay, c(
        "maxID",
        "secondID", "margin", "classification", "classification.global"
    ),
    sep = "_"
    )
    object <- AddMetaData(object = object, metadata = classification.metadata)
    Idents(object) <- paste0(assay, "_classification")
    doublets <- rownames(x = object[[]])[which(object[[paste0(
        assay,
        "_classification.global"
    )]] == "Doublet")]
    Idents(object = object, cells = doublets) <- "Doublet"
    object$hash.ID <- Idents(object = object)
    object
}

# Input: count matrix of HTO counts, name for output plots
# Output: dataframe and pdf
#' @title deMULTIplex
#'
#' @param counts
#' @param plot.name
#'
#' @importFrom deMULTIplex classifyCells findThresh findQ rescueCells findReclassCells
#' @importFrom ggplot2 ggplot aes geom_line theme labs geom_vline scale_color_manual coord_fixed xlim ylim geom_hline geom_errorbar
#' @importFrom data.table foverlaps data.table setkey
#' @importFrom magrittr use_serise
#' @importFrom grDevices pdf dev
#' @importFrom cowplot plot_grid
#'
#' @return
#' @export
#'
#' @examples
deMULTIplex <- function(counts, plot.name) {
    bar.table <- as.data.frame(t(counts))
    bar.table.full <- bar.table
    gg.list <- list()

    ## Set iter
    iter <- 1

    ## Perform Quantile Sweep
    bar.table_sweep.list <- list()
    n <- 0
    for (q in seq(0.01, 0.99, by = 0.02)) {
        n <- n + 1
        invisible(
            capture.output(
                bar.table_sweep.list[[n]] <- classifyCells(bar.table, q = q)
            )
        )
        names(bar.table_sweep.list)[n] <- paste("q=", q, sep = "")
    }

    ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
    threshold.results <- findThresh(call.list = bar.table_sweep.list)
    gg.list[[iter]] <- ggplot(data = threshold.results$res, aes(x = q, y = Proportion, color = Subset)) +
        geom_line() +
        theme(legend.position = "bottom") +
        labs(x = "quantile", y = "Proportion of cells", title = paste0("Round ", iter)) +
        geom_vline(xintercept = threshold.results$extrema, lty = 2) +
        scale_color_manual(values = c("red", "black", "blue")) +
        coord_fixed()

    ## Finalize classifications, remove negative cells
    round.calls <- classifyCells(bar.table, q = findQ(threshold.results$res, threshold.results$extrema))
    neg.cells <- names(round.calls)[which(round.calls == "Negative")]

    ## Check number of negative cells called, take them out of table
    num.neg.cells.in.round <- length(which(round.calls == "Negative"))
    if (num.neg.cells.in.round != 0) {
        message(paste0("Round ", iter, ": ", num.neg.cells.in.round, " Negative cells removed"))
        bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
    } else {
        message(paste0("Done! Took ", iter, " rounds to classify all cells."))
    }

    ## Repeat until all cells are classified
    while (num.neg.cells.in.round != 0) {
        ## Increase iter
        iter <- iter + 1

        ## Perform Quantile Sweep
        bar.table_sweep.list <- list()
        n <- 0
        for (q in c(seq(0.01, 0.99, by = 0.02))) {
            n <- n + 1
            invisible(
                capture.output(
                    bar.table_sweep.list[[n]] <- classifyCells(bar.table, q = q)
                )
            )
            names(bar.table_sweep.list)[n] <- paste("q=", q, sep = "")
        }

        ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
        threshold.results <- findThresh(call.list = bar.table_sweep.list)
        gg.list[[iter]] <- ggplot(data = threshold.results$res, aes(x = q, y = Proportion, color = Subset)) +
            geom_line() +
            theme(legend.position = "bottom") +
            labs(x = "quantile", y = "Proportion of cells", title = paste0("Round ", iter)) +
            geom_vline(xintercept = threshold.results$extrema, lty = 2) +
            scale_color_manual(values = c("red", "black", "blue")) +
            coord_fixed()

        if (length(threshold.results$extrema) == 0) {
            message(paste0("Round ", iter, ": Stopping, as no threshold extrema could be found."))
            break
        }

        ## Finalize classifications, remove negative cells
        round.calls <- classifyCells(bar.table, q = findQ(threshold.results$res, threshold.results$extrema))
        neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

        ## Check number of negative cells called, take them out of table
        num.neg.cells.in.round <- length(which(round.calls == "Negative"))
        if (num.neg.cells.in.round != 0) {
            message(paste0("Round ", iter, ": ", num.neg.cells.in.round, " Negative cells removed"))
            bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
        } else {
            message(paste0("Round ", iter, ": Done! All cells classified."))
        }
    }

    ## Finalize calls
    final.calls <- c(round.calls[which(round.calls != "Negative")], rep("Negative", length(neg.cells)))
    names(final.calls) <- c(names(round.calls[which(round.calls != "Negative")]), neg.cells)
    final.calls.rescued <- final.calls

    ## Make sure each tag has cells assigned to it, and rescue those tags that don't
    for (tag in colnames(bar.table.full)) {
        if (length(which(final.calls == tag)) == 0) {
            ## Perform Quantile Sweep
            bar.table_sweep.list <- list()
            n <- 0
            for (q in c(seq(0.01, 0.99, by = 0.02))) {
                n <- n + 1
                invisible(
                    capture.output(
                        bar.table_sweep.list[[n]] <- classifyCells(bar.table.full[neg.cells, ], q = q)
                    )
                )
                names(bar.table_sweep.list)[n] <- paste("q=", q, sep = "")
            }
            ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
            threshold.results <- findThresh(call.list = bar.table_sweep.list)
            round.calls <- classifyCells(bar.table.full[neg.cells, ], q = findQ(threshold.results$res, threshold.results$extrema))
            ## Add graph
            gg.list[[length(gg.list) + 1]] <- ggplot(data = threshold.results$res, aes(x = q, y = Proportion, color = Subset)) +
                geom_line() +
                theme(legend.position = "bottom") +
                labs(x = "quantile", y = "Proportion of cells", title = paste0("Rescue round for ", tag)) +
                geom_vline(xintercept = threshold.results$extrema, lty = 2) +
                scale_color_manual(values = c("red", "black", "blue")) +
                coord_fixed()
            ## Rescue these negative cells
            final.calls.rescued[names(round.calls[which(round.calls == tag)])] <- tag
            message(paste0(tag, " cells rescued: ", length(which(round.calls == tag))))
        } else {
            message(paste0(tag, " cells: ", length(which(final.calls == tag))))
            next
        }
    }

    ## Perform semi-supervised negative cell reclassification
    reclass.cells <- findReclassCells(bar.table.full, names(final.calls.rescued)[which(final.calls.rescued == "Negative")])
    invisible(
        capture.output(
            reclass.res <- rescueCells(bar.table.full, final.calls.rescued, reclass.cells)
        )
    )
    cutoff <- tryCatch(
        expr = foverlaps(
            x = data.table(
                ClassStability = reclass.res$ClassStability[-1],
                start = reclass.res$MatchRate_mean[-1] - 3 * reclass.res$MatchRate_sd[-1],
                end = reclass.res$MatchRate_mean[-1] + 3 * reclass.res$MatchRate_sd[-1]
            ) |>  na.omit(),
            y = data.table(
                ClassStability = reclass.res$ClassStability[1],
                start = reclass.res$MatchRate_mean[1] - 3 * reclass.res$MatchRate_sd[1],
                end = reclass.res$MatchRate_mean[1] + 3 * reclass.res$MatchRate_sd[1]
            ) |> setkey(start, end),
            type = "any",
            which = TRUE
        ) |>
            use_series(yid) |>
            is.na() |>
            which() |>
            max() + 1,
        error = function(e) {
            # If this throws an error, set cutoff to class stability value with a mean match rate closest to the true match rate
            dt <- data.table(
                ClassStability = reclass.res$ClassStability[-1],
                MatchRate_mean = reclass.res$MatchRate_mean[-1]
            ) |>
                setkey(MatchRate_mean)
            dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability
        }
    )
    # If cutoff is -Inf (because all the ranges overlap and max returns -Inf), take the closest value to the true mean
    if (cutoff == -Inf) {
        dt <- data.table(
            ClassStability = reclass.res$ClassStability[-1],
            MatchRate_mean = reclass.res$MatchRate_mean[-1]
        ) |>
            setkey(MatchRate_mean)
        cutoff <- dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability
    }
    # If cutoff is the max + 1 (because the ranges are so narrow that there are no overlaps), take the closest value to the true mean
    if (cutoff == max(reclass.res$ClassStability[-1]) + 1) {
        dt <- data.table(
            ClassStability = reclass.res$ClassStability[-1],
            MatchRate_mean = reclass.res$MatchRate_mean[-1]
        ) |>
            setkey(MatchRate_mean)
        cutoff <- dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability
    }

    ## Visualize results
    gg.list[[length(gg.list) + 1]] <- ggplot(reclass.res[-1, ], aes(x = ClassStability, y = MatchRate_mean)) +
        geom_point() +
        xlim(c(nrow(reclass.res) - 1, 1)) +
        ylim(c(0, reclass.res$MatchRate_mean[1] + 3 * reclass.res$MatchRate_sd[1] + 0.05)) +
        geom_errorbar(aes(ymin = MatchRate_mean - MatchRate_sd, ymax = MatchRate_mean + MatchRate_sd), width = .1) +
        geom_vline(xintercept = cutoff, col = "blue") +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1], color = "red") +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1] + 3 * reclass.res$MatchRate_sd[1], color = "red", lty = 2) +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1] - 3 * reclass.res$MatchRate_sd[1], color = "red", lty = 2) +
        labs(title = paste0("Class Stability Cutoff = ", cutoff))

    ## Confirm the cutoff (taken out for full automation)
    # print(gg.list[[length(gg.list)]])
    # cutoff <- as.integer(readline("Please confirm Class Stability Cutoff: "))

    ## Put new cutoff in place
    gg.list[[length(gg.list)]] <- ggplot(reclass.res[-1, ], aes(x = ClassStability, y = MatchRate_mean)) +
        geom_point() +
        xlim(c(nrow(reclass.res) - 1, 1)) +
        ylim(c(0, reclass.res$MatchRate_mean[1] + 3 * reclass.res$MatchRate_sd[1] + 0.05)) +
        geom_errorbar(aes(ymin = MatchRate_mean - MatchRate_sd, ymax = MatchRate_mean + MatchRate_sd), width = .1) +
        geom_vline(xintercept = cutoff, col = "blue") +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1], color = "red") +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1] + 3 * reclass.res$MatchRate_sd[1], color = "red", lty = 2) +
        geom_hline(yintercept = reclass.res$MatchRate_mean[1] - 3 * reclass.res$MatchRate_sd[1], color = "red", lty = 2) +
        labs(title = paste0("Class Stability Cutoff = ", cutoff))

    ## Finalize negative cell rescue results
    rescue.ind <- which(reclass.cells$ClassStability >= cutoff) ## Note: Value will be dataset-specific
    final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

    # Report of quantile selection
    pdf(file = plot.name, width = 13.333, height = 7.5)
    print(
        x = plot_grid(
            plotlist = gg.list
        )
    )
    dev.off()

    data.frame(
        deMULTIplex.calls = final.calls,
        deMULTIplex.calls.rescued = final.calls.rescued
    )
}

# Group small clusters into larger ones by setting a threshold
#' Title
#'
#' @param ids
#' @param SNN
#' @param threshold
#' @param group.singletons
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
GroupSmallClusters <- function(ids, SNN, threshold = 20, group.singletons = TRUE, verbose = TRUE) {
    # identify singletons
    singletons <- c()
    singletons <- names(x = which(x = table(ids) < threshold))
    singletons <- intersect(x = unique(x = ids), singletons)
    if (!group.singletons) {
        ids[which(ids %in% singletons)] <- "singleton"
        return(ids)
    }
    # calculate connectivity of singletons to non-singleton clusters, add singleton
    # to cluster it is most connected to
    cluster_names <- as.character(x = unique(x = ids))
    cluster_names <- setdiff(x = cluster_names, y = singletons)
    connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
    names(x = connectivity) <- cluster_names
    new.ids <- ids
    for (i in singletons) {
        i.cells <- names(which(ids == i))
        for (j in cluster_names) {
            j.cells <- names(which(ids == j))
            subSNN <- SNN[i.cells, j.cells]
            set.seed(1) # to match previous behavior, random seed being set in WhichCells
            if (is.object(x = subSNN)) {
                connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
            } else {
                connectivity[j] <- mean(x = subSNN)
            }
        }
        m <- max(connectivity, na.rm = T)
        mi <- which(x = connectivity == m, arr.ind = TRUE)
        closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
        ids[i.cells] <- closest_cluster
    }
    if (length(x = singletons) > 0 && verbose) {
        message(paste(
            length(x = singletons),
            "clusters with less than",
            threshold,
            "cells identified.",
            length(x = unique(x = ids)),
            "final clusters."
        ))
    }
    ids
}
