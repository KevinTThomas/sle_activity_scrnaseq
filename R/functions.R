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
#' @importFrom magrittr `%<>%`
#' @return Seurat object
#' @export
#'
#' @examples pbmc <- PurgeXenoGenes(object = pbmc, primary_species_prefix = "hg38-", secondary_species_prefix = "mm10-")
PurgeXenoGenes <- function(object,
                           primary_species_prefix,
                           secondary_species_prefix){
  # clean assays
  for (i in names(object@assays)){

    counts <- object[[i]]@counts[str_detect(string = rownames(object[[i]]@counts),
                                            pattern = secondary_species_prefix,
                                            negate = TRUE),]
    data <- object[[i]]@data[str_detect(string = rownames(object[[i]]@data),
                                        pattern = secondary_species_prefix,
                                        negate = TRUE),]
    meta.features <- object[[i]]@meta.features[str_detect(string = rownames(object[[i]]@meta.features),
                                                          pattern = secondary_species_prefix,
                                                          negate = TRUE),]
    rownames(counts) %<>% str_remove(pattern = primary_species_prefix)
    rownames(data) %<>% str_remove(pattern = primary_species_prefix)
    rownames(meta.features) %<>% str_remove(pattern = primary_species_prefix)

    replacement_assay <- CreateAssayObject(counts = counts)
    replacement_assay@data <- data
    replacement_assay@meta.features <- meta.features
    object[[i]] <- replacement_assay
  }

  # clean reductions
  for (j in names(object@reductions)){
    feature.loadings <- object[[j]]@feature.loadings[str_detect(string = rownames(object[[j]]@feature.loadings),
                                                                pattern = secondary_species_prefix,
                                                                negate = TRUE),]
    rownames(feature.loadings) %<>% str_remove(pattern = primary_species_prefix)
  }

  object <- SeuratObject::LogSeuratCommand(object = object)
  return(object)
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
#' @importFrom Seurat GetAssayData
#'
#' @return A copy of the object with the
#' largest significant principal component stored
#' as `max_pca_dim` in the object@misc slot
#' @export
#'
#' @examples
FindPCAElbow <- function(object,
                         assay = "RNA",
                         slot = "data",
                         reduction = "pca",
                         name = paste0("max_", reduction, "_dim"),
                         n_pcs = 100,
                         elbow_th = 0.025,
                         perform_new_pca = FALSE){

  if (isTRUE(perform_new_pca) | !(reduction %in% names(object@reductions))){
    expr_data <- GetAssayData(object,
                              assay = assay,
                              slot = slot) %>%
      as.matrix()
    pca <- rpca(A = t(expr_data),
                k = n_pcs)

    pca_var_drop <- c(diff(pca$sdev^2), 0) / -pca$sdev^2
  } else{
    pca_var_drop <- c(diff(object[[reduction]]@stdev^2), 0) / -object[[reduction]]@stdev^2
  }
  # eval(parse(text = paste0("object@misc$", name, " <- rev(which(pca_sdev_drop > elbow_th))[1]")))

  pc_below_th <- which(pca_var_drop-elbow_th<0)
  pc_th <- split(pc_below_th, cumsum(c(1, diff(pc_below_th) != 1))) %>%
    .[[which(map_dbl(.,length) > 1)[1]]] %>%
    .[1]-1
  eval(parse(text = paste0("object@misc$", name, " <- ", pc_th)))
  return(object)
}

#' @title bottom95_rna_count
#'
#' @description Calculate the bottom 95% percentile for nCount_RNA of a Seurat object
#' and store it in the `misc` slot
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
bottom99_rna_count <- function(object){
  object@misc$nCount_RNA99 <- quantile(object@meta.data[,'nCount_RNA'],
                                       seq(0,1,0.01))[['99%']]
  return(object)
}


#' @title CleanAssayFeatureNames
#'
#' @description Remove a pattern from all rownames within a Seurat Assay object
#'
#' @param object
#' @param pattern
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
CleanAssayFeatureNames <- function(object,
                                   pattern,
                                   assay){
  # clean assays
  rownames(object[[assay]]@counts) <- str_remove(string = rownames(object[[assay]]@counts),
                                                 pattern = pattern)
  rownames(object[[assay]]@data) <- str_remove(string = rownames(object[[assay]]@data),
                                               pattern = pattern)
  rownames(object[[assay]]@scale.data) <- str_remove(string = rownames(object[[assay]]@scale.data),
                                                     pattern = pattern)
  rownames(object[[assay]]@meta.features) <- str_remove(string = rownames(object[[assay]]@meta.features),
                                                        pattern = pattern)

  return(object)
}

#' @title adv_grouped_violins
#'
#' @description group subset of violin ensemble
#'
#' @param object Seurat object to plot
#' @param groups_show subset of clusters to plot
#' @param ... parameters to pass to ViolinEnsemble
#'
#' @return
#' @export
#'
#' @examples
adv_grouped_violins <- function(object, groups_show, alpha = 0.1, pt_size = 0.1, split_facet = "cols",...){
  plot_data <- ViolinEnsemble(object = object, ...) %>%
    `[[`("data") %>%
    filter(cluster %in% groups_show) %>%
    arrange(cluster)

  plots <- map(unique(plot_data$cluster), function(i){
    plot_data %>%
      filter(cluster == i) %>%
      ggplot(aes(x = ident,
                 y = value,
                 fill = feature)) +
      theme(legend.position = "none",
            strip.text = element_text(size = 9),
            axis.text = element_text(size = 9),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(.05,
                                 "lines"),
            panel.border = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(face = "italic", angle = 45),
            plot.title = element_text(hjust = 0.5, size = 9)) +
      coord_flip() +
      labs(x = "", y = "", title = i) +
      geom_violin(scale = "width",
                  trim = TRUE) +
      geom_jitter(
        alpha = alpha,
        size = pt_size) +
      facet_row(facets = vars(feature), scales = "free_x")
  })

  widths <- plot_data %>%
    select(cluster, feature) %>%
    group_by(cluster) %>%
    summarise(feature_count = length(unique(feature))) %>%
    pull(feature_count)

  if (split_facet == "cols"){
    pg <- plot_grid(plotlist = plots, nrow = 1, rel_widths = widths)
  } else if (split_facet == "rows"){
    pg <- plot_grid(plotlist = plots, ncol = 1, rel_widths = widths)
  }

  return(pg)
}

#' @title citeNormalize
#' @description A different method of normalizing cite_seq data which differs between runs. This method takes raw CITE-seq counts for each feature and divides them by the total HTO counts per cell.
#' @param cite.assay Assay variable pointing to the cite-seq data
#' @param hto.assay Assay variable pointing to the hash-tagging data
#' @param object A seurat object
#'
#' @return
#' @export
citeNormalize <- function (object, cite.assay, hto.assay) {
  cite.counts = GetAssayData(
    object = object,
    slot = "counts",
    assay = cite.assay
  )
  hto.counts <- GetAssayData(
    object = object,
    slot = "counts",
    assay = hto.assay
  ) %>% colSums()
  normalized.cite <- cite.counts/hto.counts[col(cite.counts)]
  object <- SetAssayData(
    object = object,
    slot = "counts",
    assay = cite.assay,
    new.data = normalized.cite
  )
  meta.data = data.frame(
    row.names = colnames(normalized.cite),
    nCITE.HTO = colSums(normalized.cite)
  ) %>%
    cbind(object@meta.data$nCount_CITE) %>%
    cbind(t(cite.counts)) %>%
    `colnames<-`(c("nCITE.HTO", "nRawCount_CITE", paste0("raw_", rownames(cite.counts))))
  object <- AddMetaData(object = object, metadata = meta.data)
  return(object)
}

#' @title clusterScan
#' @description Iterate through a series of clustering resolutions
#' and return the minimum number of differentially expressed genes
#' between neighboring clusters
#'
#' @param obj Seurat object
#' @param starting_res Initial resolution to start clustering at. Default: 0
#' @param increment Amount by which to raise the resolution each round of testing.
#' Default: 0.1
#' @param ending_res Final resolution to test.  Default: 1
#' @param differing_cluster_DEGs Minimum number of differentially expressed genes each
#' cluster must have as compared to its closest neighbor.  If any clusters have below
#' this number of DEGs, testing stops and an object with clustering performed using
#' the previous resolution is returned. Default: 10
#' @param fdr_thresh False discovery rate to use when testing DEGs.  Default: 0.05
#' @param assay Object assay to use for determining cluster DEGs.  Default: uses `SCT` and
#' that is not present, `RNA`
#' @param clustering_reduction Reduction to use when clustering.  Default: `harmony`,
#' falling back on `pca` if that is not available.
#' @param dims_use Number of dimensions to use for clustering.  Will use the maximum significant
#' dimensions if present in the `misc` slot, else all dimensions are used.  Default: the
#' value for `max_pca_dim` stored in the `object@misc` slot or all of the dimensions of
#' the `clustering_reduction`
#' @param verbose Be vocal about progress
#'
#' @return a named list of the minimum DEGs between neighboring clusters
#' @export
#'
#' @examples
clusterScan <- function(object,
                        starting_res = 0.1,
                        ending_res = 1,
                        increment = 0.1,
                        assay = NULL,
                        dims_use = NULL,
                        clustering_reduction = NULL,
                        verbose = FALSE){

  sc <- import(module = "scanpy", delay_load = TRUE)

  assay <- assay %||% ifelse(test = "SCT" %in% names(object@assays),
                             yes = "SCT",
                             no = "RNA")

  clustering_reduction <- clustering_reduction %||% ifelse(test = "harmony" %in% Reductions(object),
                                                           yes = "harmony",
                                                           no = "pca")

  dims_use <- dims_use %||% ifelse(test = "max_pca_dim" %in% object@misc,
                                   yes = object@misc[["max_pca_dim"]],
                                   no = ncol(object[[clustering_reduction]]@cell.embeddings))

  if (verbose) message("Converting object")
  adata <- convert_to_anndata(object,
                              assay = assay)

  if (verbose) message("Finding neighbors")
  sc$pp$neighbors(adata = adata,
                  use_rep = str_glue("X_{clustering_reduction}"),
                  n_pcs = dims_use)

  if (verbose) message("Finding clusters")
  clusters <- map_dfc(seq(from = starting_res,
                          to = ending_res,
                          by = increment),
                      function(i){
                        current_clusters <- sc$tl$leiden(adata,
                                                         resolution = i,
                                                         copy = TRUE)$obs["leiden"] %>%
                          `names<-`(paste0(tolower(assay),".res.", i)) %>%
                          as_tibble() %>%
                          mutate_if(is.factor,
                                    .funs = ~as_factor(as.integer(.))
                          )
                      })

  # We need to prevent any resolutions where there was only one cluster from being tested
  clusters_use <- clusters %>%
    pivot_longer(cols = everything(),
                 names_to = "res",
                 values_to = "cluster") %>%
    group_by(res) %>%
    summarize(count = length(unique(cluster))) %>%
    filter(count > 1) %>%
    pull(res)

  if (verbose) message(str_glue('Using: {glue::glue_collapse(clusters_use, sep = " ")}'))
  object@meta.data <- object@meta.data %>%
    as_tibble(rownames = "cell") %>%
    left_join(mutate(clusters[,clusters_use],
                     cell = rownames(object@meta.data)),
              by = "cell") %>%
    column_to_rownames("cell")

  object@misc$clusters_use <- clusters_use

  return(object)
}

# autoRes <- function(object,
#                     assay = NULL,
#                     clustering_reduction = NULL,
#                     differing_cluster_DEGs = 10,
#                     fdr_thresh = 0.05,
#                     verbose = TRUE){
#
#   assay <- assay %||% ifelse(test = "SCT" %in% names(object@assays),
#                              yes = "SCT",
#                              no = "RNA")
#
#   clustering_reduction <- clustering_reduction %||% ifelse(test = "harmony" %in% Reductions(object),
#                                                            yes = "harmony",
#                                                            no = "pca")
#
#   if (verbose) message("Beginning DEG testing")
#   deg_counts <- future_map_int(seq_along(object@misc$clusters_use),
#                                .progress = TRUE,
#                         function(i){
#                           if (verbose) message(paste("Testing ", i))
#                           Idents(object) <- object@meta.data[[object@misc$clusters_use[[i]]]]
#
#                           sCVdata <- CalcSCV(
#                             inD = object,
#                             cl = Idents(object),
#                             assayType = assay,
#                             DRforClust = clustering_reduction,
#                             exponent = exp(1),
#                             pseudocount = 1,
#                             DRthresh = 0.1,
#                             calcSil = FALSE,
#                             calcDEvsRest = TRUE,
#                             calcDEcombn = TRUE
#                           )
#
#                           DE_bw_NN <- sapply(DEneighb(sCVdata, fdr_thresh), nrow)
#                           # ^ counts # of DE genes between neighbouring clusters at 5% FDR
#                           min(DE_bw_NN)
#                         })
#
#   if (verbose) message("Finished DEG testing")
#
#   optimal_res <- object@misc$clusters_use[[max(which(deg_counts > differing_cluster_DEGs))]]
#
#   Idents(object = object) <- object@meta.data[[optimal_res]]
#   object@misc$optimal_res <- str_remove(string = optimal_res, pattern = "^res.")
#   object@misc$clusters_use <- NULL
#   return(object)
#
# }

autoRes <- function(object,
                    deg_counts,
                    differing_cluster_DEGs = 20,
                    fdr_thresh = 0.05,
                    verbose = TRUE){
  first_below_threshold = min(which(deg_counts < differing_cluster_DEGs))
  if (is.finite(first_below_threshold)){
    optimal_res <- object@misc$clusters_use[[max(first_below_threshold-1,1)]]

    Idents(object = object) <- object@meta.data[[optimal_res]]
    object@misc$optimal_res <- str_remove(string = optimal_res, pattern = "^res.")
    return(object)
  } else {
    message("No resolutions were found matching your criteria")
    stop()
  }
}

#' @title constructGEMs
#' @description Constructs the Gene Expression Matrices (GEMs) with the given parameters for multimodal data. Returns a final seurat object.
#' @param sample Character string to append onto each cell barcode.
#' @param gem_matrix Matrix containing gene expression data only. Assumes that rows are genes and that columns are cells.
#' @param adt_matrix Matrix containing antibody capture data. If NULL, assumes that no antibody capture data exists.
#' @param hashtags Character vector containing cell hashing tags (e.g. "HTO1") found in the antibody capture data. If FALSE, assumes that no cell hashing has taken place.
#'
#' @return
#' @export
constructGEMs <- function(sample=NULL, gem_matrix, adt_matrix=NULL, norm_adt_matrix=NULL, hashtags=FALSE) {
  # Add sample names to matrix
  if (!is.null(sample)) {
    rna_counts <- `colnames<-`(gem_matrix, paste(sample, colnames(gem_matrix), sep = "_"))
  }else{
    rna_counts <- gem_matrix
  }

  # If there's antibody data, make a separate matrix
  if (!is.null(adt_matrix)) {
    if (!is.null(sample)) {
      adt_counts <- `colnames<-`(adt_matrix,paste(sample,colnames(adt_matrix),sep = "_"))
      # If there's a normalized matrix, load that in, too
      if (!is.null(norm_adt_matrix)) {
        norm_adt_counts <- `colnames<-`(norm_adt_matrix,paste(sample,colnames(norm_adt_matrix),sep = "_"))
      }
    }else{
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
    }else{
      hto_counts <- NULL
    }
  }else{
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
      adt_obj <- SetAssayData(adt_obj, slot = "data", new.data = norm_adt_counts)
    }else{
      adt_obj <- CreateAssayObject(counts = adt_counts)
    }
  }else{
    adt_obj <- NULL
  }

  # Make the hto assay object, if it exists
  if (!is.null(hto_counts)) {
    # If there are normalized counts, add them into the assay
    if (!is.null(norm_hto_counts)) {
      hto_obj <- CreateAssayObject(counts = hto_counts)
      hto_obj <- SetAssayData(hto_obj, slot = "data", new.data = norm_hto_counts)
    }else{
      hto_obj <- CreateAssayObject(counts = hto_counts)
    }
  }else{
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
      message(paste0("Demuxing hashtags with HTODemux to maximize singlets"))
      seurat_obj <- HTODemux2(object = seurat_obj, assay = "HTO", positive.quantile.range = c(seq(0.80, 0.95, 0.05),0.99,0.995,0.999))
      message(paste0("Demuxing hashtags with deMULTIplex to rescue singlets"))
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
      }else{
        seurat_obj$final.HTO.ID <- seurat_obj$hash.ID
      }
    }
    return(seurat_obj)
  }else{
    return(NULL)
  }
}

#' @title assignMetaData
#' @description Adds metadata to a seurat object from a samplesheet. Does this if the samples are singular or contain multiplex data. Returns a seurat object with metadata.
#' @param seurat_obj Seurat object for a sample to be annotated with metadata.
#' @param samplesheet Samplesheet containing per sample data to be added.
#' @param multiplex Character vector indicating what kind of multiplexing, if any, should be used. "none" assumes no multiplexing is occuring. "hash" looks for hashtags in the samplesheet and adds them according to exisiting metadata in the object for cell hashing. "demuxlet" indicates that samples have been demuxed in some other fashion (e.g. by genetics) and you are providing a dataframe indicating which cells are assigned to which identity.
#' @param barcode_df A dataframe containing barcode-sample assignments, such as the \code{[prefix].best} file generated by demuxlet.
#'
#' @return
#' @export
assignMetaData <- function(seurat_obj, samplesheet, multiplex=c("none", "hash", "demuxlet"), barcode_df=NULL){
  # Determine the sample name, hooked at the beginning of the cell names
  sample = gsub("_.*$", "", seurat_obj[["RNA"]]@counts@Dimnames[[2]][[1]])

  # If no multiplexing, assign by sample name
  if (multiplex == "none") {
    sample_meta <- samplesheet %>%
      dplyr::filter(sample_name == sample) %>%
      dplyr::select(-sample_name)
    for (x in names(sample_meta)) {
      seurat_obj[[x]] <- sample_meta[[x]]
    }
    seurat_obj@meta.data$run <- paste0(
      seurat_obj@meta.data$run_batch,
      seurat_obj@meta.data$sub_batch
    )
    return(seurat_obj)
  }
  # If hashing, assign by hashtag
  if (multiplex == "hash") {
    seurat_obj@meta.data$run <- purrr::reduce(samplesheet$run[str_detect(sample, samplesheet$run)], intersect)
    sample_meta <- samplesheet %>%
      dplyr::filter(run == unique(seurat_obj@meta.data$run)) %>%
      dplyr::select(-run)
    for (x in names(sample_meta)) {
      tag_meta = character(length = ncol(seurat_obj))
      names(tag_meta) = colnames(seurat_obj)
      for (y in sample_meta$hashtag) {
        tag_meta[seurat_obj$final.HTO.ID == y] <- sample_meta %>%
          dplyr::filter(hashtag == y) %>%
          dplyr::select(x)
      }
      seurat_obj[[x]] <- unlist(tag_meta)
    }
    return(seurat_obj)
  }
}

#' @title getClosestFactor
#' @description Gets the closest factor pair for an integer
#' @param int A positive integer as a numeric value.
#'
#' @return
#' @export
getClosestFactors <- function(int){
  if (length(int) != 1) {
    message("Argument must be length 1")
    stop()
  }
  if (int%%1 != 0) {
    message("Argument must be a numeric integer")
    stop()
  }
  testNum = floor(sqrt(int))
  while (int%%testNum != 0) {
    testNum <- testNum-1
  }
  return(c(testNum, int/testNum))
}

#' @title calcGEMClassification
#' @param list A list of gene expression matrices to feed the function
#' @param primary_species_prefix
#' @param secondary_species_prefix
#' @return
#' @export
calcGEMClassification <- function(list, primary_species_prefix, secondary_species_prefix) {
  # create df of all cells with both species counts
  df1 <- map(
    list,
    function(x){
      left_join(
        x = x[str_detect(string = rownames(x),primary_species_prefix, negate = FALSE), ] %>%
          colSums() %>%
          enframe(name = "cells", value = paste0(str_remove_all(primary_species_prefix, "_"), "_counts")),
        y = x[str_detect(string = rownames(x),secondary_species_prefix, negate = FALSE), ] %>%
          colSums() %>%
          enframe(name = "cells", value = paste0(str_remove_all(secondary_species_prefix, "_"), "_counts")),
        by = "cells"
      )
    }
  ) %>%
    bind_rows(.id = "run") %>%
    mutate(id = row_number()) %>%
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
  ) %>%
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
  ) %>%
    bind_rows(.id = "run")
  # Generate cell dataframes
  mm_cells <- map(
    .x = unique(df1$run),
    .f = function (r) {
      mm_indices <- sp::point.in.polygon(
        point.x = df1 %>%
          filter(run == r) %>%
          pull(eval(parse(text=paste0("`", str_remove_all(primary_species_prefix, "_"), "_counts`")))),
        point.y = df1 %>%
          filter(run == r) %>%
          pull(eval(parse(text=paste0("`", str_remove_all(secondary_species_prefix, "_"), "_counts`")))),
        pol.x = mouse_cutoffs %>%
          filter(run == r) %>%
          pull(x),
        pol.y = mouse_cutoffs %>%
          filter(run == r) %>%
          pull(y)
      ) %>%
        (function(x) which(x != 0))
      barcodes <- df1%>%
        filter(run == r) %>%
        .[mm_indices, c("id","run","cells")]
      return(barcodes)
    }
  ) %>%
    `names<-`(unique(df1$run)) %>%
    bind_rows(.id = "run") %>% unique()
  hs_cells <- map(
    .x = unique(df1$run),
    .f = function (r) {
      hs_indices <- sp::point.in.polygon(
        point.x = df1 %>%
          filter(run == r) %>%
          pull(eval(parse(text=paste0("`", str_remove_all(primary_species_prefix, "_"), "_counts`")))),
        point.y = df1 %>%
          filter(run == r) %>%
          pull(eval(parse(text=paste0("`", str_remove_all(secondary_species_prefix, "_"), "_counts`")))),
        pol.x = human_cutoffs %>%
          filter(run == r) %>%
          pull(x),
        pol.y = human_cutoffs %>%
          filter(run == r) %>%
          pull(y)
      ) %>%
        (function(x) which(x != 0))
      barcodes <- df1 %>%
        filter(run == r) %>%
        .[hs_indices, c("id","run","cells")]
      return(barcodes)
    }
  ) %>%
    `names<-`(unique(df1$run)) %>%
    bind_rows(.id = "run")
  # Join with original df
  df1 <- left_join(
    x = df1,
    y = bind_rows(
      `names<-`(
        x = list(hs_cells, mm_cells),
        value = c(str_remove_all(primary_species_prefix, "_"), str_remove_all(secondary_species_prefix, "_"))
      ),
      .id = "call"
    )
  ) %>%
    mutate(call = replace_na(call, "Multiplet"))
  return(df1)

}

ScoreIgDiffs <- function (object, assay = "RNA", ig.hc, ig.lc, conserve.memory = FALSE) {
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
    .[[grep(pattern = "[.][0-9]+$", x = colnames(x = .[[]]), value = TRUE)]]

  # Calculate difference between all heavy chain module scores
  hc.diff <- combn(
    x = colnames(ig.hc.scores),
    m = 2,
    FUN = function (x) {ig.hc.scores[,x[1]] - ig.hc.scores[,x[2]]}
  )
  colnames(hc.diff) <- combn(
    x = colnames(ig.hc.scores),
    m = 2,
    FUN = function (x) {paste0(str_remove(string = x[1], pattern = "[.][0-9]+$"), ".", str_remove(string = x[2], pattern = "[.][0-9]+$"), ".diff")}
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
    .[[grep(pattern = "[.][0-9]+$", x = colnames(x = .[[]]), value = TRUE)]]

  # Calculate difference between each light chain module scores
  lc.diff <- combn(
    x = colnames(ig.lc.scores),
    m = 2,
    FUN = function (x) {ig.lc.scores[,x[1]] - ig.lc.scores[,x[2]]}
  )
  colnames(lc.diff) <- combn(
    x = colnames(ig.lc.scores),
    m = 2,
    FUN = function (x) {paste0(str_remove(string = x[1], pattern = "[.][0-9]+$"), ".", str_remove(string = x[2], pattern = "[.][0-9]+$"), ".diff")}
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
  ) %>%
    .[[c("S.Score", "G2M.Score", "Phase")]]

  # Find difference in modules
  cc.scores$cc.diff <- cc.scores$S.Score - cc.scores$G2M.Score

  # Add to the metadata
  object@meta.data <- bind_cols(object@meta.data, as.data.frame(hc.diff), as.data.frame(lc.diff), as.data.frame(cc.scores[,c("cc.diff", "Phase")]))

  return(object)
}

HTODemux2 <- function (object, assay = "HTO", positive.quantile.range = seq(0.8, 0.95, 0.05), init = NULL,
                       nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42,
                       verbose = TRUE)
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[,
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object,
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- cluster::clara(x = t(x = GetAssayData(object = object,
                                                           assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering),
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object,
                                          assays = assay, verbose = FALSE)[[assay]]
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
    fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values.use, distr = "nbinom"))
    # Get cutoffs for all positive quantiles for the hashtag
    cutoffs[[iter]] <- as.numeric(x = quantile(x = fit, probs = positive.quantile.range)$quantiles)
  }
  # List each combination of cutoffs as a dataframe
  cutoffs <- expand.grid(cutoffs)
  # Calculate the number of singlets for each of the cutoff values
  if (verbose) {message("Calculating cutoffs for maximum singlet recovery")}
  singlets <- future.apply::future_apply(
    X = cutoffs,
    MARGIN = 1,
    FUN = function (co) {
      discrete[discrete > 0] <- 0
      for (iter in colnames(cutoffs)) {
        discrete[iter, names(x = which(x = counts[iter, colnames(object)] > co[[iter]]))] <- 1
      }
      npositive <- colSums(x = discrete)
      singlets <- table(npositive)[["1"]]
      return(singlets)
    }
  )
  quantiles <- rep(x = list(positive.quantile.range), times = nrow(x = data)) %>%
    `names<-`(value = rownames(x = data)) %>%
    expand.grid()
  if (verbose) {
    message(paste0("Maximum singlets = ", max(singlets)))
    for (i in seq_along(cutoffs[which.max(singlets),][1,])) {
      message(paste0("Cutoff for ", names(cutoffs)[i], ": ", cutoffs[which.max(singlets),][1,i], " at quantile = ", quantiles[which.max(singlets),][1,i]))
    }
  }
  # Calculate the number of positive hashtags for each cell at the cutoff values yielding the maximum number of singlets
  discrete[discrete > 0] <- 0
  for (iter in colnames(cutoffs[which.max(singlets),][1,])) {
    discrete[iter, names(x = which(x = counts[iter, colnames(object)] > cutoffs[which.max(singlets),][1,iter]))] <- 1
  }
  npositive <- colSums(x = discrete)
  # Assign global classifications
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = Seurat:::MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data),
                                                 FUN = function(x) {
                                                   return(which(x = data[, x] == hash.max[x])[1])
                                                 })])
  hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data),
                                                    FUN = function(x) {
                                                      return(which(x = data[, x] == hash.second[x])[1])
                                                    })])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])),
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global ==
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global ==
                                                                           "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID,
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID",
                                                          "secondID", "margin", "classification", "classification.global"),
                                                 sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay,
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  object$hash.ID <- Idents(object = object)
  return(object)
}

# Input: count matrix of HTO counts, name for output plots
# Output: dataframe and pdf
deMULTIplex <- function(counts, plot.name) {
  bar.table <- as.data.frame(t(counts))
  bar.table.full <- bar.table
  gg.list <- list()

  ## Set iter
  iter <- 1

  ## Perform Quantile Sweep
  bar.table_sweep.list <- list()
  n <- 0
  for (q in seq(0.01, 0.99, by=0.02)) {
    n <- n + 1
    invisible(
      capture.output(
        bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
      )
    )
    names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
  }

  ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
  threshold.results <- findThresh(call.list=bar.table_sweep.list)
  gg.list[[iter]] <- ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) +
    geom_line() +
    theme(legend.position = "bottom") +
    labs(x = "quantile", y = "Proportion of cells", title = paste0("Round ", iter)) +
    geom_vline(xintercept=threshold.results$extrema, lty=2) +
    scale_color_manual(values=c("red","black","blue")) +
    coord_fixed()

  ## Finalize classifications, remove negative cells
  round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
  neg.cells <- names(round.calls)[which(round.calls == "Negative")]

  ## Check number of negative cells called, take them out of table
  num.neg.cells.in.round <- length(which(round.calls == "Negative"))
  if (num.neg.cells.in.round != 0) {
    message(paste0("Round ", iter, ": ", num.neg.cells.in.round, " Negative cells removed"))
    bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
  }else{
    message(paste0("Done! Took ", iter, " rounds to classify all cells."))
  }

  ## Repeat until all cells are classified
  while (num.neg.cells.in.round != 0) {
    ## Increase iter
    iter <- iter + 1

    ## Perform Quantile Sweep
    bar.table_sweep.list <- list()
    n <- 0
    for (q in c(seq(0.01, 0.99, by=0.02))) {
      n <- n + 1
      invisible(
        capture.output(
          bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
        )
      )
      names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
    }

    ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
    threshold.results <- findThresh(call.list=bar.table_sweep.list)
    gg.list[[iter]] <- ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) +
      geom_line() +
      theme(legend.position = "bottom") +
      labs(x = "quantile", y = "Proportion of cells", title = paste0("Round ", iter)) +
      geom_vline(xintercept=threshold.results$extrema, lty=2) +
      scale_color_manual(values=c("red","black","blue")) +
      coord_fixed()

    if(length(threshold.results$extrema) == 0) {
      message(paste0("Round ", iter, ": Stopping, as no threshold extrema could be found."))
      break
    }

    ## Finalize classifications, remove negative cells
    round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
    neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

    ## Check number of negative cells called, take them out of table
    num.neg.cells.in.round <- length(which(round.calls == "Negative"))
    if (num.neg.cells.in.round != 0) {
      message(paste0("Round ", iter, ": ", num.neg.cells.in.round, " Negative cells removed"))
      bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
    }else{
      message(paste0("Round ", iter, ": Done! All cells classified."))
    }
  }

  ## Finalize calls
  final.calls <- c(round.calls[which(round.calls != "Negative")], rep("Negative",length(neg.cells)))
  names(final.calls) <- c(names(round.calls[which(round.calls != "Negative")]), neg.cells)
  final.calls.rescued <- final.calls

  ## Make sure each tag has cells assigned to it, and rescue those tags that don't
  for (tag in colnames(bar.table.full)) {
    if (length(which(final.calls == tag)) == 0) {
      ## Perform Quantile Sweep
      bar.table_sweep.list <- list()
      n <- 0
      for (q in c(seq(0.01, 0.99, by=0.02))) {
        n <- n + 1
        invisible(
          capture.output(
            bar.table_sweep.list[[n]] <- classifyCells(bar.table.full[neg.cells,], q=q)
          )
        )
        names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
      }
      ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
      threshold.results <- findThresh(call.list=bar.table_sweep.list)
      round.calls <- classifyCells(bar.table.full[neg.cells,], q=findQ(threshold.results$res, threshold.results$extrema))
      ## Add graph
      gg.list[[length(gg.list) + 1]] <- ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) +
        geom_line() +
        theme(legend.position = "bottom") +
        labs(x = "quantile", y = "Proportion of cells", title = paste0("Rescue round for ", tag)) +
        geom_vline(xintercept=threshold.results$extrema, lty=2) +
        scale_color_manual(values=c("red","black","blue")) +
        coord_fixed()
      ## Rescue these negative cells
      final.calls.rescued[names(round.calls[which(round.calls == tag)])] <- tag
      message(paste0(tag, " cells rescued: ", length(which(round.calls == tag))))
    }else{
      message(paste0(tag, " cells: ", length(which(final.calls == tag))))
      next
    }
  }

  ## Perform semi-supervised negative cell reclassification
  reclass.cells <- findReclassCells(bar.table.full, names(final.calls.rescued)[which(final.calls.rescued=="Negative")])
  invisible(
    capture.output(
      reclass.res <- rescueCells(bar.table.full, final.calls.rescued, reclass.cells)
    )
  )
  cutoff <- tryCatch(
    expr = foverlaps(
      x = data.table(
        ClassStability = reclass.res$ClassStability[-1],
        start = reclass.res$MatchRate_mean[-1] - 3*reclass.res$MatchRate_sd[-1],
        end = reclass.res$MatchRate_mean[-1] + 3*reclass.res$MatchRate_sd[-1]
      ) %>% na.omit(),
      y = data.table(
        ClassStability = reclass.res$ClassStability[1],
        start = reclass.res$MatchRate_mean[1] - 3*reclass.res$MatchRate_sd[1],
        end = reclass.res$MatchRate_mean[1] + 3*reclass.res$MatchRate_sd[1]
      ) %>% setkey(start, end),
      type = "any",
      which = TRUE
    ) %>%
      .$yid %>%
      is.na() %>%
      which() %>%
      max() + 1,
    error = function (e) {
      # If this throws an error, set cutoff to class stability value with a mean match rate closest to the true match rate
      dt <- data.table(
        ClassStability = reclass.res$ClassStability[-1],
        MatchRate_mean = reclass.res$MatchRate_mean[-1]
      ) %>%
        setkey(MatchRate_mean)
      return(dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability)
    }
  )
  # If cutoff is -Inf (because all the ranges overlap and max returns -Inf), take the closest value to the true mean
  if (cutoff == -Inf) {
    dt <- data.table(
      ClassStability = reclass.res$ClassStability[-1],
      MatchRate_mean = reclass.res$MatchRate_mean[-1]
    ) %>%
      setkey(MatchRate_mean)
    cutoff <- dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability
  }
  # If cutoff is the max + 1 (because the ranges are so narrow that there are no overlaps), take the closest value to the true mean
  if (cutoff == max(reclass.res$ClassStability[-1])+1) {
    dt <- data.table(
      ClassStability = reclass.res$ClassStability[-1],
      MatchRate_mean = reclass.res$MatchRate_mean[-1]
    ) %>%
      setkey(MatchRate_mean)
    cutoff <- dt[J(reclass.res$MatchRate_mean[1]), roll = "nearest"]$ClassStability
  }

  ## Visualize results
  gg.list[[length(gg.list) + 1]] <- ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) +
    geom_point() +
    xlim(c(nrow(reclass.res)-1,1)) +
    ylim(c(0, reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1]+0.05)) +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_vline(xintercept = cutoff, col = "blue") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    labs(title = paste0("Class Stability Cutoff = ", cutoff))

  ## Confirm the cutoff (taken out for full automation)
  # print(gg.list[[length(gg.list)]])
  # cutoff <- as.integer(readline("Please confirm Class Stability Cutoff: "))

  ## Put new cutoff in place
  gg.list[[length(gg.list)]] <- ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) +
    geom_point() +
    xlim(c(nrow(reclass.res)-1,1)) +
    ylim(c(0, reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1]+0.05)) +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_vline(xintercept = cutoff, col = "blue") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
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

  df <- data.frame(
    deMULTIplex.calls = final.calls,
    deMULTIplex.calls.rescued = final.calls.rescued
  )

  return(df)
}

# Group small clusters into larger ones by setting a threshold
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
  return(ids)
}
