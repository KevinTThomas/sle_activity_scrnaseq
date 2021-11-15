Figure 5
================
Kevin Thomas
15 November, 2021

## 5A

``` r
mono_clusters <- sc@meta.data %>%
  filter(coarse_cell_type == "Monocytes") %>%
  pull(clusters_annotated) %>%
  unique()
DimPlot(
  object = sc,
  cells = WhichCells(object = sc, expression = `coarse_cell_type` %in% c("Monocytes") & `wnnUMAP_1` < 0 & `wnnUMAP_2` < 3),
  reduction = "umap_wnn",
  dims = c(1,2),
  group.by = "clusters_annotated",
  cols = `names<-`(pals::cols25(n = length(levels(sc$clusters_annotated))), levels(sc$clusters_annotated))[levels(sc$clusters_annotated) %in% mono_clusters]
) +
  labs(title = NULL, x = NULL, y = NULL) +
  coord_fixed()
```

![](figure_5_s19_files/figure-gfm/5A-1.png)<!-- -->

## 5B

``` r
# Generate CITEseq markers for Monocytes
mono_markers_cite <-
  ## Look only at Monocytes
  subset(sc, `coarse_cell_type` %in% c("Monocytes")) %>%
  ## Define CITEseq markers of each annotated cluster by fast wilcoxon test
  wilcoxauc(
    X = .,
    seurat_assay = "CITE",
    group_by = "clusters_annotated"
  ) %>%
  as_tibble() %>%
  dplyr::select(feature,
                cluster = group,
                p_val = pval,
                p_val_adj = padj,
                avg_logFC = logFC,
                pct.1 = pct_in,
                pct.2 = pct_out,
                auc
  ) %>%
  ## Significance filter
  filter(p_val_adj <= 0.1) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster)[c(5, 10, 7, 11, 9, 6, 1, 2, 8, 4, 3)])) %>%
  group_by(cluster) %>%
  ## Positive markers only, remove the isotype control
  filter(avg_logFC > 0, feature != "isotype-control") %>%
  ## Rank markers in each cluster by average logFC
  arrange(-avg_logFC, by.group = TRUE) %>%
  ## If two clusters share a marker, the cluster with a higher average logFC is prioritized
  filter(
    avg_logFC == map_dbl(split(., .$feature), ~max(.x$avg_logFC))[feature]
  ) %>%
  ## Top 3 features for each cluster
  slice_head(n = 3) %>%
  .$feature
# Bubbleplot
bubbleplot(
  object = subset(sc, `coarse_cell_type` %in% c("Monocytes")),
  assay = "CITE",
  slot = "data",
  features_plot = c("CD14", mono_markers_cite),
  preserve_feature_order = TRUE,
  grouping_var = "clusters_annotated",
  filter_exp_pct_thresh = 1,
  avg_func = "median",
  do_return = TRUE
) +
  scale_color_distiller(palette = "Reds", direction = 1) +
  labs(title = NULL, x = NULL, y = NULL, col = "Scaled median expression", size = "Percent cells expressing") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour', which will replace the existing scale.

![](figure_5_s19_files/figure-gfm/5B-1.png)<!-- -->

## 5C

``` r
# Generate gene markers for Monocytes
mono_markers_genes <-
  subset(sc, `coarse_cell_type` %in% c("Monocytes")) %>%
  wilcoxauc(
    X = .,
    seurat_assay = "SCT",
    group_by = "clusters_annotated"
  ) %>%
  as_tibble() %>%
  dplyr::select(feature,
                cluster = group,
                p_val = pval,
                p_val_adj = padj,
                avg_logFC = logFC,
                pct.1 = pct_in,
                pct.2 = pct_out,
                auc
  ) %>%
  filter(p_val_adj <= 0.1) %>%
  mutate(cluster = factor(cluster, levels = unique(cluster)[c(5, 10, 7, 11, 9, 6, 1, 2, 8, 4, 3)])) %>%
  group_by(cluster) %>% 
  filter(auc > 0.6) %>%
  arrange(-auc, by.group = TRUE) %>%
  filter(
    auc == map_dbl(split(., .$feature), ~max(.x$auc))[feature]
  ) %>%
  slice_head(n = 5) %>%
  .$feature
# Heatmap
scale <- function(x) {
  (x-min(x))/(max(x)-min(x))
}
mat <- FetchData(
  object = subset(sc, `coarse_cell_type` %in% c("Monocytes")),
  vars = c("clusters_annotated", rev(mono_markers_genes))
) %>%
  as_tibble(rownames = "cell") %>%
  group_by(clusters_annotated) %>%
  summarize_if(is.numeric, mean) %>%
  mutate_at(mono_markers_genes, scale) %>%
  as.data.frame() %>%
  column_to_rownames("clusters_annotated")
pheatmap(
  mat = t(mat),
  scale = "none",
  breaks = seq.int(from = 0, to = 1, length.out = 100),
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Blues"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 90
)
```

![](figure_5_s19_files/figure-gfm/5C-1.png)<!-- -->

## 5D-5G, S19A-S19G

``` r
# Collect percent of Monocytes in each cluster per patient
vln_data <- sc %>%
  FetchData(
    object = .,
    vars = c("clusters_annotated", "subject_id", "ancestry", "classification"),
    cells = WhichCells(sc, expression = `coarse_cell_type` %in% c("Monocytes"))
  ) %>%
  group_by(ancestry, classification, subject_id, clusters_annotated) %>%
  summarize(total_cells = n()) %>%
  mutate(
    percent_total = 100*total_cells/sum(total_cells),
    group = paste(ancestry, classification)
  )
```

    ## `summarise()` has grouped output by 'ancestry', 'classification', 'subject_id'. You can override using the `.groups` argument.

``` r
# Re-order factor levels
vln_data$group <- factor(
  x = vln_data$group,
  levels = c(
    "EA Control",
    "EA SLE INACT",
    "EA SLE ACT",
    "AA Control",
    "AA SLE INACT",
    "AA SLE ACT"
  )
)
# Violin plots for each cluster
vln_data %>%
  ggplot(aes(x = group, y = percent_total)) +
  geom_violin(aes(fill = clusters_annotated), scale = "width") +
  scale_fill_manual(
    values = `names<-`(
      x = pals::cols25(n = length(levels(vln_data$clusters_annotated))),
      values = levels(vln_data$clusters_annotated)
    )
  ) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.4) +
  ## Show statistics for Wilcoxon comparisons
  ggpubr::stat_compare_means(
    aes(label = ..p.format..),
    method = "wilcox.test",
    method.args = list(
      formula = percent_total ~ group,
      p.adjust.method = "BH"
    ),
    comparisons = list(
      c("EA Control", "AA Control"),
      c("EA Control", "EA SLE INACT"),
      c("EA SLE INACT", "EA SLE ACT"),
      c("EA Control", "EA SLE ACT"),
      c("AA Control", "AA SLE INACT"),
      c("AA SLE INACT", "AA SLE ACT"),
      c("AA Control", "AA SLE ACT")
    ),
    hide.ns = TRUE,
    show.legend = TRUE,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  ) +
  xlab(NULL) +
  ylab(NULL) +
  facet_wrap(~clusters_annotated, scales = "free", ncol = 6) +
  cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
```

    ## Warning in wilcox.test.default(c(5.26315789473684, 2.42261103633917, 13.3879781420765, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(6.41803989592368, 8, 7.69230769230769, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(6.41803989592368, 8, 7.69230769230769, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(4.42324371205551, 21.6, 15.3846153846154, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1.99479618386817, 2.4, 1.92307692307692, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1.561144839549, 11.2, 21.6346153846154, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(4.96031746031746, 2.50696378830084, 1.93236714975845, : cannot compute exact p-value with ties

![](figure_5_s19_files/figure-gfm/5D-5G-1.png)<!-- -->

## 5H

``` r
monos_scored <- subset(
  sc_scored,
  `clusters_annotated` %in% c("Classical monos", "Intermediate monos", "Non-classical monos")
)
grid <- expand.grid(levels(monos_scored$classification), str_sort(unique(monos_scored$clusters_annotated)))
monos_scored$cluster_class <- factor(
  paste(monos_scored$clusters_annotated, monos_scored$classification),
  levels = paste(grid$Var2, grid$Var1)
)
cowplot::plot_grid(
  bubbleplot(
  object = subset(monos_scored, `ancestry` == "EA"),
  features_plot = c("M1.2", "M3.4", "M5.12"),
  preserve_feature_order = TRUE,
  grouping_var = "cluster_class",
  colors_use = "magma",
) +
  scale_size_continuous(breaks = c(20, 40, 60, 80)) +
  labs(title = "EA", y = NULL, x = "Interferon Response Modules"),
bubbleplot(
  object = subset(monos_scored, `ancestry` == "AA"),
  features_plot = c("M1.2", "M3.4", "M5.12"),
  preserve_feature_order = TRUE,
  grouping_var = "cluster_class",
  colors_use = "magma",
) +
  scale_size_continuous(breaks = c(20, 40, 60)) +
  labs(title = "AA", y = NULL, x = "Interferon Response Modules")
)
```

    ## Scale for 'size' is already present. Adding another scale for 'size', which will replace the existing scale.
    ## Scale for 'size' is already present. Adding another scale for 'size', which will replace the existing scale.

![](figure_5_s19_files/figure-gfm/5H-1.png)<!-- -->

## 5I-5K

``` r
# DC Violin plots
# Collect percent of Dendritic cells in each cluster per patient
vln_data <- sc %>%
  FetchData(
    object = .,
    vars = c("clusters_annotated", "subject_id", "ancestry", "classification"),
    cells = WhichCells(sc, expression = `coarse_cell_type` %in% c("cDCs", "pDCs"))
  ) %>%
  group_by(ancestry, classification, subject_id, clusters_annotated) %>%
  summarize(total_cells = n()) %>%
  mutate(
    percent_total = 100*total_cells/sum(total_cells),
    group = paste(ancestry, classification)
  )
```

    ## `summarise()` has grouped output by 'ancestry', 'classification', 'subject_id'. You can override using the `.groups` argument.

``` r
# Re-order factor levels
vln_data$group <- factor(
  x = vln_data$group,
  levels = c(
    "EA Control",
    "EA SLE INACT",
    "EA SLE ACT",
    "AA Control",
    "AA SLE INACT",
    "AA SLE ACT"
  )
)
# Violin plots for each cluster
vln_data %>%
  ggplot(aes(x = group, y = percent_total)) +
  geom_violin(aes(fill = clusters_annotated), scale = "width") +
  scale_fill_manual(
    values = `names<-`(
      x = pals::cols25(n = length(levels(vln_data$clusters_annotated))),
      values = levels(vln_data$clusters_annotated)
    )
  ) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.4) +
  ## Show statistics for Wilcoxon comparisons
  ggpubr::stat_compare_means(
    aes(label = ..p.format..),
    method = "wilcox.test",
    method.args = list(
      formula = percent_total ~ group,
      p.adjust.method = "BH"
    ),
    comparisons = list(
      c("EA Control", "AA Control"),
      c("EA Control", "EA SLE INACT"),
      c("EA SLE INACT", "EA SLE ACT"),
      c("EA Control", "EA SLE ACT"),
      c("AA Control", "AA SLE INACT"),
      c("AA SLE INACT", "AA SLE ACT"),
      c("AA Control", "AA SLE ACT")
    ),
    hide.ns = TRUE,
    show.legend = TRUE,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  ) +
  xlab(NULL) +
  ylab(NULL) +
  facet_wrap(~clusters_annotated, scales = "free", ncol = 6) +
  cowplot::theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
```

    ## Warning in wilcox.test.default(c(1.47058823529412, 5.26315789473684, 12.5, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(1.47058823529412, 5.26315789473684, 12.5, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5, 5.88235294117647, 6.66666666666667, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(2.94117647058824, 13.3333333333333, 5.26315789473684, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(2.94117647058824, 13.3333333333333, 5.26315789473684, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(87.7551020408163, 85.7142857142857, 85.4545454545455, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(95.5882352941177, 100, 86.6666666666667, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(95.5882352941177, 100, 86.6666666666667, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(90, 93.3333333333333, 88.2352941176471, : cannot compute exact p-value with ties

![](figure_5_s19_files/figure-gfm/5I-5K-1.png)<!-- -->

## 5L

``` r
# DC HLA Plots
hla_genes <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-E", "HLA-F")
data_to_plot <- FetchData(
  object = sc,
  vars = c("classification", "ancestry", hla_genes),
  cells =  WhichCells(sc, expression = `coarse_cell_type` == "cDCs")
) %>%
  as_tibble(rownames = "cell") %>%
  pivot_longer(cols = hla_genes, names_to = "feature", values_to = "count") %>%
  mutate(feature = factor(feature, levels = rev(hla_genes))) %>%
  group_by(
    ident = factor(
      x = paste0(get("ancestry"), " ", get("classification")),
      levels = c("EA Control", "EA SLE INACT", "EA SLE ACT", "AA Control", "AA SLE INACT", "AA SLE ACT")
    ),
    feature
  ) %>%
  summarize(
    avg_exp = mean(expm1(x = count)),
    pct_exp = 100*sum(count>0)/length(count)
  ) %>%
  ungroup() %>%
  group_by(feature) %>%
  mutate(avg_exp_scale = compositions::normalize(avg_exp))
```

    ## Note: Using an external vector in selections is ambiguous.
    ## ℹ Use `all_of(hla_genes)` instead of `hla_genes` to silence this message.
    ## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This message is displayed once per session.

    ## `summarise()` has grouped output by 'ident'. You can override using the `.groups` argument.

``` r
ggplot(
  data = data_to_plot,
  aes(
    x = ident,
    y = feature,
    size = pct_exp,
    color = avg_exp_scale
  )
) +
  geom_point() +
  geom_rect(
    aes(xmin = 0.5, xmax = 6.5, ymin = 10.5, ymax = 14.55),
    size = 0.5,
    col = "black",
    fill = "#00000000"
  ) +
  geom_rect(
    aes(xmin = 0.5, xmax = 6.5, ymin = 2.5, ymax = 10.45),
    size = 0.5,
    col = "red",
    fill = "#00000000"
  ) +
  scale_color_gradientn(
    colors = viridisLite::inferno(n = 100)
  ) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, size = "Percent cells expressing", col = "Scaled average expression") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
```

![](figure_5_s19_files/figure-gfm/5L-1.png)<!-- -->
