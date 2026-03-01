#' @title Pseudo-bulk Gene Expression Visualization (Violin Plot)
#'
#' @description
#' This function aggregates single-cell RNA-seq expression data into pseudo-bulk 
#' profiles at the sample (patient) level. By treating each sample as a biological 
#' replicate (N), this method provides a more robust statistical framework for 
#' comparing groups (e.g., Healthy vs. MDS vs. sAML) and reduces the inflated 
#' p-values often seen in cell-level analyses. 
#' 
#' Features:
#' 1. Supports cell-type specific filtering.
#' 2. Offers multiple normalization methods: Raw, Log-CPM, and Seurat's Log-Normalization.
#' 3. Integrates customizable statistical tests (T-test, Wilcoxon, etc.).
#' 4. Customizable color palettes for publication-quality figures.
#'
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param target_gene String. The name of the gene to be visualized.
#' @param group_var String. Metadata column name for grouping (default: "condition").
#' @param sample_var String. Metadata column name for sample/patient IDs (default: "orig.ident").
#' @param celltype_var String. Metadata column name for cell type annotations (default: "detailed.celltypes").
#' @param filter_celltype String vector or NULL. List of cell types to include. If NULL, all cells are used.
#' @param method String. Normalization method: "log.cpm" (default), "log.norm", or "raw.expression".
#' @param comparisons List of character vectors. Group pairs for statistical testing (e.g., list(c("MDS", "sAML"))).
#' @param stats_method String. Statistical test: "t.test", "wilcox.test", "anova", or "kruskal.test".
#' @param stats_label String. P-value format: "p.signif" (stars) or "p.format" (numeric value).
#' @param color_pal Named character vector. Custom colors for the groups.
#' 
#' @return A ggplot2 object.
#' 
#' @author Hyundong Yoon
#' @export
#'
#' @import Seurat dplyr ggplot2 ggpubr
sc.pseudo.expr.vln <- function(seurat_obj,
                               target_gene,
                               group_var = "condition",
                               sample_var = "orig.ident",
                               celltype_var = "detailed.celltypes",
                               filter_celltype = NULL,
                               method = c("log.cpm", "log.norm", "raw.expression"),
                               comparisons = NULL,
                               stats_method = c("t.test", "wilcox.test", "anova", "kruskal.test"),
                               stats_label = c("p.signif", "p.format"),
                               color_pal = NULL) {
  
  # 1. Argument Validation
  method <- match.arg(method)
  stats_method <- match.arg(stats_method)
  stats_label <- match.arg(stats_label)
  
  if (!target_gene %in% rownames(seurat_obj)) {
    stop(paste0("Error: Gene '", target_gene, "' was not found in the Seurat object."))
  }
  
  # 2. Data Extraction and Cell Filtering
  meta_data <- seurat_obj@meta.data
  
  if (!is.null(filter_celltype)) {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% filter_celltype, ])
    if (length(valid_cells) == 0) stop("No cells found matching the specified 'filter_celltype'.")
  } else {
    valid_cells <- rownames(meta_data)
  }
  
  # Prepare dataframe for aggregation
  raw_counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
  
  df <- data.frame(
    sample = meta_data[valid_cells, sample_var],
    group = meta_data[valid_cells, group_var],
    counts = as.numeric(raw_counts[target_gene, valid_cells]),
    total_umi = seurat_obj$nCount_RNA[valid_cells]
  )
  
  # 3. Pseudo-bulk Aggregation (Summation by sample)
  pb_df <- df %>%
    group_by(sample, group) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      .groups = "drop"
    ) %>%
    filter(!is.na(group))
  
  # 4. Normalization Logic
  if (method == "raw.expression") {
    pb_df$Expression <- pb_df$pb_gene_counts
    y_label <- "Total Sum of Raw Counts"
  } else if (method == "log.cpm") {
    # Calculation: log2((counts / total) * 1,000,000 + 1)
    pb_df$Expression <- log2((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 1e6 + 1)
    y_label <- expression(bold(log[2] * (CPM + 1)))
  } else if (method == "log.norm") {
    # Calculation: ln((counts / total) * 10,000 + 1)
    pb_df$Expression <- log1p((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 10000)
    y_label <- "ln(CP10k + 1)"
  }
  
  # 5. Visualization
  p <- ggplot(pb_df, aes(x = group, y = Expression, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 3, shape = 21, fill = "white", color = "black") +
    theme_classic() +
    labs(
      title = paste0("Pseudo-bulk: ", target_gene),
      subtitle = paste0("Normalization: ", method, 
                        if (!is.null(filter_celltype)) " (Filtered Cells)" else " (All Cells)"),
      y = y_label, x = ""
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkblue"),
      axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
      axis.title.y = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  # Apply custom color palette
  if (!is.null(color_pal)) {
    p <- p + scale_fill_manual(values = color_pal)
  }
  
  # Add statistical layers (Supports multiple pairs)
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(
      comparisons = comparisons,
      method = stats_method,
      label = stats_label
    )
  }
  
  return(p)
}
