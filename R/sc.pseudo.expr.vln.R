#' Pseudo-bulk Violin Plot for Single-Cell RNA-seq Data
#'
#' @description
#' Aggregates single-cell expression data into pseudo-bulk profiles at the sample (patient) level.
#' This approach effectively reduces false-positive rates in differential expression analysis
#' by treating each patient as a biological replicate (N). It supports cell-type specific
#' filtering and various normalization methods for robust visualization.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene Character. The gene to visualize (e.g., "LAIR1").
#' @param group_var Character. The metadata column name for group conditions. Default is "condition".
#' @param sample_var Character. The metadata column name for patient/sample IDs. Default is "orig.ident".
#' @param celltype_var Character. The metadata column name for cell types. Default is "detailed.celltypes".
#' @param filter_celltype Character vector or NULL. Specific cell types to filter and aggregate. If NULL, uses all cells. Default is NULL.
#' @param method Character. Normalization method: "log.cpm" (default), "log.norm", or "raw.expression".
#' @param comparisons List of character vectors. Group pairs for t-test. e.g., list(c("MDS", "sAML")).
#' @param color_pal Named vector of colors for the groups.
#' 
#' @return A ggplot object.
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
                                 color_pal = NULL) {
  
  # 1. Argument validation
  method <- match.arg(method)
  if (!target_gene %in% rownames(seurat_obj)) stop("Target gene not found in the Seurat object.")
  
  # 2. Extract metadata and raw counts
  meta_data <- seurat_obj@meta.data
  
  # Cell filtering logic
  if (!is.null(filter_celltype)) {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% filter_celltype, ])
    if (length(valid_cells) == 0) stop("No cells found for the specified filter_celltype.")
  } else {
    valid_cells <- rownames(meta_data) # Use all cells if NULL
  }
  
  # Create dataframe for pseudo-bulk
  df <- data.frame(
    sample = meta_data[valid_cells, sample_var],
    group = meta_data[valid_cells, group_var],
    counts = as.numeric(GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[target_gene, valid_cells]),
    total_umi = seurat_obj$nCount_RNA[valid_cells]
  )
  
  # 3. Pseudo-bulk Aggregation (Summing at the sample level)
  pb_df <- df %>%
    group_by(sample, group) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      .groups = "drop"
    ) %>%
    filter(!is.na(group)) # Remove NA groups
  
  # 4. Apply selected normalization method
  if (method == "raw.expression") {
    pb_df$Expression <- pb_df$pb_gene_counts
    y_label <- "Total Sum of Raw Counts"
    sub_title <- "Method: Raw Expression (No Normalization)"
    
  } else if (method == "log.cpm") {
    pb_df$Expression <- log2((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 1e6 + 1)
    y_label <- expression(bold(log[2] * (CPM + 1)))
    sub_title <- "Method: Standardized log2(CPM + 1)"
    
  } else if (method == "log.norm") {
    pb_df$Expression <- log1p((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 10000)
    y_label <- "ln(CP10k + 1)"
    sub_title <- "Method: Seurat Default LogNormalization"
  }
  
  # 5. Visualization
  p <- ggplot(pb_df, aes(x = group, y = Expression, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 2.5, shape = 21, fill = "white", color = "black") +
    theme_classic() +
    labs(
      title = paste0("Pseudo-bulk Profile: ", target_gene),
      subtitle = paste0(sub_title, ifelse(is.null(filter_celltype), " (All Cells)", paste0(" (", paste(filter_celltype, collapse=", "), ")"))),
      y = y_label,
      x = ""
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
      axis.title.y = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  # Add specific colors if provided
  if (!is.null(color_pal)) {
    p <- p + scale_fill_manual(values = color_pal)
  }
  
  # Add statistical comparisons if provided
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")
  }
  
  return(p)
}