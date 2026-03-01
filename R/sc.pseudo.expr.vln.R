#' @title Flexible Pseudo-bulk Gene Expression Visualization
#'
#' @description
#' This function visualizes gene expression using pseudo-bulk aggregation. 
#' It allows for "Global" analysis (all cells) or "Side-by-side" comparison of 
#' multiple specific cell types using facets.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene String. Gene name.
#' @param group_var String. Grouping variable (default: "condition").
#' @param sample_var String. Sample/Patient ID variable (default: "orig.ident").
#' @param celltype_var String. Cell type annotation variable (default: "detailed.celltypes").
#' @param filter_celltype String or vector. Use "all" for global aggregation, 
#' or a vector of cell types (e.g., c("Progenitor", "GMP")) for specific analysis.
#' @param method Normalization: "log.cpm" (default), "log.norm", "raw.expression".
#' @param comparisons List of pairs for stats (e.g., list(c("MDS", "sAML"))).
#' @param stats_method Statistical test (default: "t.test").
#' @param stats_label P-value format: "p.signif" or "p.format".
#' @param color_pal Custom color palette for groups.
#' 
#' @return A ggplot2 object with optional faceting.
#' @author Hyundong Yoon
#' @export
sc.pseudo.expr.vln <- function(seurat_obj,
                               target_gene,
                               group_var = "condition",
                               sample_var = "orig.ident",
                               celltype_var = "detailed.celltypes",
                               filter_celltype = "all", 
                               method = c("log.cpm", "log.norm", "raw.expression"),
                               comparisons = NULL,
                               stats_method = c("t.test", "wilcox.test", "anova", "kruskal.test"),
                               stats_label = c("p.signif", "p.format"),
                               color_pal = NULL) {
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(Seurat)

  method <- match.arg(method)
  stats_method <- match.arg(stats_method)
  stats_label <- match.arg(stats_label)

  if (!target_gene %in% rownames(seurat_obj)) stop("Gene not found.")

  # 1. Prepare Data & Logic for Multi-celltype vs Global
  meta_data <- seurat_obj@meta.data
  
  if (length(filter_celltype) == 1 && filter_celltype == "all") {
    valid_cells <- rownames(meta_data)
    split_mode <- FALSE
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% filter_celltype, ])
    if (length(valid_cells) == 0) stop("No matching cell types found.")
    split_mode <- TRUE # Activate faceting
  }

  # 2. Extract Counts
  df <- data.frame(
    sample = meta_data[valid_cells, sample_var],
    group = meta_data[valid_cells, group_var],
    celltype = meta_data[valid_cells, celltype_var],
    counts = as.numeric(GetAssayData(seurat_obj, layer = "counts")[target_gene, valid_cells]),
    total_umi = seurat_obj$nCount_RNA[valid_cells]
  )

  # 3. Pseudo-bulk Aggregation
  # If split_mode is TRUE, we group by celltype as well
  grouping_cols <- if(split_mode) c("sample", "group", "celltype") else c("sample", "group")
  
  pb_df <- df %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      .groups = "drop"
    )

  # 4. Normalization
  if (method == "raw.expression") {
    pb_df$Expression <- pb_df$pb_gene_counts
    y_lab <- "Raw Counts"
  } else if (method == "log.cpm") {
    pb_df$Expression <- log2((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 1e6 + 1)
    y_lab <- expression(bold(log[2] * (CPM + 1)))
  } else if (method == "log.norm") {
    pb_df$Expression <- log1p((pb_df$pb_gene_counts / pb_df$pb_total_umi) * 10000)
    y_lab <- "ln(CP10k + 1)"
  }

  # 5. Plotting
  p <- ggplot(pb_df, aes(x = group, y = Expression, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 2, shape = 21, fill = "white") +
    theme_classic() +
    labs(title = paste0("Gene: ", target_gene), y = y_lab, x = "") +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          axis.text.x = element_text(angle=45, hjust=1, face="bold"))

  # Apply Faceting if multiple cell types are selected
  if (split_mode) {
    p <- p + facet_wrap(~celltype, scales = "free_y") +
      theme(strip.background = element_blank(), strip.text = element_text(face="bold", size=12))
  }

  if (!is.null(color_pal)) p <- p + scale_fill_manual(values = color_pal)
  
  if (!is.null(comparisons)) {
    p <- p + stat_compare_means(comparisons = comparisons, method = stats_method, label = stats_label)
  }

  return(p)
}
