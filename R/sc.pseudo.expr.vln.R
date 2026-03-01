#' @title Comprehensive Pseudo-bulk Gene Expression Analyzer & Plotter
#'
#' @description
#' This function provides a robust pipeline for visualizing single-cell RNA-seq data 
#' using pseudo-bulk aggregation at the sample (patient) level. By treating patients 
#' as biological replicates, it resolves the false-positive p-value inflation common 
#' in cell-level differential expression analyses. 
#' 
#' Key Capabilities:
#' 1. Flexible Filtering: Analyze the global landscape ("all" cells) or compare specific 
#'    lineages side-by-side using facets (e.g., celltypes = c("Progenitor", "Monocyte")).
#' 2. Multiple Normalizations: Choose between 'Average RNA Expression' (raw counts scaled 
#'    by cell number), 'log2(CPM + 1)', or Seurat's default 'Log Normalization'. 
#'    Use method = "all" to visualize all three side-by-side.
#' 3. Advanced Statistics: Seamlessly compute significance between specific pairs 
#'    (or all pairs) using T-test, Wilcoxon, ANOVA, etc.
#' 4. Publication Ready: Customizable colors, dynamic legend control, and Nature-style layouts.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene String. Gene name to analyze.
#' @param group_var String. Metadata column for groups (default: "condition").
#' @param sample_var String. Metadata column for sample IDs (default: "orig.ident").
#' @param celltype_var String. Metadata column for cell types (default: "detailed.celltypes").
#' @param celltypes String or vector. "all" for global analysis, or a vector of specific cell types.
#' @param method Normalization: "all", "log.cpm" (default), "log.norm", or "raw.expression".
#' @param comparisons String or list. "all" for comprehensive pairwise stats, or specific pairs.
#' @param stats_method Statistical test (default: "t.test").
#' @param stats_label P-value format: "p.signif" or "p.format".
#' @param color_pal Custom color palette for groups.
#' @param show_legend Logical. Toggle legend display (default: FALSE).
#' 
#' @return A ggplot2 object, or a patchwork object if method = "all".
#' @author Hyundong Yoon
#' @export
sc.pseudo.expr.vln <- function(seurat_obj,
                               target_gene,
                               group_var = "condition",
                               sample_var = "orig.ident",
                               celltype_var = "detailed.celltypes",
                               celltypes = "all", 
                               method = c("log.cpm", "log.norm", "raw.expression", "all"),
                               comparisons = NULL,
                               stats_method = c("t.test", "wilcox.test", "anova", "kruskal.test"),
                               stats_label = c("p.signif", "p.format"),
                               color_pal = NULL,
                               show_legend = FALSE) {
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(Seurat)
  require(patchwork)

  # Parameter handling
  if (length(method) > 1) method <- method[1] 
  stats_method <- match.arg(stats_method)
  stats_label <- match.arg(stats_label)

  if (!target_gene %in% rownames(seurat_obj)) stop("Gene not found in Seurat object.")

  # 1. Prepare Data & Logic for Cell Types
  meta_data <- seurat_obj@meta.data
  
  if (length(celltypes) == 1 && celltypes == "all") {
    valid_cells <- rownames(meta_data)
    split_mode <- FALSE
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% celltypes, ])
    if (length(valid_cells) == 0) stop("No matching cell types found in the object.")
    split_mode <- TRUE 
  }

  # 2. Extract Counts
  df <- data.frame(
    sample = meta_data[valid_cells, sample_var],
    group = meta_data[valid_cells, group_var],
    celltype = meta_data[valid_cells, celltype_var],
    counts = as.numeric(GetAssayData(seurat_obj, layer = "counts")[target_gene, valid_cells]),
    total_umi = seurat_obj$nCount_RNA[valid_cells]
  ) %>% filter(!is.na(group))

  # 3. Pseudo-bulk Aggregation (Adding cell count tracking for raw averages)
  grouping_cols <- if(split_mode) c("sample", "group", "celltype") else c("sample", "group")
  
  pb_df <- df %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      n_cells = n(), # Track number of cells to calculate average expression
      .groups = "drop"
    )

  # 4. Handle Comparisons ("all" logic)
  if (!is.null(comparisons) && length(comparisons) == 1 && comparisons == "all") {
    group_levels <- as.character(unique(pb_df$group))
    if (length(group_levels) > 1) {
      comparisons <- combn(group_levels, 2, simplify = FALSE)
    } else {
      comparisons <- NULL 
    }
  }

  # 5. Helper function to generate a single plot
  generate_plot <- function(data, method_type) {
    
    if (method_type == "raw.expression") {
      # [Updated] Average expression per cell instead of massive sum
      data$Expression <- data$pb_gene_counts / data$n_cells
      y_lab <- "Average RNA Expression"
      
    } else if (method_type == "log.cpm") {
      data$Expression <- log2((data$pb_gene_counts / data$pb_total_umi) * 1e6 + 1)
      y_lab <- expression(bold(log[2] * (CPM + 1)))
      
    } else if (method_type == "log.norm") {
      data$Expression <- log1p((data$pb_gene_counts / data$pb_total_umi) * 10000)
      y_lab <- "Log Normalization" # [Updated] Simplified label
    }

    p <- ggplot(data, aes(x = group, y = Expression, fill = group)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = "black") +
      geom_jitter(width = 0.1, size = 2, shape = 21, fill = "white", color = "black") +
      theme_classic() +
      labs(title = paste0(target_gene, " (", method_type, ")"), y = y_lab, x = "") +
      theme(plot.title = element_text(face="bold", hjust=0.5, size=13),
            axis.text.x = element_text(angle=45, hjust=1, face="bold", size=11))

    if (split_mode) {
      p <- p + facet_wrap(~celltype, scales = "free_y") +
        theme(strip.background = element_blank(), strip.text = element_text(face="bold", size=12))
    }

    if (!is.null(color_pal)) p <- p + scale_fill_manual(values = color_pal)
    
    if (!is.null(comparisons)) {
      p <- p + stat_compare_means(comparisons = comparisons, method = stats_method, label = stats_label)
    }
    
    if (!show_legend) {
      p <- p + theme(legend.position = "none")
    }
    
    return(p)
  }

  # 6. Logic for multiple methods vs single method
  if (method == "all") {
    p_raw <- generate_plot(pb_df, "raw.expression")
    p_cpm <- generate_plot(pb_df, "log.cpm")
    p_norm <- generate_plot(pb_df, "log.norm")
    
    combined_plot <- p_raw | p_cpm | p_norm
    
    if (show_legend) {
      combined_plot <- combined_plot + plot_layout(guides = "collect")
    }
    
    return(combined_plot)
  } else {
    return(generate_plot(pb_df, method))
  }
}
