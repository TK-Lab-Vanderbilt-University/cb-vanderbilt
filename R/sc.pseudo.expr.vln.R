#' @title Comprehensive Pseudo-bulk Gene Expression Analyzer & Plotter
#'
#' @description
#' This function provides a robust, end-to-end pipeline for visualizing single-cell 
#' RNA-seq data using pseudo-bulk aggregation at the sample (patient) level. 
#' Treating patients as biological replicates resolves the false-positive p-value 
#' inflation commonly seen in single-cell differential expression analyses.
#' 
#' Key Features:
#' 1. Flexible Cell-Type Filtering: Analyze the global landscape ("all" cells) or 
#'    compare specific cell lineages side-by-side using faceting.
#' 2. Multiple Normalizations: Choose between 'Average RNA Expression' (raw counts 
#'    scaled by cell number), 'log2(CPM + 1)', or Seurat's 'Log Normalization'. 
#'    Set method = "all" to visualize all three side-by-side.
#' 3. Enhanced Visualization: Features Nature-style aesthetics, including customized 
#'    colors, individual patient jitter points, and a prominent median line with 
#'    clean margins inside each violin plot.
#' 4. Advanced Statistics: Seamlessly compute significance between any specific 
#'    pairs (or all combinations) using T-test, Wilcoxon, ANOVA, or Kruskal-Wallis.
#' 5. Publication Ready: Dynamic legend control and optimized layout formatting.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene String or vector. Gene name(s) to analyze.
#' @param group_var String. Metadata column for groups (default: "condition").
#' @param sample_var String. Metadata column for sample IDs (default: "orig.ident").
#' @param celltype_var String. Metadata column for cell types (default: "detailed.celltypes").
#' @param celltypes String or vector. "all" for global analysis, or a vector of specific cell types.
#' @param method Normalization: "all", "log.cpm" (default), "log.norm", or "raw.expression".
#' @param comparisons String or list. "all" for comprehensive pairwise stats, or specific pairs.
#' @param stats_method Statistical test (default: "t.test").
#' @param stats_label P-value format: "p.signif" or "p.format".
#' @param color_pal Custom color palette for groups.
#' @param dot_color String. Fill color for jitter points (default: "white").
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
                               dot_color = "white",
                               show_legend = FALSE) {
  
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggpubr)
  require(Seurat)
  require(patchwork)

  # Parameter handling
  if (length(method) > 1) method <- method[1] 
  stats_method <- match.arg(stats_method)
  stats_label <- match.arg(stats_label)

  # Check if all target genes exist
  missing_genes <- target_gene[!target_gene %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    stop(paste0("Genes not found in object: ", paste(missing_genes, collapse = ", ")))
  }

  # 1. Prepare Data & Logic for Cell Types
  meta_data <- seurat_obj@meta.data
  
  if (length(celltypes) == 1 && celltypes == "all") {
    valid_cells <- rownames(meta_data)
    split_mode_cell <- FALSE
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% celltypes, ])
    if (length(valid_cells) == 0) stop("No matching cell types found in the object.")
    split_mode_cell <- TRUE 
  }

  # 2. Extract Counts for multiple genes
  raw_counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[target_gene, valid_cells, drop = FALSE]
  
  # Reshape data to long format
  df_list <- lapply(target_gene, function(g) {
    data.frame(
      sample = meta_data[valid_cells, sample_var],
      group = meta_data[valid_cells, group_var],
      celltype = meta_data[valid_cells, celltype_var],
      gene = g,
      counts = as.numeric(raw_counts[g, ]),
      total_umi = seurat_obj$nCount_RNA[valid_cells]
    )
  })
  df <- do.call(rbind, df_list) %>% filter(!is.na(group))

  # 3. Pseudo-bulk Aggregation
  grouping_cols <- c("sample", "group", "gene")
  if (split_mode_cell) grouping_cols <- c(grouping_cols, "celltype")
  
  pb_df <- df %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      n_cells = n(), 
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
      data$Expression <- data$pb_gene_counts / data$n_cells
      y_lab <- "Average RNA Expression"
    } else if (method_type == "log.cpm") {
      data$Expression <- log2((data$pb_gene_counts / data$pb_total_umi) * 1e6 + 1)
      y_lab <- expression(bold(log[2] * (CPM + 1)))
    } else if (method_type == "log.norm") {
      data$Expression <- log1p((data$pb_gene_counts / data$pb_total_umi) * 10000)
      y_lab <- "Log Normalization" 
    }

    p <- ggplot(data, aes(x = group, y = Expression, fill = group)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.7, color = "black") +
      geom_jitter(width = 0.1, size = 1.5, shape = 21, fill = dot_color, color = "black") +
      stat_summary(fun = median, fun.min = median, fun.max = median, 
                   geom = "crossbar", width = 0.4, color = "black", linewidth = 0.6) +
      theme_classic() +
      labs(title = paste0("Method: ", method_type), y = y_lab, x = "") +
      theme(plot.title = element_text(face="bold", hjust=0.5, size=12),
            axis.text.x = element_text(angle=45, hjust=1, face="bold", size=10))

    # Faceting Logic for Genes and Celltypes
    if (length(target_gene) > 1 && split_mode_cell) {
      p <- p + facet_grid(gene ~ celltype, scales = "free_y")
    } else if (length(target_gene) > 1) {
      p <- p + facet_wrap(~gene, scales = "free_y")
    } else if (split_mode_cell) {
      p <- p + facet_wrap(~celltype, scales = "free_y")
    }

    if (!is.null(color_pal)) p <- p + scale_fill_manual(values = color_pal)
    if (!is.null(comparisons)) {
      p <- p + stat_compare_means(comparisons = comparisons, method = stats_method, label = stats_label)
    }
    if (!show_legend) p <- p + theme(legend.position = "none")
    
    return(p + theme(strip.background = element_blank(), strip.text = element_text(face="bold", size=11)))
  }

  # 6. Final Layout
  if (method == "all") {
    p_raw <- generate_plot(pb_df, "raw.expression")
    p_cpm <- generate_plot(pb_df, "log.cpm")
    p_norm <- generate_plot(pb_df, "log.norm")
    combined_plot <- p_raw | p_cpm | p_norm
    if (show_legend) combined_plot <- combined_plot + plot_layout(guides = "collect")
    return(combined_plot)
  } else {
    return(generate_plot(pb_df, method))
  }
}
