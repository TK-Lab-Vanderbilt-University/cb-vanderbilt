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
#' @param celltypes String or vector. Use "all" for global aggregation, 
#' or a vector of cell types (e.g., c("Progenitor", "GMP")) for specific analysis.
#' @param method Normalization: "all", "log.cpm" (default), "log.norm", or "raw.expression".
#' @param comparisons String or list. Use "all" for all pairwise comparisons, or a list of pairs (e.g., list(c("MDS", "sAML"))).
#' @param stats_method Statistical test (default: "t.test").
#' @param stats_label P-value format: "p.signif" or "p.format".
#' @param color_pal Custom color palette for groups.
#' 
#' @return A ggplot2 object or a patchwork object if multiple methods are selected.
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
                               color_pal = NULL) {
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(Seurat)
  require(patchwork)

  # Parameter handling
  if (length(method) > 1) method <- method[1] # Default if not explicitly provided
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
  ) %>% filter(!is.na(group)) # Ensure no NA groups

  # 3. Pseudo-bulk Aggregation
  grouping_cols <- if(split_mode) c("sample", "group", "celltype") else c("sample", "group")
  
  pb_df <- df %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      pb_gene_counts = sum(counts),
      pb_total_umi = sum(total_umi),
      .groups = "drop"
    )

  # 4. Handle Comparisons ("all" logic)
  if (!is.null(comparisons) && length(comparisons) == 1 && comparisons == "all") {
    # Generate all pairwise combinations of the groups
    group_levels <- as.character(unique(pb_df$group))
    if (length(group_levels) > 1) {
      comparisons <- combn(group_levels, 2, simplify = FALSE)
    } else {
      comparisons <- NULL # Cannot compare if only 1 group exists
    }
  }

  # 5. Helper function to generate a single plot
  generate_plot <- function(data, method_type) {
    
    if (method_type == "raw.expression") {
      data$Expression <- data$pb_gene_counts
      y_lab <- "Raw Counts"
    } else if (method_type == "log.cpm") {
      data$Expression <- log2((data$pb_gene_counts / data$pb_total_umi) * 1e6 + 1)
      y_lab <- expression(bold(log[2] * (CPM + 1)))
    } else if (method_type == "log.norm") {
      data$Expression <- log1p((data$pb_gene_counts / data$pb_total_umi) * 10000)
      y_lab <- "ln(CP10k + 1)"
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
    return(p)
  }

  # 6. Logic for multiple methods vs single method
  if (method == "all") {
    p_raw <- generate_plot(pb_df, "raw.expression")
    p_cpm <- generate_plot(pb_df, "log.cpm")
    p_norm <- generate_plot(pb_df, "log.norm")
    
    # Return patchwork combined plot
    return(p_raw | p_cpm | p_norm)
  } else {
    # Return single plot
    return(generate_plot(pb_df, method))
  }
}
