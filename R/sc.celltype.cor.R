#' @title Systemic Cell-Type Co-expression Correlation Analysis
#'
#' @description
#' This function analyzes whether the expression of a target gene is systemically 
#' coordinated between different cell types across a cohort of samples. It calculates 
#' the percentage of gene-positive cells per sample/cell-type and performs 
#' Pearson correlation analysis.
#'
#' @param seurat_obj A Seurat object. Ensure the default assay is set to "RNA".
#' @param target_gene Character. The gene of interest (e.g., "LAIR1", "VSIR").
#' @param A.cell Character. The first cell type for the scatter plot (Y-axis).
#' @param B.cell Character. The second cell type for the scatter plot (X-axis).
#' @param celltype_col Character. Metadata column containing cell type labels. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column containing experimental groups (e.g., "condition").
#' @param sample_col Character. Metadata column containing sample IDs (e.g., "orig.ident").
#' @param min_cells Integer. Minimum number of cells required per sample/cell-type to be included. Default: 0.
#' @param min_samples Integer. Minimum number of biological replicates (samples) required to compute correlation. Prevents overfitting. Default: 4.
#' @param only.significant Logical. If TRUE, performs an all-to-all cell type pair screening. Default: FALSE.
#' @param color.use Named vector. Custom colors for the experimental groups.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cor_table}: A data frame of Pearson correlation results (R, P-value, N).
#'   \item \code{plot}: A ggplot2 object showing the linear regression between A.cell and B.cell.
#' }
#' @export
sc.celltype.cor <- function(seurat_obj,
                            target_gene,
                            A.cell = NULL,
                            B.cell = NULL,
                            celltype_col = "detailed.celltypes",
                            group_col = "condition",
                            sample_col = "orig.ident",
                            min_cells = 0,
                            min_samples = 4,
                            only.significant = FALSE,
                            color.use = NULL) {
  
  # Load dependencies
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggpubr)
    library(rlang)
  })
  
  cat(sprintf("[1/4] Extracting expression data for: %s\n", target_gene))
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # Fetch necessary data from Seurat object
  vars_to_fetch <- c(sample_col, group_col, celltype_col, target_gene)
  data.df <- Seurat::FetchData(seurat_obj, vars = vars_to_fetch)
  colnames(data.df) <- c("sample_id", "condition", "celltype", "expression")
  
  # Clean data
  data.df <- data.df %>% filter(!is.na(celltype), !is.na(condition))
  valid_cells <- setdiff(unique(as.character(data.df$celltype)), NA)
  
  cat("[2/4] Calculating gene-positive proportions per sample...\n")
  stats.df <- data.df %>%
    dplyr::group_by(condition, sample_id, celltype) %>%
    dplyr::summarise(
      total_count = dplyr::n(),
      pos_count = sum(expression > 0),
      pct_pos = (pos_count / total_count) * 100,
      .groups = "drop"
    ) %>%
    dplyr::filter(total_count >= min_cells)
  
  if (nrow(stats.df) == 0) {
    stop("Error: No data remains after applying min_cells filter.")
  }

  # Pivot to wide format for correlation analysis
  wide.df <- stats.df %>%
    dplyr::select(condition, sample_id, celltype, pct_pos) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = pct_pos)
  
  cat("[3/4] Computing Pearson correlations...\n")
  results_list <- list()
  conditions <- unique(wide.df$condition)
  
  # Determine cell type pairs to analyze
  if (only.significant) {
    active_cells <- intersect(valid_cells, colnames(wide.df))
    if(length(active_cells) < 2) stop("Insufficient cell types available for correlation.")
    cell_pairs <- combn(active_cells, 2, simplify = FALSE)
  } else {
    if (is.null(A.cell) | is.null(B.cell)) stop("A.cell and B.cell must be specified.")
    cell_pairs <- list(c(A.cell, B.cell))
  }
  
  # Iterate through conditions and cell type pairs
  for (cond in conditions) {
    sub_data <- wide.df %>% dplyr::filter(condition == cond)
    for (pair in cell_pairs){
      c1 <- pair[1]; c2 <- pair[2]
      
      if (c1 %in% names(sub_data) && c2 %in% names(sub_data)) {
        test_data <- sub_data %>% dplyr::select(dplyr::all_of(c(c1, c2))) %>% na.omit()
        
        # Filter by min_samples to prevent overfitting bias
        if (nrow(test_data) >= min_samples) {
          sd1 <- sd(test_data[[c1]], na.rm = TRUE)
          sd2 <- sd(test_data[[c2]], na.rm = TRUE)
          
          if (!is.na(sd1) && !is.na(sd2) && sd1 > 0 && sd2 > 0) {
            test_res <- cor.test(test_data[[c1]], test_data[[c2]], method = "pearson")
            results_list[[length(results_list) + 1]] <- data.frame(
              Condition = cond, 
              Cell_1 = c1, 
              Cell_2 = c2,
              R_value = round(as.numeric(test_res$estimate), 3),
              P_value = as.numeric(test_res$p.value),
              N_samples = nrow(test_data)
            )
          }
        }
      }
    }
  }
  
  if (length(results_list) == 0) stop("No correlations met the min_samples criteria.")
  all_results <- do.call(rbind, results_list) %>% dplyr::arrange(P_value)
  
  # [4/4] Visualization Logic
  cat(sprintf("[4/4] Generating regression plot: %s vs %s\n", B.cell, A.cell))
  
  plot.df <- stats.df %>%
    dplyr::filter(celltype %in% c(A.cell, B.cell)) %>%
    dplyr::select(condition, sample_id, celltype, pct_pos) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = pct_pos) %>%
    na.omit()
  
  # Check if enough data exists for plotting
  if (nrow(plot.df) < min_samples) {
    p <- ggplot() + 
      labs(title = paste("Insufficient Data (N <", min_samples, ")")) + 
      theme_void()
  } else {
    if (!is.null(color.use) && !is.null(names(color.use))) {
      plot.df$condition <- factor(plot.df$condition, levels = names(color.use))
    }
    
    p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = !!rlang::sym(B.cell), y = !!rlang::sym(A.cell))) +
      ggplot2::geom_point(ggplot2::aes(color = condition), size = 3, alpha = 0.8) +
      ggplot2::geom_smooth(method = "lm", color = "black", fill = "lightgray", linewidth = 1, formula = y ~ x) +
      ggplot2::facet_wrap(~condition, scales = "fixed") +
      ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
                       size = 4, fontface = "bold") +
      ggplot2::labs(
        title = paste0("Cohort-Level Co-expression: ", target_gene),
        subtitle = paste0(B.cell, " (X) vs ", A.cell, " (Y)"),
        x = paste0(target_gene, "+ in ", B.cell, " (%)"),
        y = paste0(target_gene, "+ in ", A.cell, " (%)"),
        color = "Group"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom"
      )
    
    if (!is.null(color.use)) p <- p + ggplot2::scale_color_manual(values = color.use)
  }
  
  cat("Analysis complete.\n")
  return(list(cor_table = all_results, plot = p))
}
