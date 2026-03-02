#' @title Cell-Type Co-expression Correlation Analysis (sc.celltype.cor)
#'
#' @description 
#' This function performs a systemic co-expression network analysis of a target gene 
#' across diverse cell types within the tumor microenvironment (TME) or clinical cohorts. 
#' It calculates the proportion of cells expressing the target gene per biological replicate 
#' and computes Pearson correlations between cell type pairs to identify coordinated 
#' expression dynamics.
#' 
#' @details 
#' **Dual Modes of Operation:**
#' 1. **Targeted Mode (`only.significant = FALSE`):** Validates a specific hypothesis by 
#'    visualizing the linear regression between two predefined cell types (`A.cell` and `B.cell`).
#' 2. **Auto-Screening Mode (`only.significant = TRUE`):** Acts as a powerful discovery tool. 
#'    It scans all possible pairwise combinations of cell types, calculates their correlation 
#'    coefficients and p-values across conditions, and extracts the significant networks 
#'    (`p_value <= p_cutoff`). To prevent graphical memory overload, it saves each significant 
#'    plot as an individual element within a named list.
#' 
#' **Biological Significance:**
#' Identifying strong positive correlations of inhibitory receptors (e.g., LAIR1, VSIR) 
#' between specific lineages (e.g., Progenitors and Regulatory T cells) exclusively within 
#' a disease state provides robust evidence for synchronized, systemic immunosuppressive 
#' cross-talk.
#'
#' @param seurat_obj A Seurat object. Ensure the default assay is set to "RNA".
#' @param target_gene Character. The gene of interest to analyze (e.g., "VSIR", "LAIR1").
#' @param A.cell Character. The first cell type (Y-axis). Ignored if only.significant = TRUE. Default: NULL.
#' @param B.cell Character. The second cell type (X-axis). Ignored if only.significant = TRUE. Default: NULL.
#' @param celltype_col Character. Metadata column name for cell types. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param sample_col Character. Metadata column name for biological replicates. Default: "orig.ident".
#' @param min_cells Integer. Minimum total cells required in a sample to be included in the analysis. Default: 0.
#' @param only.significant Logical. If TRUE, scans all pairs and saves plots meeting the p_cutoff into a list. Default: FALSE.
#' @param p_cutoff Numeric. The p-value threshold for filtering significant pairs. Default: 0.06.
#' @param color.use Vector. Custom color palette for experimental groups. Default: NULL.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cor_table}: A comprehensive data frame of pairwise Pearson correlations and significance flags.
#'   \item \code{plot}: If `only.significant = FALSE`, a single ggplot2 object. 
#'         If `only.significant = TRUE`, a named list of ggplot2 objects accessible via `$plot$CellA_CellB`.
#' }
#' @author Hyundong Yoon
#' @export
#'
#' @examples
#' # Example 1: Auto-Screening for significant VSIR networks
#' # auto_results <- sc.celltype.cor(
#' #   seurat_obj = Integ_object, 
#' #   target_gene = "VSIR",
#' #   only.significant = TRUE,
#' #   color.use = c("Healthy.donor"="#003a6f", "MDS"="#d67043", "sAML"="darkred")
#' # )
#' # 
#' # # Access a specific significant pair's plot:
#' # auto_results$plot$Progenitor_Regulatory.Tcell
#' #
#' # Example 2: Targeted validation of a specific pair
#' # specific_results <- sc.celltype.cor(
#' #   seurat_obj = Integ_object, 
#' #   target_gene = "LAIR1",
#' #   A.cell = "Progenitor",
#' #   B.cell = "Classical.Monocyte"
#' # )

sc.celltype.cor <- function(seurat_obj,
                            target_gene,
                            A.cell = NULL,
                            B.cell = NULL,
                            celltype_col = "detailed.celltypes",
                            group_col = "condition",
                            sample_col = "orig.ident",
                            min_cells = 0,
                            only.significant = FALSE,
                            p_cutoff = 0.06,
                            color.use = NULL) {
  
  # Check required dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Please install 'ggpubr'.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Please install 'rlang'.")
  
  cat(sprintf("[1] Extracting expression data for %s...\n", target_gene))
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  vars_to_fetch <- c(sample_col, group_col, celltype_col, target_gene)
  data.df <- Seurat::FetchData(seurat_obj, vars = vars_to_fetch)
  colnames(data.df) <- c("sample_id", "condition", "celltype", "expression")
  
  valid_cells <- unique(as.character(data.df$celltype))
  
  # Validation for targeted mode
  if (!only.significant) {
    if (is.null(A.cell) || is.null(B.cell)) stop("A.cell and B.cell must be provided if only.significant = FALSE.")
    if (!(A.cell %in% valid_cells) || !(B.cell %in% valid_cells)) stop("Target cell types not found in the specified celltype column.")
  }
  
  cat("[2] Calculating target gene positive proportions per sample...\n")
  stats.df <- data.df %>%
    dplyr::filter(celltype %in% valid_cells) %>%
    dplyr::group_by(condition, sample_id, celltype) %>%
    dplyr::summarise(
      total_count = dplyr::n(),
      pos_count = sum(expression > 0),
      pct_pos = pos_count / total_count * 100,
      .groups = "drop"
    ) %>%
    dplyr::filter(total_count >= min_cells)
  
  wide.df <- stats.df %>%
    dplyr::select(condition, sample_id, celltype, pct_pos) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = pct_pos)
  
  cat("[3] Computing all pairwise Pearson correlations...\n")
  results_list <- list()
  conditions <- unique(wide.df$condition)
  cell_pairs <- combn(valid_cells, 2, simplify = FALSE) 
  
  for (cond in conditions) {
    sub_data <- wide.df %>% dplyr::filter(condition == cond)
    for (pair in cell_pairs){
      cell_1 <- pair[1]
      cell_2 <- pair[2]
      
      if (cell_1 %in% names(sub_data) && cell_2 %in% names(sub_data)) {
        test_data <- sub_data %>%
          dplyr::select(dplyr::all_of(c(cell_1, cell_2))) %>%
          na.omit()
        
        if (nrow(test_data) >= 3) {
          if(sd(test_data[[cell_1]]) > 0 && sd(test_data[[cell_2]]) > 0) {
            test_res <- cor.test(test_data[[cell_1]], test_data[[cell_2]], method = "pearson")
            results_list[[length(results_list) + 1]] <- data.frame(
              Condition = cond,
              Cell_1 = cell_1,
              Cell_2 = cell_2,
              R_value = as.numeric(test_res$estimate),
              P_value = as.numeric(test_res$p.value),
              N_samples = nrow(test_data)
            )
          }
        }
      }
    }
  }
  
  all_results <- do.call(rbind, results_list)
  clean_view <- all_results %>%
    dplyr::mutate(
      Correlation_Type = dplyr::case_when(
        R_value > 0.3  ~ "Positive",
        R_value < -0.3 ~ "Negative",
        TRUE ~ "Weak/None"
      ),
      Significance = dplyr::case_when(
        P_value < 0.05 ~ "Significant (<0.05)",
        P_value <= p_cutoff ~ sprintf("Marginal (<=%s)", p_cutoff),
        TRUE ~ "Insignificant"
      )
    ) %>%
    dplyr::arrange(Condition, P_value, dplyr::desc(abs(R_value)))
  
  # --- Visualization Helper Function ---
  create_scatter <- function(y_cell, x_cell) {
    plot.df <- stats.df %>%
      dplyr::filter(celltype %in% c(y_cell, x_cell)) %>%
      dplyr::select(condition, sample_id, celltype, pct_pos) %>%
      tidyr::pivot_wider(names_from = celltype, values_from = pct_pos) %>%
      na.omit()
    
    p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = !!rlang::sym(x_cell), y = !!rlang::sym(y_cell))) +
      ggplot2::geom_point(ggplot2::aes(color = condition), size = 2.5, alpha = 0.8) +
      ggplot2::geom_smooth(method = "lm", color = "black", fill = "lightgray", linewidth = 0.8, formula = y ~ x) +
      ggplot2::facet_wrap(~condition, scales = "fixed") +
      ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5) +
      ggplot2::labs(
        title = paste0(x_cell, " vs ", y_cell),
        x = paste0(target_gene, "+ in ", x_cell, " (%)"),
        y = paste0(target_gene, "+ in ", y_cell, " (%)"),
        color = "Group"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        strip.text = ggplot2::element_text(size = 11, face = "bold"),
        legend.position = "bottom"
      )
    
    # Apply user-defined color palette if provided
    if (!is.null(color.use)) {
      p <- p + ggplot2::scale_color_manual(values = color.use)
    }
    return(p)
  }

  # --- Plot Generation Logic ---
  if (only.significant) {
    cat(sprintf("[4] only.significant = TRUE. Screening for pairs with P-value <= %s in any condition...\n", p_cutoff))
    
    sig_pairs_df <- all_results %>%
      dplyr::filter(P_value <= p_cutoff) %>%
      dplyr::select(Cell_1, Cell_2) %>%
      dplyr::distinct()
    
    if (nrow(sig_pairs_df) == 0) {
      cat("No significant pairs found under the defined p-value cutoff.\n")
      final_plot <- NULL
    } else {
      cat(sprintf("Found %d significant pair(s). Storing plots in a named list...\n", nrow(sig_pairs_df)))
      
      plot_list <- list()
      plot_names <- character(nrow(sig_pairs_df))
      
      for (i in seq_len(nrow(sig_pairs_df))) {
        c1 <- sig_pairs_df$Cell_1[i]
        c2 <- sig_pairs_df$Cell_2[i]
        
        # Create individual plot and store in list
        plot_list[[i]] <- create_scatter(c1, c2)
        # Construct the requested naming convention: Cell1_Cell2
        plot_names[i] <- paste0(c1, "_", c2) 
      }
      
      # Assign names to the list elements
      names(plot_list) <- plot_names
      final_plot <- plot_list
    }
    
    # Filter the returned cor_table to highlight only these significant pairs
    clean_view <- clean_view %>%
      dplyr::inner_join(sig_pairs_df, by = c("Cell_1", "Cell_2"))
    
  } else {
    cat(sprintf("[4] Generating linear regression scatter plot for %s vs %s...\n", B.cell, A.cell))
    final_plot <- create_scatter(A.cell, B.cell)
  }
  
  cat("[Done] Successfully computed network and generated plots.\n")
  
  return(list(
    cor_table = clean_view,
    plot = final_plot
  ))
}
