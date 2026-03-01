#' @title Cell-Type Co-expression Correlation Analysis (sc.celltype.cor)
#'
#' @description 
#' Analyzes the systemic co-expression network of a target gene across different cell types 
#' within the tumor microenvironment (TME) or clinical cohorts. It calculates the proportion 
#' of cells expressing the target gene per biological replicate and computes Pearson 
#' correlations between all possible cell type pairs. It also generates a publication-ready 
#' linear regression scatter plot for a specific pair of interest.
#'
#' @details 
#' **Biological Significance:**
#' This function is crucial for identifying coordinated immune evasion mechanisms or 
#' systemic cross-talk. For instance, a strong positive correlation of an inhibitory 
#' receptor (e.g., VSIR, LAIR1) between Progenitors and Tregs specifically in a disease 
#' state (like MDS) provides strong evidence for a synchronized immunosuppressive network.
#'
#' @param seurat_obj A Seurat object. Ensure the default assay is set to "RNA".
#' @param target_gene Character. The gene of interest (e.g., "VSIR", "LAIR1").
#' @param A.cell Character. The first cell type of interest for the scatter plot.
#' @param B.cell Character. The second cell type of interest for the scatter plot.
#' @param celltype_col Character. Metadata column for cell types. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column for experimental groups. Default: "condition".
#' @param sample_col Character. Metadata column for biological replicates. Default: "orig.ident".
#' @param min_cells Integer. Minimum number of total cells required in a sample to be included. Default: 0.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{cor_table}: A data frame of all significant pairwise Pearson correlations.
#'   \item \code{plot}: A ggplot2 object displaying the linear regression of the target cell pair.
#' }
#' @export
#'
#' @examples
#' # 1. Execute the function
#' # VSIR_correlation <- sc.celltype.cor(
#' #   seurat_obj = Integ_object, 
#' #   target_gene = "VSIR", 
#' #   A.cell = "Progenitor", 
#' #   B.cell = "Regulatory.Tcell"
#' # )
#' # 
#' # # 2. Display the visualization graph (can be saved directly)
#' # VSIR_correlation$plot
#' # 
#' # # 3. Check the significant network table between all cells
#' # View(VSIR_correlation$cor_table)
#' # 
#' # # 4. Export the table to a file
#' # # write.table(VSIR_correlation$cor_table, file = '~/clean_view.txt', quote = FALSE, row.names = FALSE, sep = '\t')

sc.celltype.cor <- function(seurat_obj,
                                   target_gene,
                                   A.cell,
                                   B.cell,
                                   celltype_col = "detailed.celltypes",
                                   group_col = "condition",
                                   sample_col = "orig.ident",
                                   min_cells = 0) {
  
  # Check required dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Please install 'ggpubr' for stat_cor().")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Please install 'rlang'.")
  
  cat(sprintf("[1] Extracting expression data for %s...\n", target_gene))
  
  # Set default assay to RNA safely
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # Safely fetch data and standardize column names
  vars_to_fetch <- c(sample_col, group_col, celltype_col, target_gene)
  data.df <- Seurat::FetchData(seurat_obj, vars = vars_to_fetch)
  colnames(data.df) <- c("sample_id", "condition", "celltype", "expression")
  
  valid_cells <- unique(as.character(data.df$celltype))
  
  if (!(A.cell %in% valid_cells) || !(B.cell %in% valid_cells)) {
    stop("Target cell types (A.cell or B.cell) not found in the specified celltype column.")
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
        P_value < 0.05 ~ "Significant",
        TRUE ~ "Insignificant"
      )
    ) %>%
    dplyr::arrange(Condition, Significance, dplyr::desc(abs(R_value))) %>%
    dplyr::filter(Correlation_Type != "Weak/None")
  
  cat(sprintf("[4] Generating linear regression scatter plot for %s vs %s...\n", B.cell, A.cell))
  
  plot.df <- stats.df %>%
    dplyr::filter(celltype %in% c(A.cell, B.cell)) %>%
    dplyr::select(condition, sample_id, celltype, pct_pos) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = pct_pos) %>%
    na.omit()
  
  # Ensure standard naming if condition matches the user's specific colors
  color_mapping <- c("MDS" = "#d67043", "Healthy.donor" = "#003a6f")
  
  p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = !!rlang::sym(B.cell), y = !!rlang::sym(A.cell))) +
    ggplot2::geom_point(ggplot2::aes(color = condition), size = 3, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", color = "black", fill = "lightgray", linewidth = 0.8) +
    ggplot2::facet_wrap(~condition, scales = "fixed") +
    ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
    ggplot2::labs(
      title = paste0("Correlation of ", target_gene, "+ Proportion"),
      subtitle = paste0(B.cell, " vs ", A.cell),
      x = paste0(target_gene, "+ in ", B.cell, " (%)"),
      y = paste0(target_gene, "+ in ", A.cell, " (%)"),
      color = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
  
  # Apply specific colors if MDS and Healthy.donor are present
  if (all(unique(plot.df$condition) %in% names(color_mapping))) {
    p <- p + ggplot2::scale_color_manual(values = color_mapping)
  }
  
  cat("[Done] Successfully computed network and generated plot.\n")
  
  return(list(
    cor_table = clean_view,
    plot = p
  ))
}