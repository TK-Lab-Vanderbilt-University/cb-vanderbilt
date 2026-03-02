#' @title Visualize AUCell Scores with Statistical Tests (sc.PlotAUCell)
#'
#' @description 
#' Generates publication-ready violin plots overlaid with boxplots to visualize 
#' AUCell scores across different groups/clusters. It automatically computes and 
#' adds pairwise statistical comparisons (e.g., p-values, significance stars) 
#' using `ggpubr`.
#'
#' @details 
#' **Statistical Considerations:**
#' Single-cell RNA-seq data often violates the assumptions of normality. Therefore, 
#' the default statistical test is set to the non-parametric `wilcox.test` 
#' (Mann-Whitney U test). You can change it to `t.test` if preferred.
#'
#' @param seurat_obj A Seurat object containing AUCell scores in its metadata.
#' @param target_pathway Character. The name of the pathway/score to visualize (e.g., "HALLMARK_P53_PATHWAY").
#' @param group_col Character. The metadata column used for grouping cells (e.g., "condition", "seurat_clusters").
#' @param comparisons A list of character vectors specifying which groups to compare (e.g., `list(c("GroupA", "GroupB"))`).
#' @param method Character. The statistical test to perform. Options: "wilcox.test" (default) or "t.test".
#' @param palette Character. Color palette for `ggpubr` (default: "npg").
#'
#' @return A `ggplot` object.
#' @export
#'
#' @import ggplot2 ggpubr
#'
#' @examples
#' # [Visualization Example]
#' # my_pairs <- list(c("Healthy.donor", "MDS"), c("MDS", "sAML"), c("Healthy.donor", "sAML"))
#' # p <- sc.PlotAUCell(seurat_obj = Integ_object,
#' #                    target_pathway = "HALLMARK_P53_PATHWAY",
#' #                    group_col = "condition",
#' #                    comparisons = my_pairs)
#' # print(p)
sc.PlotAUCell <- function(seurat_obj, 
                          target_pathway, 
                          group_col, 
                          comparisons, 
                          method = "wilcox.test", 
                          palette = "npg") {
  
  # 1. Extract metadata
  plot_data <- seurat_obj@meta.data
  
  # 2. Check if the pathway exists in metadata
  if (!target_pathway %in% colnames(plot_data)) {
    stop(sprintf("Error: '%s' not found in Seurat metadata. Did you run sc.RunAUCell first?", target_pathway))
  }
  
  # 3. Check if the grouping column exists
  if (!group_col %in% colnames(plot_data)) {
    stop(sprintf("Error: Grouping column '%s' not found in Seurat metadata.", group_col))
  }
  
  # 4. Generate ggviolin plot with stat_compare_means
  p <- ggpubr::ggviolin(
    data = plot_data, 
    x = group_col, 
    y = target_pathway, 
    fill = group_col,              
    palette = palette,               
    add = "boxplot",               
    add.params = list(fill = "white", width = 0.1)
  ) + 
    ggpubr::stat_compare_means(
      comparisons = comparisons, 
      method = method, 
      label = "p.signif" # Use "p.format" to show exact p-values
    ) +
    ggplot2::labs(
      title = gsub("_", " ", target_pathway), # Replaces underscores with spaces for cleaner titles
      x = "Group",
      y = "AUCell Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) # Prevents overlapping x-axis labels
    )
  
  return(p)
}