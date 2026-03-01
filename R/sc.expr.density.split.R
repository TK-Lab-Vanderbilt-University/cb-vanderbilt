#' @title Split Density Plot for Single-Cell Data (Nebulosa Based)
#'
#' @description
#' This function generates Kernel Density Estimation (KDE) plots for gene expression, 
#' split by a metadata variable (e.g., condition). It overcomes the "overplotting" 
#' issue of standard FeaturePlots and provides a clearer view of expression hot-spots 
#' across different experimental groups.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene String. The gene to visualize.
#' @param split_by String. Metadata column to split the plots by (e.g., "condition").
#' @param nrow Number of rows for the layout. Default is 1 (Horizontal layout).
#' @param ncol Number of columns for the layout. Default is NULL.
#' @param viridis_opt String. Viridis color scale: "magma", "plasma", "inferno", "viridis", "cividis".
#' @param show_legend Logical. Whether to collect and show a single master legend.
#' 
#' @return A patchwork object combining density plots for each group.
#' @author Hyundong Yoon
#' @export
sc.expr.density.split <- function(seurat_obj,
                                  target_gene,
                                  split_by = "condition",
                                  nrow = 1,
                                  ncol = NULL,
                                  viridis_opt = "magma",
                                  show_legend = TRUE) {
  
  require(Nebulosa)
  require(patchwork)
  require(ggplot2)
  require(Seurat)
  
  # 1. Validation
  if (!target_gene %in% rownames(seurat_obj)) stop("Gene not found in the object.")
  if (!split_by %in% colnames(seurat_obj@meta.data)) stop("split_by column not found in metadata.")
  
  # 2. Get Group Levels (following factor levels if available)
  if (is.factor(seurat_obj@meta.data[[split_by]])) {
    group_levels <- levels(seurat_obj@meta.data[[split_by]])
  } else {
    group_levels <- sort(unique(seurat_obj@meta.data[[split_by]]))
  }
  
  # 3. Generate Density Plot for each group
  plot_list <- lapply(group_levels, function(group) {
    
    # Subset for the specific group
    sub_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[split_by]] == group, ]))
    
    # Create Nebulosa density plot
    p <- plot_density(sub_obj, target_gene, joint = FALSE) +
      labs(title = group) + # Set group name as panel title
      theme_classic() +
      scale_color_viridis_c(option = viridis_opt) + # High contrast colors
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )
    
    return(p)
  })
  
  # 4. Combine plots using patchwork
  combined <- wrap_plots(plot_list, nrow = nrow, ncol = ncol)
  
  # 5. Global Layout adjustments
  if (show_legend) {
    combined <- combined + plot_layout(guides = "collect") & 
      theme(legend.position = "bottom")
  } else {
    combined <- combined & theme(legend.position = "none")
  }
  
  # Remove main title as requested
  combined <- combined + plot_annotation(title = NULL)
  
  return(combined)
}