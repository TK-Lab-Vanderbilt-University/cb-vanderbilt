#' @title Gene Expression Proportion Barplot (Percentage)
#'
#' @description
#' This function calculates and visualizes the percentage of cells expressing 
#' a target gene (counts > 0) across different groups. It provides a clear 
#' summary of gene prevalence within specific cell types or the entire population.
#' 
#' Key Features:
#' 1. Percentage Calculation: Computes (Number of cells with counts > 0 / Total cells) * 100.
#' 2. Multi-gene Support: Automatically creates faceted plots for multiple genes.
#' 3. Automatic Annotation: Displays the exact percentage value on top of each bar.
#' 4. Group Comparison: Visualizes the spread of expression across defined conditions.
#'
#' @param seurat_obj A Seurat object.
#' @param target_gene String or vector. Gene name(s) to analyze.
#' @param group_var String. Metadata column for groups (default: "condition").
#' @param celltype_var String. Metadata column for cell types (default: "detailed.celltypes").
#' @param celltypes String or vector. "all" for global analysis, or a vector of specific cell types.
#' @param color_pal Custom color palette for groups.
#' @param show_legend Logical. Toggle legend display (default: TRUE).
#' 
#' @return A ggplot2 object.
#' @author Hyundong Yoon
#' @export
sc.expr.prop.bar <- function(seurat_obj,
                             target_gene,
                             group_var = "condition",
                             celltype_var = "detailed.celltypes",
                             celltypes = "all",
                             color_pal = NULL,
                             show_legend = TRUE) {
  
  require(dplyr)
  require(ggplot2)
  require(Seurat)
  
  # 1. Parameter handling & Validation
  missing_genes <- target_gene[!target_gene %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    stop(paste0("Genes not found in object: ", paste(missing_genes, collapse = ", ")))
  }
  
  # 2. Filter Cells
  meta_data <- seurat_obj@meta.data
  if (length(celltypes) == 1 && celltypes == "all") {
    valid_cells <- rownames(meta_data)
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% celltypes, ])
    if (length(valid_cells) == 0) stop("No matching cell types found.")
  }
  
  # 3. Extract and Process Data
  raw_counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[target_gene, valid_cells, drop = FALSE]
  
  prop_list <- lapply(target_gene, function(g) {
    data.frame(
      group = meta_data[valid_cells, group_var],
      gene = g,
      expressed = as.numeric(raw_counts[g, ]) > 0
    ) %>%
      group_by(group, gene) %>%
      summarise(
        pct = sum(expressed) / n() * 100,
        .groups = "drop"
      )
  })
  
  prop_df <- do.call(rbind, prop_list) %>% filter(!is.na(group))
  
  # 4. Visualization
  p <- ggplot(prop_df, aes(x = gene, y = pct, fill = group)) +
    geom_col(color = "black", size = 0.5, width = 0.8, 
             position = position_dodge(width = 0.85)) +
    
    # Numerical annotation on top
    geom_text(aes(label = sprintf("%.1f%%", pct)), 
              position = position_dodge(width = 0.85),
              vjust = -0.5, 
              fontface = "bold", 
              size = 3.5) +
    
    theme_classic() +
    labs(
      title = NULL,     # 제목 제거
      subtitle = NULL,  # 부제목 제거
      y = "Percent of Cells (%)",
      x = "Genes"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    theme(
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      legend.title = element_blank()
    )
  
  # Apply color palette
  if (!is.null(color_pal)) {
    p <- p + scale_fill_manual(values = color_pal)
  }
  
  # Legend control
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}