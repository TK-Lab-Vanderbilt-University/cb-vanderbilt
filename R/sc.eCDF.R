#' @title Quantitative Visualization of Expression Magnitude Shifts via eCDF
#'
#' @description 
#' This function generates an Empirical Cumulative Distribution Function (eCDF) plot 
#' to rigorously evaluate gene expression magnitude shifts across biological conditions. 
#' By utilizing $\log_2(\text{CPM} + 1)$ transformation directly from raw counts, it provides 
#' a standardized metric that ensures linear interpretability of fold changes, 
#' bypassing the stochastic noise and normalization artifacts often present in 
#' standard density-based visualizations like Violin plots.
#'
#' @details 
#' **Statistical Rationale:**
#' In clinical cohorts (e.g., Healthy vs. MDS vs. sAML), traditional Violin plots often 
#' obscure subtle expression shifts due to kernel density estimation and the high 
#' frequency of dropouts. This function isolates expressing cells (count > 0) to 
#' characterize the true distribution of the transcriptional magnitude. 
#' The eCDF is a non-parametric approach that captures the entire spectrum of 
#' population-level shifts, making it significantly more sensitive for detecting 
#' progressive up/down-regulation in longitudinal or disease-progression studies.
#'
#' @param seurat_obj A Seurat object containing a 'counts' layer in the RNA assay.
#' @param target_gene Character. The gene of interest for distribution analysis.
#' @param target_celltype Character. The specific cell population to analyze. 
#' Set to "All" to perform a global analysis across the entire dataset. Default: "All".
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param group_order Character vector. The desired ordering of experimental groups 
#' to visualize progression (e.g., c("Healthy.donor", "MDS", "sAML")).
#' @param custom_colors Character vector. Optional hex codes for group colors. 
#' If NULL, an optimized clinical progression palette is applied. Default: NULL.
#'
#' @return A publication-quality ggplot2 object representing the eCDF distributions.
#' @export
#'
#' @examples
#' # [Example 1: Analyzing global VSIR shift across disease states]
#' # p_global <- sc.eCDF(
#' #   seurat_obj = Integ_object, 
#' #   target_gene = "VSIR", 
#' #   target_celltype = "All",
#' #   group_order = c("Healthy.donor", "MDS", "sAML")
#' # )
#' # 
#' # [Example 2: Analyzing cell-type specific LAIR1 shift]
#' # p_celltype <- sc.eCDF(
#' #   seurat_obj = Integ_object, 
#' #   target_gene = "LAIR1", 
#' #   target_celltype = "Progenitor"
#' # )

sc.eCDF <- function(seurat_obj,
                    target_gene,
                    target_celltype = "All",
                    celltype_col = "detailed.celltypes",
                    group_col = "condition",
                    group_order = c("Healthy.donor", "MDS", "sAML"),
                    custom_colors = NULL) {
  
  # Dependency validation
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")
  
  message(sprintf("Computing log2(CPM + 1) for %s...", target_gene))
  
  # Ensure target assay is available
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  raw_counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
  
  if (!(target_gene %in% rownames(raw_counts))) {
    stop(sprintf("Gene '%s' not found in the counts layer of the Seurat object.", target_gene))
  }
  
  # 1. On-the-fly Gold Standard Normalization: log2(CPM + 1)
  # This avoids distorting the distribution with natural log or arbitrary scaling factors.
  gene_vector <- raw_counts[target_gene, ]
  total_counts <- seurat_obj$nCount_RNA
  log2_cpm_val <- log2(((gene_vector / total_counts) * 1e6) + 1)
  
  # 2. Data Integration and Pre-processing
  plot_data <- data.frame(
    expression = log2_cpm_val,
    group = seurat_obj@meta.data[[group_col]],
    stringsAsFactors = FALSE
  )
  
  # Apply cell-type filter or global aggregate
  if (target_celltype != "All") {
    plot_data$celltype <- seurat_obj@meta.data[[celltype_col]]
    
    sub_data <- plot_data %>%
      dplyr::filter(celltype == target_celltype, 
                    group %in% group_order,
                    expression > 0) %>% 
      dplyr::mutate(group = factor(group, levels = group_order))
    
    subtitle_txt <- sprintf("Magnitude shift in expressing %s cells (counts > 0)", target_celltype)
  } else {
    sub_data <- plot_data %>%
      dplyr::filter(group %in% group_order,
                    expression > 0) %>% 
      dplyr::mutate(group = factor(group, levels = group_order))
    
    subtitle_txt <- "Global magnitude shift across all expressing cells (counts > 0)"
  }
  
  if (nrow(sub_data) == 0) {
    stop("Insufficient expressing cells (count > 0) identified for the given criteria.")
  }
  
  # 3. Esthetic Configuration (Clinical Progression Palette)
  if (is.null(custom_colors)) {
    # Blue (Healthy) -> Orange (Disease) -> Red (Severe/Progression)
    progression_palette <- c("#003a6f", "#d67043", "#a50f15", "#4a1486", "#006d2c")
    my_colors <- progression_palette[1:length(group_order)]
    names(my_colors) <- group_order
  } else {
    if (length(custom_colors) < length(group_order)) {
      stop("Custom color vector length must match the number of unique groups in group_order.")
    }
    my_colors <- custom_colors
    names(my_colors) <- group_order
  }
  
  # 4. Visualization Construction
  message("Generating eCDF visualization...")
  
  p <- ggplot2::ggplot(sub_data, ggplot2::aes(x = expression, color = group)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = my_colors) +
    ggplot2::labs(
      title = paste0("Empirical Cumulative Distribution: ", target_gene),
      subtitle = subtitle_txt,
      y = "Cumulative Fraction of Cells", 
      x = expression(bold(log[2] * (CPM + 1))),
      color = "Condition"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray30"),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.title.y = ggplot2::element_text(face = "bold", size = 12),
      axis.text = ggplot2::element_text(size = 11),
      legend.position = c(0.75, 0.25),
      # Explicit namespacing with scales::alpha to prevent function masking conflicts
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.7), 
                                                color = "black", linewidth = 0.5),
      legend.title = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}