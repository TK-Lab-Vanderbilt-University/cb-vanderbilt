#' @title Custom Polar Volcano Plot with Radial Labels
#' @description Generates a high-quality polar volcano plot for single-cell DE results.
#'              Features outward-facing up-regulation, inward-facing down-regulation, 
#'              and radial cell type labeling in the center hole.
#' 
#' @param diffData Dataframe of DE results (must contain 'genes', 'celltype', 'avg_log2FC', 'p_val')
#' @param target_genes Character vector of genes to highlight (e.g., LAIR1, VSIR).
#' @param color.use Named character vector of colors for cell types.
#' @param log2FC_cut Numeric. log2FC threshold for background bars. Default is 0.25.
#' @param pval_cut Numeric. p-value threshold for background bars. Default is 0.05.
#' @param hole_size Numeric. Controls the radius of the center hole for labels. Default is 5.0.
#' @param dot_size Numeric. Size of the highlighted target gene points. Default is 5.
#' @return A ggplot object.
#' @author Hyundong Yoon
#' @export
sc.Polar.Volcano <- function(diffData, 
                             target_genes, 
                             color.use, 
                             log2FC_cut = 0.25, 
                             pval_cut = 0.05,
                             hole_size = 5.0,
                             dot_size = 5) {
  
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  
  # 1. Data Processing and Mapping
  # Rename columns and factorize cluster levels to match the provided color palette order
  plot_data <- diffData %>%
    rename(gene = genes, cluster = celltype) %>%
    mutate(cluster = factor(cluster, levels = names(color.use)))
  
  # 2. Calculate Background Bar Ranges
  # Define the range for the grey background bars based on significant genes
  sig_data <- plot_data %>%
    filter(abs(avg_log2FC) >= log2FC_cut & p_val < pval_cut)
  
  back_data <- sig_data %>%
    group_by(cluster) %>%
    summarise(
      min_val = min(avg_log2FC) - 0.5,
      max_val = max(avg_log2FC) + 0.5,
      .groups = 'drop'
    )
  
  # 3. Radial Label Calculation
  # Calculate rotation angles and horizontal justification for labels in polar coordinates
  label_data <- back_data %>%
    mutate(
      id = as.numeric(cluster),
      n_tot = n(),
      # Convert position to degrees for rotation
      angle_base = 90 - 360 * (id - 0.5) / n_tot,
      # Flip text on the left side of the plot to keep it readable (not upside down)
      angle = ifelse(angle_base < -90, angle_base + 180, angle_base),
      hjust = ifelse(angle_base < -90, 0, 1),
      # Position labels slightly inward from the minimum log2FC value
      y_pos = min_val - 0.3 
    )
  
  # 4. Filter Target Gene Data
  # Isolate the specific genes requested for dot highlighting
  target_data <- plot_data %>% filter(gene %in% target_genes)
  
  # 5. Build Plot Layers
  p <- ggplot() +
    # Background Bars: Representing the spread of DE genes for each cluster
    geom_col(data = back_data, aes(x = cluster, y = max_val), fill = "grey93", width = 0.9) +
    geom_col(data = back_data, aes(x = cluster, y = min_val), fill = "grey93", width = 0.9) +
    
    # Baseline: Indicates log2FC = 0
    geom_hline(yintercept = 0, color = "white", linewidth = 1.2) +
    
    # Center Tiles: Colored band representing each cell type
    geom_tile(data = back_data, aes(x = cluster, y = 0, fill = cluster),
              height = 0.6, color = "black", alpha = 0.6, show.legend = FALSE) +
    
    # Radial Labels: Cell type names positioned in the center hole
    geom_text(data = label_data, 
              aes(x = cluster, y = y_pos, label = cluster, angle = angle, hjust = hjust),
              size = 2.8, color = "black") +
    
    # Target Points: Large colored dots for LAIR1 and VSIR
    geom_point(data = target_data, aes(x = cluster, y = avg_log2FC, color = cluster),
               size = dot_size) +
    
    # Gene Labels: Non-overlapping text labels for target genes
    ggrepel::geom_text_repel(data = target_data,
                             aes(x = cluster, y = avg_log2FC, label = gene),
                             size = dot_size, fontface = "bold.italic",
                             box.padding = 1, min.segment.length = 0) +
    
    # Mapping custom colors to scales
    scale_fill_manual(values = color.use) +
    scale_color_manual(values = color.use) +
    
    # 6. Coordinate and Theme Settings
    # Adjust y-axis limits to expand the center hole for labels
    scale_y_continuous(limits = c(min(back_data$min_val) - hole_size, max(back_data$max_val) + 1),
                       expand = c(0, 0)) +
    coord_polar(clip = "off") +
    theme_void() + # Clean background for publication
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none" # Labels are directly on plot, so legend is redundant
    )
  
  return(p)
}