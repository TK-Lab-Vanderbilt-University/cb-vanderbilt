#' @title Custom Polar Volcano Plot with Radial Labels and log2FC Scale
#' @description Generates a polar volcano plot with manual grid lines to indicate log2FC magnitude.
#' 
#' @param diffData Dataframe of DE results.
#' @param target_genes Character vector of genes to highlight.
#' @param color.use Named character vector of colors for cell types.
#' @param grid_breaks Numeric vector. Values to draw scale circles (e.g., c(-2, 2, 4)).
#' @param hole_size Numeric. Controls the center hole radius. Default is 5.5.
#' @param dot_size Numeric. Size of the target gene points. Default is 3.
sc.Polar.Volcano <- function(diffData, 
                             target_genes, 
                             color.use, 
                             grid_breaks = c(-2, 2, 4, 6),
                             log2FC_cut = 0.25, 
                             pval_cut = 0.05,
                             hole_size = 5.5,
                             dot_size = 3) {
  
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  
  # 1. Data Processing
  plot_data <- diffData %>%
    rename(gene = genes, cluster = celltype) %>%
    mutate(cluster = factor(cluster, levels = names(color.use)))
  
  # 2. Background Bar Range Calculation
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
  label_data <- back_data %>%
    mutate(
      id = as.numeric(cluster),
      n_tot = n(),
      angle_base = 90 - 360 * (id - 0.5) / n_tot,
      angle = ifelse(angle_base < -90, angle_base + 180, angle_base),
      hjust = ifelse(angle_base < -90, 0, 1),
      y_pos = min_val - 0.5 
    )
  
  # 4. Filter Target Gene Data
  target_data <- plot_data %>% filter(gene %in% target_genes)
  
  # 5. Build Plot
  p <- ggplot() +
    # Background Bars
    geom_col(data = back_data, aes(x = cluster, y = max_val), fill = "grey93", width = 0.9) +
    geom_col(data = back_data, aes(x = cluster, y = min_val), fill = "grey93", width = 0.9) +
    
    # [NEW] Manual Grid Lines (Concentric Circles)
    geom_hline(yintercept = grid_breaks, color = "grey80", linetype = "dashed", linewidth = 0.3) +
    
    # [NEW] Scale Labels (avg_log2FC numbers)
    # Placed at the position of the first cluster for consistency
    annotate("text", x = 0.5, y = grid_breaks, label = grid_breaks, 
             size = 2.5, color = "grey40", fontface = "italic") +
    
    # Baseline (log2FC = 0)
    geom_hline(yintercept = 0, color = "white", linewidth = 1.2) +
    
    # Center Cluster Tiles
    geom_tile(data = back_data, aes(x = cluster, y = 0, fill = cluster),
              height = 0.6, color = "black", alpha = 0.6, show.legend = FALSE) +
    
    # Radial Cell Type Labels
    geom_text(data = label_data, 
              aes(x = cluster, y = y_pos, label = cluster, angle = angle, hjust = hjust),
              size = 2.8, color = "black") +
    
    # Target Gene Points
    geom_point(data = target_data, aes(x = cluster, y = avg_log2FC, color = cluster),
               size = dot_size) +
    
    # Target Gene Text Labels
    ggrepel::geom_text_repel(data = target_data,
                             aes(x = cluster, y = avg_log2FC, label = gene),
                             size = 4, fontface = "bold.italic",
                             box.padding = 0.8, min.segment.length = 0) +
    
    # Color Scales
    scale_fill_manual(values = color.use) +
    scale_color_manual(values = color.use) +
    
    # Coordinates & Theme
    scale_y_continuous(limits = c(min(back_data$min_val) - hole_size, max(back_data$max_val) + 1),
                       expand = c(0, 0)) +
    coord_polar(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "none"
    )
  
  return(p)
}
