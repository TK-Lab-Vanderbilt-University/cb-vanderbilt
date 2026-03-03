#' @title Custom Polar Volcano Plot with Outside Labels
#' @description Generates a polar volcano plot where cell type names are placed on the outer rim.
#'              Target genes are color-coded: Up-regulated (darkred) and Down-regulated (darkblue).
#' 
#' @param diffData Dataframe of DE results (must contain 'genes', 'celltype', 'avg_log2FC', 'p_val').
#' @param target_genes Character vector of genes to highlight.
#' @param color.use Named character vector of colors for the center cluster tiles.
#' @param grid_breaks Numeric vector for log2FC guide circles.
#' @param hole_size Numeric. Controls the inner hole radius. Default is 2.0.
#' @param dot_size Numeric. Size of the target gene points. Default is 3.
sc.Polar.Volcano <- function(diffData, 
                             target_genes, 
                             color.use, 
                             grid_breaks = c(-2, 2, 4, 6),
                             log2FC_cut = 0.25, 
                             pval_cut = 0.05,
                             hole_size = 2.0,
                             dot_size = 3) {
  
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  
  # 1. Data Processing
  plot_data <- diffData %>%
    rename(gene = genes, cluster = celltype) %>%
    mutate(cluster = factor(cluster, levels = names(color.use)))
  
  # 2. Calculate Background Bar Ranges
  sig_data <- plot_data %>%
    filter(abs(avg_log2FC) >= log2FC_cut & p_val < pval_cut)
  
  back_data <- sig_data %>%
    group_by(cluster) %>%
    summarise(
      min_val = min(avg_log2FC) - 0.5,
      max_val = max(avg_log2FC) + 0.5,
      .groups = 'drop'
    )
  
  # 3. Outside Label Calculation
  # Positions labels beyond the maximum up-regulation bars
  label_data <- back_data %>%
    mutate(
      id = as.numeric(cluster),
      n_tot = n(),
      angle_base = 90 - 360 * (id - 0.5) / n_tot,
      angle = ifelse(angle_base < -90, angle_base + 180, angle_base),
      # For outside labels, hjust is flipped compared to inside labels
      hjust = ifelse(angle_base < -90, 1, 0),
      y_pos = max_val + 0.8 # Positioned outside the max bar
    )
  
  # 4. Filter and Color-Code Target Genes
  target_data <- plot_data %>% 
    filter(gene %in% target_genes) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))
  
  # 5. Build Plot
  p <- ggplot() +
    # Background Bars (Grey distribution)
    geom_col(data = back_data, aes(x = cluster, y = max_val), fill = "grey93", width = 0.9) +
    geom_col(data = back_data, aes(x = cluster, y = min_val), fill = "grey93", width = 0.9) +
    
    # log2FC Scale Guides (Dashed Circles)
    geom_hline(yintercept = grid_breaks, color = "grey85", linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = 0.5, y = grid_breaks, label = grid_breaks, 
             size = 2.2, color = "grey50", fontface = "italic") +
    
    # Baseline (log2FC = 0)
    geom_hline(yintercept = 0, color = "white", linewidth = 1.2) +
    
    # Center Cluster Tiles (The Color Ring)
    geom_tile(data = back_data, aes(x = cluster, y = 0, fill = cluster),
              height = 0.5, color = "black", alpha = 0.8, show.legend = FALSE) +
    
    # [UPDATED] Outside Cell Type Labels
    geom_text(data = label_data, 
              aes(x = cluster, y = y_pos, label = cluster, angle = angle, hjust = hjust),
              size = 2.8, color = "black", fontface = "bold") +
    
    # [UPDATED] Target Points (Colored by Up/Down)
    geom_point(data = target_data, 
               aes(x = cluster, y = avg_log2FC, color = direction),
               size = dot_size, show.legend = FALSE) +
    
    # [UPDATED] Target Gene Labels (Color matched to dots)
    ggrepel::geom_text_repel(data = target_data,
                             aes(x = cluster, y = avg_log2FC, label = gene, color = direction),
                             size = 4, fontface = "bold.italic",
                             box.padding = 0.8, min.segment.length = 0, 
                             show.legend = FALSE) +
    
    # Color Scales
    scale_fill_manual(values = color.use) +
    scale_color_manual(values = c("Up" = "darkred", "Down" = "darkblue")) +
    
    # Coordinates & Theme
    # Adjusted limits to leave room for outside labels
    scale_y_continuous(limits = c(min(back_data$min_val) - hole_size, max(back_data$max_val) + 3.5),
                       expand = c(0, 0)) +
    coord_polar(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(20, 20, 20, 20),
      legend.position = "none"
    )
  
  return(p)
}
