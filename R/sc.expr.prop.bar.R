#' @title Gene Expression Proportion Barplot (Percentage)
#' @description
#' Percentage of cells expressing target gene(s) (counts > 0) by group.
#' Error.bar = TRUE  : mean +/- SE (or SD) across samples (sample_var) + sample points
#' statistic = TRUE  : pairwise t-test between groups per gene (on per-sample percentages)
#'
#' @param Error.bar Logical. Add error bars (and sample points) based on sample-level percentages.
#' @param statistic Logical. Pairwise t-test between groups for each gene, with brackets.
#' @param sample_var Metadata column defining biological replicates (default "orig.ident").
#' @param error_type "se" or "sd".
#' @param p_label "star" (*, **, ***, ns) or "p" (numeric).
#' @param p_adjust Method for p.adjust across pairs within a gene (only matters when >2 groups).
#' @param min_cells Minimum cells for a sample to be used in sample-level mode.
#' @param show_points Show sample points when Error.bar or statistic is TRUE.
#' @return A ggplot2 object. Test table: attr(p, "stats").
#' @author Hyundong Yoon
#' @export
sc.expr.prop.bar <- function(seurat_obj,
                             target_gene,
                             group_var = "condition",
                             celltype_var = "detailed.celltypes",
                             celltypes = "all",
                             color_pal = NULL,
                             show_legend = TRUE,
                             Error.bar = FALSE,
                             statistic = FALSE,
                             sample_var = "orig.ident",
                             error_type = c("se", "sd"),
                             p_label = c("star", "p"),
                             p_adjust = "BH",
                             min_cells = 10,
                             show_points = TRUE) {

  require(dplyr)
  require(ggplot2)
  require(Seurat)
  error_type <- match.arg(error_type)
  p_label    <- match.arg(p_label)
  sample_mode <- isTRUE(Error.bar) || isTRUE(statistic)

  # 1. Validation
  missing_genes <- target_gene[!target_gene %in% rownames(seurat_obj)]
  if (length(missing_genes) > 0) {
    stop(paste0("Genes not found in object: ", paste(missing_genes, collapse = ", ")))
  }
  meta_data <- seurat_obj@meta.data
  if (sample_mode && !sample_var %in% colnames(meta_data)) {
    stop(paste0("sample_var '", sample_var, "' not found in meta.data (needed for Error.bar / statistic)."))
  }

  # 2. Filter cells
  if (length(celltypes) == 1 && celltypes == "all") {
    valid_cells <- rownames(meta_data)
  } else {
    valid_cells <- rownames(meta_data[meta_data[[celltype_var]] %in% celltypes, ])
    if (length(valid_cells) == 0) stop("No matching cell types found.")
  }

  # 3. Extract data
  raw_counts <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")[target_gene, valid_cells, drop = FALSE]

  grp <- meta_data[valid_cells, group_var]
  grp_levels <- if (is.factor(grp)) levels(droplevels(grp)) else sort(unique(as.character(grp)))
  grp <- as.character(grp)
  smp <- if (sample_mode) as.character(meta_data[valid_cells, sample_var]) else NA_character_

  cell_df <- do.call(rbind, lapply(target_gene, function(g) {
    data.frame(group = grp, sample = smp, gene = g,
               expressed = as.numeric(raw_counts[g, ]) > 0, stringsAsFactors = FALSE)
  })) %>% filter(!is.na(group))

  # 4. Summaries
  sample_df <- NULL
  stat_df <- NULL
  if (!sample_mode) {
    # original behaviour: pooled percentage per group
    bar_df <- cell_df %>%
      group_by(group, gene) %>%
      summarise(pct = sum(expressed) / n() * 100, .groups = "drop") %>%
      mutate(err = 0)
  } else {
    sample_df <- cell_df %>%
      group_by(group, gene, sample) %>%
      summarise(n_cells = n(), pct = sum(expressed) / n_cells * 100, .groups = "drop") %>%
      filter(n_cells >= min_cells)
    bar_df <- sample_df %>%
      group_by(group, gene) %>%
      summarise(n_s = n(), m = mean(pct), s = sd(pct), .groups = "drop") %>%
      mutate(err = if (error_type == "se") s / sqrt(n_s) else s,
             err = ifelse(is.na(err), 0, err),
             pct = m)
  }
  bar_df$group <- factor(bar_df$group, levels = grp_levels)
  bar_df$gene  <- factor(bar_df$gene,  levels = target_gene)
  if (!is.null(sample_df)) {
    sample_df$group <- factor(sample_df$group, levels = grp_levels)
    sample_df$gene  <- factor(sample_df$gene,  levels = target_gene)
  }
  use_err <- isTRUE(Error.bar)
  top_of_bar <- bar_df$pct + if (use_err) bar_df$err else 0

  # 5. Pairwise t-test per gene (on per-sample percentages)
  if (isTRUE(statistic)) {
    K <- length(grp_levels)
    dodge_w <- 0.85
    pair_idx <- combn(seq_len(K), 2)
    stat_df <- do.call(rbind, lapply(seq_along(target_gene), function(i) {
      g <- target_gene[i]
      d <- sample_df[sample_df$gene == g, ]
      ymax_g <- max(c(top_of_bar[bar_df$gene == g],
                      if (show_points) d$pct else 0), na.rm = TRUE)
      out <- do.call(rbind, lapply(seq_len(ncol(pair_idx)), function(k) {
        a <- grp_levels[pair_idx[1, k]]; b <- grp_levels[pair_idx[2, k]]
        xa <- d$pct[d$group == a]; xb <- d$pct[d$group == b]
        p <- if (length(xa) >= 2 && length(xb) >= 2) {
          tryCatch(t.test(xa, xb)$p.value, error = function(e) NA_real_)   # Welch
        } else NA_real_
        data.frame(gene = g, group1 = a, group2 = b, n1 = length(xa), n2 = length(xb),
                   mean1 = mean(xa), mean2 = mean(xb), p = p,
                   x1 = i + (pair_idx[1, k] - (K + 1) / 2) * (dodge_w / K),
                   x2 = i + (pair_idx[2, k] - (K + 1) / 2) * (dodge_w / K),
                   ymax_g = ymax_g, stringsAsFactors = FALSE)
      }))
      out$padj <- p.adjust(out$p, method = p_adjust)
      out
    }))
    stat_df$label <- if (p_label == "star") {
      ifelse(is.na(stat_df$padj), "n.t.",
        ifelse(stat_df$padj < 0.001, "***",
          ifelse(stat_df$padj < 0.01, "**",
            ifelse(stat_df$padj < 0.05, "*", "ns"))))
    } else {
      ifelse(is.na(stat_df$padj), "n.t.", paste0("p = ", signif(stat_df$padj, 2)))
    }
    # stack brackets: narrower pairs lower
    stat_df$span <- abs(stat_df$x2 - stat_df$x1)
    stat_df <- stat_df %>% group_by(gene) %>%
      mutate(rank = rank(span, ties.method = "first"),
             ybar = ymax_g * 1.10 + (rank - 1) * ymax_g * 0.12) %>% ungroup()
    stat_df <- stat_df[!is.na(stat_df$padj) | TRUE, ]
    # drop brackets for untestable pairs (keep in table)
    plot_stat <- stat_df[!is.na(stat_df$padj), ]
  }

  # 6. Visualization
  p <- ggplot(bar_df, aes(x = gene, y = pct, fill = group)) +
    geom_col(color = "black", linewidth = 0.5, width = 0.8,
             position = position_dodge(width = 0.85))

  if (use_err) {
    p <- p + geom_errorbar(aes(ymin = pmax(pct - err, 0), ymax = pct + err),
                           position = position_dodge(width = 0.85), width = 0.25, linewidth = 0.5)
  }
  if (sample_mode && show_points) {
    p <- p + geom_point(data = sample_df, aes(x = gene, y = pct, group = group),
                        inherit.aes = FALSE,
                        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.85),
                        size = 1.5, color = "grey20", alpha = 0.85)
  }

  # numeric annotation (above error bar when shown)
  p <- p + geom_text(aes(y = pct + if (use_err) err else 0, label = sprintf("%.1f%%", pct)),
                     position = position_dodge(width = 0.85),
                     vjust = -0.5, fontface = "bold", size = 3.5)

  n_pairs <- if (isTRUE(statistic)) choose(length(grp_levels), 2) else 0
  p <- p + theme_classic() +
    labs(title = NULL, subtitle = NULL, y = "Percent of Cells (%)", x = "Genes") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15 + (if (isTRUE(statistic)) 0.15 + 0.12 * max(0, n_pairs - 1) else 0)))) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(face = "bold", size = 12),
          legend.title = element_blank())

  if (isTRUE(statistic) && nrow(plot_stat) > 0) {
    tick <- plot_stat$ymax_g * 0.02
    p <- p +
      geom_segment(data = plot_stat, aes(x = x1, xend = x2, y = ybar, yend = ybar), inherit.aes = FALSE, linewidth = 0.5) +
      geom_segment(data = plot_stat, aes(x = x1, xend = x1, y = ybar - tick, yend = ybar), inherit.aes = FALSE, linewidth = 0.5) +
      geom_segment(data = plot_stat, aes(x = x2, xend = x2, y = ybar - tick, yend = ybar), inherit.aes = FALSE, linewidth = 0.5) +
      geom_text(data = plot_stat, aes(x = (x1 + x2) / 2, y = ybar, label = label),
                inherit.aes = FALSE, vjust = -0.3, size = 4)
  }

  if (!is.null(color_pal)) p <- p + scale_fill_manual(values = color_pal)
  if (!show_legend) p <- p + theme(legend.position = "none")

  if (isTRUE(statistic)) attr(p, "stats") <- as.data.frame(stat_df)
  return(p)
}
