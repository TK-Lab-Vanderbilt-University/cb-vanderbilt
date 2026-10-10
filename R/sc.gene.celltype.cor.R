#' Correlate a gene's sample-level expression with celltype fractions
#'
#' For one group of samples (e.g. Tumor only), computes per-sample celltype
#' fractions and a per-sample gene expression value, then runs a Spearman (or
#' Pearson/Kendall) correlation for every celltype. The unit of analysis is the
#' SAMPLE (not the cell), so p-values are not inflated by cell number.
#'
#' Filtering is minimal by default: no sample is dropped for low cell number and
#' a celltype is tested whenever at least 3 samples have a usable value and both
#' variables are not constant. Raise `min_cells_sample` / `min_sample_present`
#' for stricter behaviour.
#'
#' @param seurat_obj Seurat object.
#' @param gene Gene symbol (default "VSIR").
#' @param group_var,group_value Metadata column / value used to subset
#'   samples (default condition == "Tumor"). Set group_var = NULL to use all cells.
#' @param celltype_var Metadata column with celltype labels.
#' @param sample_var Metadata column with sample IDs (the replicate unit).
#' @param gene_def How the per-sample gene value is defined:
#'   "excl" = mean over all cells EXCEPT the celltype being tested (avoids
#'            circularity when the gene is expressed by that celltype);
#'   "epi"  = mean over `epi_celltypes` only (tumor/epithelial-cell expression);
#'   "all"  = mean over all cells (circular for celltypes expressing the gene).
#' @param epi_celltypes Celltypes used for gene_def = "epi".
#' @param denom Denominator for fractions: "all" = % of all cells,
#'   "immune" = % of immune cells (only `immune_celltypes` are tested).
#' @param immune_celltypes Celltypes counted as immune for denom = "immune".
#' @param min_cells_sample Drop samples with fewer cells than this (default 0 = keep all).
#' @param min_sample_present A celltype is tested only if it is present
#'   (count > 0) in at least this many samples (default 3 = minimum for a correlation).
#' @param method "spearman", "pearson" or "kendall".
#' @param p_adjust Method passed to p.adjust (default "BH").
#' @param plot Draw rho bar plot and scatter plots?
#' @param scatter_p Draw scatter plots for celltypes with raw p < scatter_p
#'   (default 1 = draw for every tested celltype).
#' @param label_samples Label points with sample names in scatter plots.
#' @param smooth Trend line: "loess" (wavy, default), "lm" (straight) or "none".
#' @param smooth_span loess span; smaller = more wiggly (default 0.9, 0.5-1 sensible for ~8 samples).
#' @param smooth_se Draw the confidence band around the trend line?
#' @param point_size Dot size in scatter plots (default 5).
#' @param point_color Dot fill colour. If NULL, dots are coloured by `color_by`.
#' @param color_by Optional named vector (names = sample IDs) or NULL; used when point_color is NULL.
#' @param line_color Trend line colour.
#' @param band_color Confidence band colour.
#' @param base_size Base font size for plots.
#' @param verbose Print progress messages.
#'
#' @return A list (class "sc_gene_celltype_cor") with `stats` (data.frame),
#'   `fractions` (sample x celltype %), `gene_expr` (per-sample gene value[s]),
#'   `plots` (list of ggplots), `params`.
#' @export
sc.gene.celltype.cor <- function(
    seurat_obj,
    gene = "VSIR",
    group_var = "condition",
    group_value = "Tumor",
    celltype_var = "detailed.celltypes",
    sample_var = "orig.ident",
    gene_def = c("excl", "epi", "all"),
    epi_celltypes = c("Cholangiocyte", "Hepatocyte", "Cycling.cell"),
    denom = c("all", "immune"),
    immune_celltypes = NULL,
    min_cells_sample = 0,
    min_sample_present = 3,
    method = c("spearman", "pearson", "kendall"),
    p_adjust = "BH",
    plot = TRUE,
    scatter_p = 1,
    label_samples = TRUE,
    smooth = c("loess", "lm", "none"),
    smooth_span = 0.9,
    smooth_se = TRUE,
    point_size = 5,
    point_color = "#E64B35",
    color_by = NULL,
    line_color = "#3C5488",
    band_color = "#B8C4E0",
    base_size = 13,
    verbose = TRUE) {

  gene_def <- match.arg(gene_def)
  denom    <- match.arg(denom)
  method   <- match.arg(method)
  smooth   <- match.arg(smooth)

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required for plot = TRUE.")

  md <- seurat_obj@meta.data
  for (v in c(celltype_var, sample_var)) {
    if (!v %in% colnames(md)) stop("Column not found in meta.data: ", v)
  }
  if (!gene %in% rownames(seurat_obj)) stop("Gene not found: ", gene)

  # ---- subset cells --------------------------------------------------------
  keep <- rep(TRUE, nrow(md))
  if (!is.null(group_var)) {
    if (!group_var %in% colnames(md)) stop("Column not found in meta.data: ", group_var)
    keep <- as.character(md[[group_var]]) %in% as.character(group_value)
  }
  cells <- rownames(md)[keep]
  if (length(cells) == 0) stop("No cells left after subsetting.")
  md <- md[cells, , drop = FALSE]
  expr <- as.numeric(Seurat::FetchData(seurat_obj, vars = gene, cells = cells)[, 1])
  ct <- as.character(md[[celltype_var]])
  sm <- as.character(md[[sample_var]])

  # ---- drop small samples (default: keep all) -------------------------------
  ok_s <- names(which(table(sm) >= min_cells_sample))
  sel  <- sm %in% ok_s
  expr <- expr[sel]; ct <- ct[sel]; sm <- sm[sel]
  samples <- sort(unique(sm))
  if (verbose) message("Samples used (n = ", length(samples), "): ", paste(samples, collapse = ", "))
  if (length(samples) < 3) stop("Fewer than 3 samples: correlation cannot be computed.")
  if (length(samples) < 6) warning("Only ", length(samples), " samples: correlations are very unstable; interpret rho, not p.")

  # ---- fractions (sample x celltype, %) ------------------------------------
  cnt <- table(factor(sm, levels = samples), ct)
  if (denom == "immune") {
    if (is.null(immune_celltypes)) stop("denom = 'immune' requires immune_celltypes.")
    cnt <- cnt[, colnames(cnt) %in% immune_celltypes, drop = FALSE]
  }
  frac <- unclass(cnt) / pmax(rowSums(cnt), 1) * 100

  # ---- per-sample gene values ----------------------------------------------
  sf    <- factor(sm, levels = samples)
  n_all <- tapply(expr, sf, length)
  s_all <- tapply(expr, sf, sum)
  n_ct  <- sapply(colnames(cnt), function(x) tapply(ct == x, sf, sum))
  s_ct  <- sapply(colnames(cnt), function(x) tapply(ifelse(ct == x, expr, 0), sf, sum))

  gene_for <- function(cell_type) {
    if (gene_def == "all") {
      s_all / n_all
    } else if (gene_def == "excl") {
      den <- n_all - n_ct[, cell_type]
      ifelse(den > 0, (s_all - s_ct[, cell_type]) / den, NA_real_)
    } else {
      is_e <- ct %in% epi_celltypes
      m <- tapply(ifelse(is_e, expr, NA), sf, mean, na.rm = TRUE)
      m[is.nan(m)] <- NA_real_
      m
    }
  }

  # ---- correlation per celltype --------------------------------------------
  rows <- lapply(colnames(frac), function(x) {
    g <- gene_for(x); y <- frac[, x]
    ok <- !is.na(g) & !is.na(y)
    n_present <- sum(cnt[, x] > 0)
    if (sum(ok) < 3 || n_present < min_sample_present ||
        stats::sd(y[ok]) == 0 || stats::sd(g[ok]) == 0) return(NULL)
    ct_res <- suppressWarnings(stats::cor.test(g[ok], y[ok], method = method))
    data.frame(celltype = x, n_samples = sum(ok), n_present = n_present,
               mean_pct = mean(y[ok]), rho = unname(ct_res$estimate),
               p = ct_res$p.value, stringsAsFactors = FALSE)
  })
  stats_df <- do.call(rbind, rows)
  if (is.null(stats_df)) stop("No celltype could be tested (too few samples or constant values).")
  stats_df$padj <- stats::p.adjust(stats_df$p, p_adjust)
  stats_df$sig  <- ifelse(stats_df$padj < 0.05, "FDR<0.05", ifelse(stats_df$p < 0.05, "p<0.05", "ns"))
  stats_df <- stats_df[order(-stats_df$rho), ]
  rownames(stats_df) <- NULL

  gene_tab <- if (gene_def == "excl") {
    sapply(colnames(frac), gene_for)
  } else data.frame(sample = samples, value = as.numeric(gene_for(NULL)))

  # ---- plots ---------------------------------------------------------------
  plots <- list()
  if (plot) {
    sub_txt <- paste0(if (!is.null(group_var)) paste0(group_var, " = ", paste(group_value, collapse = "/"), "  |  ") else "",
                      "sample-level, n = ", length(samples), "  |  gene: ", gene_def, ", denom: ", denom)

    # (1) rho bar plot
    d <- stats_df
    d$celltype <- factor(d$celltype, levels = d$celltype[order(d$rho)])
    d$star <- ifelse(d$padj < 0.05, "**", ifelse(d$p < 0.05, "*", ""))
    xr <- max(abs(d$rho)) * 1.15
    plots$rho <- ggplot2::ggplot(d, ggplot2::aes(rho, celltype, fill = rho)) +
      ggplot2::geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
      ggplot2::geom_text(ggplot2::aes(label = star, hjust = ifelse(rho >= 0, -0.3, 1.3)),
                         size = 5, vjust = 0.75, color = "grey15") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey25") +
      ggplot2::scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#DC0000",
                                    midpoint = 0, limits = c(-1, 1), name = "rho") +
      ggplot2::scale_x_continuous(limits = c(-xr, xr), expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::labs(x = paste0(method, " rho"), y = NULL,
                    title = paste0(gene, " expression vs celltype fraction"),
                    subtitle = sub_txt, caption = "* p < 0.05   ** FDR < 0.05") +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(color = "grey35", size = base_size - 3),
                     axis.text.y = ggplot2::element_text(color = "black"),
                     legend.position = "right")

    # (2) scatter plots
    top <- stats_df$celltype[stats_df$p < scatter_p]
    use_repel <- requireNamespace("ggrepel", quietly = TRUE)
    plots$scatter <- lapply(stats::setNames(top, top), function(x) {
      dd <- data.frame(sample = samples, g = as.numeric(gene_for(x)), y = frac[, x])
      dd <- dd[stats::complete.cases(dd), ]
      r  <- stats_df[stats_df$celltype == x, ]

      p <- ggplot2::ggplot(dd, ggplot2::aes(g, y))
      if (smooth != "none" && nrow(dd) >= 4) {
        sm_args <- list(method = smooth, formula = y ~ x, se = smooth_se,
                        color = line_color, fill = band_color, linewidth = 1.3, alpha = 0.35)
        if (smooth == "loess") sm_args$span <- smooth_span
        p <- p + do.call(ggplot2::geom_smooth, sm_args)
      }
      if (is.null(point_color) && !is.null(color_by)) {
        dd$grp <- color_by[dd$sample]
        p <- p + ggplot2::geom_point(data = dd, ggplot2::aes(fill = grp), shape = 21,
                                     size = point_size, color = "white", stroke = 1.2)
      } else {
        p <- p + ggplot2::geom_point(shape = 21, size = point_size, fill = point_color %||% "#E64B35",
                                     color = "white", stroke = 1.2)
      }
      if (label_samples) {
        p <- if (use_repel) p + ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3.2,
                                                         color = "grey25", max.overlaps = Inf,
                                                         box.padding = 0.5, segment.color = "grey70")
             else p + ggplot2::geom_text(ggplot2::aes(label = sample), size = 3.2, vjust = -1.4,
                                         color = "grey25", check_overlap = TRUE)
      }
      p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
        ggplot2::labs(x = paste0(gene, " expression (mean log-norm; ", gene_def, ")"),
                      y = paste0(x, " (% of ", denom, " cells)"),
                      title = x,
                      subtitle = sprintf("%s rho = %.2f   p = %.3g   FDR = %.3g   (n = %d)",
                                         method, r$rho, r$p, r$padj, r$n_samples)) +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                       plot.subtitle = ggplot2::element_text(color = "grey35", size = base_size - 3),
                       panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
                       axis.line = ggplot2::element_line(color = "grey30"))
    })
  }

  structure(list(stats = stats_df, fractions = frac, gene_expr = gene_tab, plots = plots,
                 params = list(gene = gene, group_var = group_var, group_value = group_value,
                               gene_def = gene_def, denom = denom, method = method,
                               smooth = smooth, n_samples = length(samples), samples = samples)),
            class = "sc_gene_celltype_cor")
}

# small helper (NULL-coalescing) kept internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.sc_gene_celltype_cor <- function(x, ...) {
  cat("sc.gene.celltype.cor: ", x$params$gene, " | n samples = ", x$params$n_samples,
      " | gene_def = ", x$params$gene_def, " | denom = ", x$params$denom, "\n", sep = "")
  print(utils::head(x$stats, 10), row.names = FALSE)
  invisible(x)
}
