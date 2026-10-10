#' Correlate a gene's sample-level expression with celltype fractions
#'
#' For one group of samples (e.g. Tumor only), computes per-sample celltype
#' fractions and a per-sample gene expression value, then runs a Spearman (or
#' Pearson/Kendall) correlation for every celltype. The unit of analysis is the
#' SAMPLE (not the cell), so p-values are not inflated by cell number.
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
#' @param min_cells_sample Drop samples with fewer cells than this.
#' @param min_sample_present A celltype is tested only if it is present
#'   (count > 0) in at least this many samples.
#' @param method "spearman", "pearson" or "kendall".
#' @param p_adjust Method passed to p.adjust (default "BH").
#' @param plot Draw rho bar plot and scatter plots?
#' @param scatter_p Draw scatter plots for celltypes with raw p < scatter_p.
#' @param label_samples Label points with sample names in scatter plots.
#' @param verbose Print progress messages.
#'
#' @return A list with `stats` (data.frame), `fractions` (sample x celltype %),
#'   `gene_expr` (per-sample gene value[s]), `plots` (list of ggplots), `params`.
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
    min_cells_sample = 200,
    min_sample_present = 4,
    method = c("spearman", "pearson", "kendall"),
    p_adjust = "BH",
    plot = TRUE,
    scatter_p = 0.1,
    label_samples = TRUE,
    verbose = TRUE) {

  gene_def <- match.arg(gene_def)
  denom    <- match.arg(denom)
  method   <- match.arg(method)

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

  # ---- drop small samples --------------------------------------------------
  ok_s <- names(which(table(sm) >= min_cells_sample))
  sel  <- sm %in% ok_s
  expr <- expr[sel]; ct <- ct[sel]; sm <- sm[sel]
  samples <- sort(unique(sm))
  if (verbose) message("Samples used (n = ", length(samples), "): ", paste(samples, collapse = ", "))
  if (length(samples) < 4) warning("Fewer than 4 samples: correlation is not meaningful.")

  # ---- fractions (sample x celltype, %) ------------------------------------
  cnt <- table(factor(sm, levels = samples), ct)
  if (denom == "immune") {
    if (is.null(immune_celltypes)) stop("denom = 'immune' requires immune_celltypes.")
    cnt <- cnt[, colnames(cnt) %in% immune_celltypes, drop = FALSE]
  }
  frac <- unclass(cnt) / pmax(rowSums(cnt), 1) * 100

  # ---- per-sample gene values ----------------------------------------------
  n_all <- tapply(expr, factor(sm, levels = samples), length)
  s_all <- tapply(expr, factor(sm, levels = samples), sum)
  n_ct  <- sapply(colnames(cnt), function(x) tapply(ct == x, factor(sm, levels = samples), sum))
  s_ct  <- sapply(colnames(cnt), function(x) tapply(ifelse(ct == x, expr, 0), factor(sm, levels = samples), sum))

  gene_for <- function(cell_type) {
    if (gene_def == "all") {
      s_all / n_all
    } else if (gene_def == "excl") {
      den <- n_all - n_ct[, cell_type]
      ifelse(den > 0, (s_all - s_ct[, cell_type]) / den, NA_real_)
    } else {
      is_e <- ct %in% epi_celltypes
      m <- tapply(ifelse(is_e, expr, NA), factor(sm, levels = samples), mean, na.rm = TRUE)
      m[is.nan(m)] <- NA_real_
      m
    }
  }

  # ---- correlation per celltype --------------------------------------------
  rows <- lapply(colnames(frac), function(x) {
    g <- gene_for(x); y <- frac[, x]
    ok <- !is.na(g) & !is.na(y)
    n_present <- sum(cnt[, x] > 0)
    if (sum(ok) < 4 || n_present < min_sample_present || sd(y[ok]) == 0 || sd(g[ok]) == 0) return(NULL)
    ct_res <- suppressWarnings(stats::cor.test(g[ok], y[ok], method = method))
    data.frame(celltype = x, n_samples = sum(ok), n_present = n_present,
               mean_pct = mean(y[ok]), rho = unname(ct_res$estimate),
               p = ct_res$p.value, stringsAsFactors = FALSE)
  })
  stats_df <- do.call(rbind, rows)
  if (is.null(stats_df)) stop("No celltype passed the filters (min_sample_present / variance).")
  stats_df$padj <- stats::p.adjust(stats_df$p, p_adjust)
  stats_df$sig <- ifelse(stats_df$padj < 0.05, "FDR<0.05", ifelse(stats_df$p < 0.05, "p<0.05", "ns"))
  stats_df <- stats_df[order(-stats_df$rho), ]
  rownames(stats_df) <- NULL

  gene_tab <- if (gene_def == "excl") {
    sapply(colnames(frac), gene_for)
  } else data.frame(sample = samples, value = as.numeric(gene_for(NULL)))

  # ---- plots ---------------------------------------------------------------
  plots <- list()
  if (plot) {
    d <- stats_df
    d$celltype <- factor(d$celltype, levels = d$celltype[order(d$rho)])
    plots$rho <- ggplot2::ggplot(d, ggplot2::aes(rho, celltype, fill = sig)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "p<0.05" = "#F4A582", ns = "grey80"), name = NULL) +
      ggplot2::labs(x = paste0(method, " rho (", gene, " vs celltype fraction)"), y = NULL,
                    title = paste0(gene, " vs celltype fraction"),
                    subtitle = paste0(if (!is.null(group_var)) paste0(group_var, " = ", paste(group_value, collapse = "/"), ", ") else "",
                                      "sample-level, n = ", length(samples), "; gene_def = ", gene_def, ", denom = ", denom)) +
      ggplot2::theme_classic(base_size = 11)

    top <- stats_df$celltype[stats_df$p < scatter_p]
    plots$scatter <- lapply(stats::setNames(top, top), function(x) {
      dd <- data.frame(sample = samples, g = as.numeric(gene_for(x)), y = frac[, x])
      r <- stats_df[stats_df$celltype == x, ]
      p <- ggplot2::ggplot(dd, ggplot2::aes(g, y)) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "grey30", fill = "grey85", linewidth = 0.6) +
        ggplot2::geom_point(size = 3, color = "#E64B35") +
        ggplot2::labs(x = paste0(gene, " (mean log-norm; ", gene_def, ")"),
                      y = paste0(x, " (% of ", denom, " cells)"),
                      subtitle = sprintf("%s rho = %.2f, p = %.3g, FDR = %.3g", method, r$rho, r$p, r$padj)) +
        ggplot2::theme_classic(base_size = 11)
      if (label_samples) p <- p + ggplot2::geom_text(ggplot2::aes(label = sample), size = 2.5, vjust = -1, check_overlap = TRUE)
      p
    })
  }

  structure(list(stats = stats_df, fractions = frac, gene_expr = gene_tab, plots = plots,
                 params = list(gene = gene, group_var = group_var, group_value = group_value,
                               gene_def = gene_def, denom = denom, method = method,
                               n_samples = length(samples), samples = samples)),
            class = "sc_gene_celltype_cor")
}

#' @export
print.sc_gene_celltype_cor <- function(x, ...) {
  cat("sc.gene.celltype.cor: ", x$params$gene, " | n samples = ", x$params$n_samples,
      " | gene_def = ", x$params$gene_def, " | denom = ", x$params$denom, "\n", sep = "")
  print(utils::head(x$stats, 10), row.names = FALSE)
  invisible(x)
}
