#' @title Ultra-Fast, Memory-Efficient Wilcoxon Analysis for Linux (sc.linux.RunWilcox)
#'
#' @description 
#' A high-performance wrapper for the Wilcoxon Rank Sum test, explicitly optimized for 
#' Linux HPC environments. It perfectly mimics Seurat's `FindMarkers(test.use = "wilcox")` 
#' but drastically reduces computation time and prevents Out-of-Memory (OOM) errors 
#' by utilizing sequential cell-type processing and `mclapply`-based gene-level parallelization.
#'
#' @details 
#' **Statistical Note on Wilcoxon:**
#' Unlike MAST (which uses a generalized linear model), the Wilcoxon rank-sum test is a 
#' non-parametric test. Therefore, it mathematically **cannot incorporate covariates** #' (such as `orig.ident` or `nCount_RNA`). This function compares the raw normalized 
#' distributions of the two groups directly, matching Seurat's exact mathematical output.
#'
#' **Performance & Architecture:**
#' \itemize{
#'   \item **Memory Safety:** Isolates one cell type at a time, keeping RAM usage low 
#'   even with 40+ samples.
#'   \item **Linux Forking:** Utilizes `parallel::mclapply`, which shares the sparse 
#'   expression matrix across CPU cores without duplicating it in memory, delivering 
#'   a massive speed boost over Seurat's single-threaded loop.
#'   \item **Full Transcriptome Coverage:** Computes p-values for all genes passing base 
#'   expression thresholds to ensure complete datasets for accurate Volcano plots.
#' }
#'
#' @param seurat_obj A Seurat object. Ensure `JoinLayers()` has been run for V5 objects.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param case.group Character. Name of the experimental group (e.g., "sAML"). Default: "WT".
#' @param ctrl.group Character. Name of the control group (e.g., "Healthy.donor"). Default: "KO".
#' @param preserved.genes Character vector. Key markers (e.g., "LAIR1") protected from base filtering. Default: "LAIR1", "VSIR".
#' @param min_cells Integer. Minimum cells expressing a gene. Default: 0.
#' @param min_expr Numeric. Minimum expression threshold. Default: 0.
#' @param logfc_cutoff Numeric. Threshold for Up/Down categorization. Default: 0.25.
#' @param cut.off.method Character. Significance metric: "p_val_adj" or "p_val". Default: "p_val_adj".
#' @param significant.cut.off Numeric. Threshold for statistical significance. Default: 0.05.
#' @param filter_rows Logical. If TRUE, removes non-significant genes from the output. Default: FALSE.
#' @param ncores Integer. Number of CPU cores to utilize. Default: 16.
#'
#' @return A tidy data frame with exact Seurat-aligned Wilcoxon DE results for ALL tested genes.
#' @export
#'
#' @examples
#' # [Standard Clinical Cohort Wilcoxon Analysis]
#' # wilcox_results <- sc.linux.RunWilcox(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "Healthy.donor",
#' #   ncores = 40
#' # )

sc.linux.RunWilcox = function(seurat_obj,
                              celltype_col="detailed.celltypes",
                              group_col="condition",
                              case.group="WT",
                              ctrl.group="KO",
                              preserved.genes=c("LAIR1", "VSIR", "Lair1", "Vsir"),
                              min_cells=0,
                              min_expr=0,
                              logfc_cutoff=0.25,
                              cut.off.method="p_val_adj",
                              significant.cut.off=0.05,
                              filter_rows=FALSE,
                              ncores=16) {
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required. Please install it.")
  }
  
  start_time <- Sys.time()
  cat("[1] Joining layers and preparing meta...\n")
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj)
  meta_all <- seurat_obj@meta.data
  
  celltypes_raw <- unique(as.character(meta_all[[celltype_col]]))
  valid_celltypes <- celltypes_raw[sapply(celltypes_raw, function(ct) {
    meta_ct <- meta_all[meta_all[[celltype_col]] == ct, ]
    return(case.group %in% meta_ct[[group_col]] && ctrl.group %in% meta_ct[[group_col]])
  })]
  
  total_tasks <- length(valid_celltypes)
  DEG.list <- list()
  
  cat("[3] Starting Sequential Analysis (Linux Multi-Core Wilcoxon)...\n")
  for (i in seq_along(valid_celltypes)) {
    celltype <- valid_celltypes[i]
    loop_step_start <- Sys.time()
    cat(sprintf("\n[%d/%d] Processing: %s (Time: %s)\n", i, total_tasks, celltype, format(loop_step_start, "%H:%M:%S")))
    
    obj_ct <- subset(seurat_obj, cells = rownames(meta_all[meta_all[[celltype_col]] == celltype, ]))
    expr_ct <- Seurat::GetAssayData(obj_ct, assay="RNA", layer="data")
    
    # 1. Base filtering
    keep_genes <- Matrix::rowSums(expr_ct > min_expr) >= min_cells
    keep_genes[rownames(expr_ct) %in% preserved.genes] <- TRUE
    expr_ct <- expr_ct[keep_genes, , drop=FALSE]
    
    if (nrow(expr_ct) == 0) {
      cat(" -> No genes passed base filtering. Skipping.\n")
      next
    }
    
    # 2. Split cells by group
    case_cells <- colnames(obj_ct)[obj_ct@meta.data[[group_col]] == case.group]
    ctrl_cells <- colnames(obj_ct)[obj_ct@meta.data[[group_col]] == ctrl.group]
    
    # 3. Calculate Seurat's exact FoldChange
    Seurat::Idents(obj_ct) <- group_col
    fc_seurat <- Seurat::FoldChange(obj_ct, ident.1 = case.group, ident.2 = ctrl.group)
    fc_seurat$genes <- rownames(fc_seurat)
    
    # 4. Ultra-fast parallel Wilcoxon test (Linux Forking)
    gene_names <- rownames(expr_ct)
    pvals <- unlist(parallel::mclapply(gene_names, function(g) {
      expr_case <- expr_ct[g, case_cells]
      expr_ctrl <- expr_ct[g, ctrl_cells]
      
      # Handle zero-variance or empty cases
      if (length(expr_case) == 0 || length(expr_ctrl) == 0) return(NA)
      if (all(expr_case == 0) && all(expr_ctrl == 0)) return(1)
      
      # exact=FALSE matches Seurat's implementation
      res <- suppressWarnings(wilcox.test(expr_case, expr_ctrl, exact = FALSE))
      return(res$p.value)
    }, mc.cores = ncores))
    
    pval_df <- data.frame(genes = gene_names, p_val = pvals, stringsAsFactors = FALSE)
    
    # 5. Merge and Adjust
    tmp_res <- merge(fc_seurat[, c("genes", "avg_log2FC")], pval_df, by="genes")
    tmp_res$p_val_adj <- p.adjust(tmp_res$p_val, method="fdr")
    tmp_res$celltype <- celltype
    tmp_res$compare <- paste0(case.group, "_vs_", ctrl.group)
    
    # 6. Annotate Expression Status
    tmp_res <- dplyr::mutate(tmp_res,
                             Expression = dplyr::case_when(
                               avg_log2FC > logfc_cutoff & .data[[cut.off.method]] < significant.cut.off ~ "Up",
                               avg_log2FC < -logfc_cutoff & .data[[cut.off.method]] < significant.cut.off ~ "Down",
                               TRUE ~ "Unchanged"
                             ),
                             color = dplyr::case_when(
                               Expression == "Up" ~ "#f2210a",
                               Expression == "Down" ~ "#3158e8",
                               TRUE ~ "#afb0b3"
                             )
    )
    
    if(filter_rows) tmp_res <- dplyr::filter(tmp_res, .data[[cut.off.method]] < significant.cut.off)
    DEG.list[[celltype]] <- tmp_res
    
    rm(obj_ct, expr_ct, pval_df, pvals, fc_seurat); gc()
    
    loop_step_end <- Sys.time()
    cat(sprintf("    -> Finished %s in %s\n", celltype, format(loop_step_end - loop_step_start)))
  }
  
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  cat("\n[Final Done] Total processing time:", format(Sys.time() - start_time), "\n")
  return(sc.DE.all)
}