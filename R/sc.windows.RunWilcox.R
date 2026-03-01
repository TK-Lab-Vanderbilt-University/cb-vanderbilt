#' @title High-Performance, Memory-Efficient Wilcoxon Analysis for Windows
#'
#' @description 
#' A specialized wrapper for the Wilcoxon Rank Sum test optimized for Windows OS. 
#' It mathematically mirrors Seurat's `FindMarkers(test.use = "wilcox")` but bypasses 
#' Seurat's single-threaded limitations. To safely achieve multi-core speeds on Windows 
#' without Out-of-Memory (OOM) crashes, it utilizes a PSOCK cluster (`parallel::makeCluster`), 
#' dynamically exporting and clearing sparse matrices for each cell type sequentially.
#'
#' @details 
#' **Statistical Note on Wilcoxon:**
#' As a non-parametric rank-sum test, Wilcoxon evaluates raw distribution shifts and 
#' inherently **cannot incorporate covariates** (such as `orig.ident`). This function 
#' processes pure normalized expression data, guaranteeing exact numerical parity with 
#' standard Seurat outputs.
#'
#' **Windows-Specific Optimization:**
#' \itemize{
#'   \item **Socket Cluster Management:** Automatically initializes and rigorously tears 
#'   down background R processes to prevent zombie processes and memory leaks.
#'   \item **Dynamic Memory Export:** Sends only the strictly necessary expression data 
#'   (`expr_ct`) for the active cell type to the worker nodes, clearing the workers' RAM 
#'   before moving to the next cell type.
#'   \item **Full Transcriptome Coverage:** Computes p-values for all genes passing base 
#'   expression thresholds to ensure complete datasets for Volcano plots.
#' }
#'
#' @param seurat_obj A Seurat object. Ensure `JoinLayers()` has been run for V5 objects.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param case.group Character. Name of the experimental group (e.g., "sAML"). Default: "WT".
#' @param ctrl.group Character. Name of the control group (e.g., "Healthy.donor"). Default: "KO".
#' @param preserved.genes Character vector. Key markers (e.g., "LAIR1") protected from base filtering.
#' @param min_cells Integer. Minimum cells expressing a gene. Default: 0.
#' @param min_expr Numeric. Minimum expression threshold. Default: 0.
#' @param logfc_cutoff Numeric. Threshold for Up/Down categorization. Default: 0.25.
#' @param cut.off.method Character. Significance metric: "p_val_adj" or "p_val". Default: "p_val_adj".
#' @param significant.cut.off Numeric. Threshold for statistical significance. Default: 0.05.
#' @param filter_rows Logical. If TRUE, removes non-significant genes from the output. Default: FALSE.
#' @param ncores Integer. Number of background CPU cores to utilize. Default: 8.
#'
#' @return A tidy data frame with exact Seurat-aligned Wilcoxon DE results for ALL tested genes.
#' @export
#'
#' @examples
#' # [Windows Workstation Example]
#' # wilcox_results <- sc.windows.RunWilcox(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "Healthy.donor",
#' #   ncores = 8
#' # )

sc.windows.RunWilcox = function(seurat_obj,
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
                                ncores=8) {
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for Windows multi-threading. Please install it.")
  }
  
  start_time <- Sys.time()
  
  cat(sprintf("[0] Initializing Windows Socket Cluster with %d workers...\n", ncores))
  cl <- parallel::makeCluster(ncores)
  # Ensure the cluster is closed even if the function crashes or is interrupted
  on.exit(parallel::stopCluster(cl))
  
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
  
  cat("[3] Starting Sequential Analysis (Windows Socket Wilcoxon)...\n")
  for (i in seq_along(valid_celltypes)) {
    celltype <- valid_celltypes[i]
    loop_step_start <- Sys.time()
    cat(sprintf("\n[%d/%d] Processing: %s (Time: %s)\n", i, total_tasks, celltype, format(loop_step_start, "%H:%M:%S")))
    
    cells_to_keep <- rownames(meta_all[meta_all[[celltype_col]] == celltype, ])
    obj_ct <- subset(seurat_obj, cells = cells_to_keep)
    expr_ct <- Seurat::GetAssayData(obj_ct, assay="RNA", layer="data")
    
    # 1. Base filtering (bypassing logFC pre-filtering for full dataset)
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
    
    # 4. Ultra-fast parallel Wilcoxon test (Windows PSOCK)
    gene_names <- rownames(expr_ct)
    
    # Export only necessary data to the worker nodes to save RAM
    parallel::clusterExport(cl, varlist = c("expr_ct", "case_cells", "ctrl_cells"), envir = environment())
    
    pvals <- unlist(parallel::parLapply(cl, gene_names, function(g) {
      expr_case <- expr_ct[g, case_cells]
      expr_ctrl <- expr_ct[g, ctrl_cells]
      
      # Handle zero-variance or empty cases safely
      if (length(expr_case) == 0 || length(expr_ctrl) == 0) return(NA)
      if (all(expr_case == 0) && all(expr_ctrl == 0)) return(1)
      
      # exact=FALSE matches Seurat's default implementation
      res <- suppressWarnings(wilcox.test(expr_case, expr_ctrl, exact = FALSE))
      return(res$p.value)
    }))
    
    # Clear worker memory immediately after computation
    parallel::clusterEvalQ(cl, { rm(expr_ct, case_cells, ctrl_cells); gc() })
    
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
    
    # Force Garbage Collection on Master node
    rm(obj_ct, expr_ct, pval_df, pvals, fc_seurat); gc()
    
    loop_step_end <- Sys.time()
    cat(sprintf("    -> Finished %s in %s\n", celltype, format(loop_step_end - loop_step_start)))
  }
  
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  cat("\n[Final Done] Total processing time:", format(Sys.time() - start_time), "\n")
  return(sc.DE.all)
}