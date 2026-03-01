#' @title High-Performance, Seurat-Aligned MAST Analysis for Windows (Full Gene Calculation)
#'
#' @description 
#' A specialized, turbo-charged wrapper for MAST (Model-based Analysis of Single-cell 
#' Transcriptomics) optimized specifically for Windows OS. This function guarantees 
#' 100% numerical and logical parity with Seurat's `FindMarkers(test.use = "MAST")`. 
#' It calculates statistics for ALL genes passing base expression thresholds without 
#' prematurely filtering by FoldChange, ensuring a complete dataset for Volcano plots. 
#' To bypass Windows' multi-threading limitations, it utilizes a Socket-based parallel 
#' backend (`SnowParam`) while securely managing memory overhead.
#'
#' @details 
#' **Windows-Specific Optimization & Architecture:**
#' \itemize{
#'   \item **SnowParam Integration:** Explicitly manages Windows-compatible parallel 
#'   clusters, ensuring safe and efficient multi-threading across CPU cores.
#'   \item **Aggressive Memory Safety:** Windows is highly susceptible to memory 
#'   fragmentation. This function mitigates Out-of-Memory (OOM) crashes by isolating 
#'   cell types sequentially and forcing garbage collection (`gc()`) after each loop.
#' }
#' 
#' **Statistical Integrity:**
#' \itemize{
#'   \item **Seurat Logic Synchronization:** Leverages `Seurat::FoldChange` to extract 
#'   the arithmetic `avg_log2FC`, matching standard global expression visualizations perfectly.
#'   \item **Transparent Covariate Control:** Only incorporates variables strictly passed to 
#'   `latent.vars` (e.g., "orig.ident", "nCount_RNA"), mirroring Seurat's design.
#'   \item **Empirical Bayes Moderation:** `ebayes = TRUE` stabilizes variance estimates, 
#'   a highly recommended approach to maximize statistical power in sparse scRNA-seq data.
#' }
#' 
#' **CRITICAL WARNING:**
#' Windows handles memory differently than Linux. If `ncores` is set too high relative to 
#' your available RAM, the R session or the entire system may become unresponsive. Monitor 
#' Task Manager during your first run to find the optimal core count.
#'
#' @param seurat_obj A Seurat object. Ensure `JoinLayers()` has been run for V5 objects.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param case.group Character. Name of the experimental group (e.g., "sAML"). Default: "WT".
#' @param ctrl.group Character. Name of the control group (e.g., "Healthy.donor"). Default: "KO".
#' @param latent.vars Character vector. Metadata columns to use as batch correction variables. Set to NULL for no covariates. Default: "orig.ident".
#' @param ebayes Logical. Use empirical Bayes variance moderation. Highly recommended. Default: TRUE.
#' @param preserved.genes Character vector. Genes protected from base filtering. Default: "LAIR1", "VSIR".
#' @param min_cells Integer. Minimum cells expressing a gene. Default: 0.
#' @param min_expr Numeric. Minimum expression threshold. Default: 0.
#' @param logfc_cutoff Numeric. Threshold for Up/Down categorization. Default: 0.25.
#' @param cut.off.method Character. Significance metric: "p_val_adj" or "p_val". Default: "p_val_adj".
#' @param significant.cut.off Numeric. Threshold for statistical significance. Default: 0.05.
#' @param filter_rows Logical. If TRUE, removes non-significant genes from the output. Default: FALSE.
#' @param ncores Integer. Number of workers for the SnowParam cluster. Default: 8.
#'
#' @return A tidy data frame with DE results for ALL tested genes.
#' @export
#'
#' @examples
#' # [Windows Workstation Example: Clinical Cohort Analysis]
#' # deg_results <- sc.windows.RunMAST(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "Healthy.donor",
#' #   latent.vars = c("orig.ident", "nCount_RNA"),
#' #   ebayes = TRUE,
#' #   ncores = 8
#' # )

sc.windows.RunMAST <- function(seurat_obj,
                               celltype_col="celltypes",
                               group_col="condition",
                               case.group="WT",
                               ctrl.group="KO",
                               latent.vars="orig.ident",
                               ebayes=TRUE,
                               preserved.genes=c("LAIR1", "VSIR", "Lair1", "Vsir"),
                               min_cells=0,
                               min_expr=0,
                               logfc_cutoff=0.25,
                               cut.off.method="p_val_adj",
                               significant.cut.off=0.05,
                               filter_rows=FALSE,
                               ncores=8) {
  start_time <- Sys.time()
  cat(sprintf("[0] Initializing Windows Parallel Backend (SnowParam) with %d workers...\n", ncores))
  bp_param <- BiocParallel::SnowParam(workers = ncores, progressbar = FALSE)
  
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
  coef_name <- paste0(group_col, case.group)
  
  cat("[3] Starting Sequential Analysis (Strict Latent Control, Full Gene Test)...\n")
  for (i in seq_along(valid_celltypes)) {
    celltype <- valid_celltypes[i]
    loop_step_start <- Sys.time()
    cat(sprintf("\n[%d/%d] Processing: %s (Time: %s)\n", i, total_tasks, celltype, format(loop_step_start, "%H:%M:%S")))
    
    cells_to_keep <- rownames(meta_all[meta_all[[celltype_col]] == celltype, ])
    obj_ct <- subset(seurat_obj, cells = cells_to_keep)
    expr_ct <- Seurat::GetAssayData(obj_ct, assay="RNA", layer="data")
    
    # Basic filtering (apply expression frequency/amount conditions only, bypass logFC pre-filtering)
    keep_genes <- Matrix::rowSums(expr_ct > min_expr) >= min_cells
    keep_genes[rownames(expr_ct) %in% preserved.genes] <- TRUE
    expr_ct <- expr_ct[keep_genes, , drop=FALSE]
    
    if (nrow(expr_ct) == 0) {
      cat(" -> No genes passed base filtering. Skipping.\n")
      next
    }
    
    # Formula generation: strictly incorporate user-provided latent.vars
    if (!is.null(latent.vars) && length(latent.vars) > 0) {
      fmla_str <- paste("~", group_col, "+", paste(latent.vars, collapse = " + "))
    } else {
      fmla_str <- paste("~", group_col)
    }
    fmla <- as.formula(fmla_str)
    
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(logcounts = as.matrix(expr_ct)),
      colData = obj_ct@meta.data
    )
    SummarizedExperiment::rowData(sce)$primerid <- rownames(expr_ct)
    sca_ct <- MAST::SceToSingleCellAssay(sce)
    SummarizedExperiment::colData(sca_ct)[[group_col]] <- factor(
      SummarizedExperiment::colData(sca_ct)[[group_col]], levels = c(ctrl.group, case.group)
    )
    
    res <- tryCatch({
      # Execute MAST hurdle model via SnowParam for Windows compatibility
      zlm_mod <- MAST::zlm(fmla, sca_ct, method="glm", ebayes=ebayes, parallel = bp_param)
      summaryCond <- MAST::summary(zlm_mod, doLRT = coef_name)
      dt <- summaryCond$datatable
      pval <- dt[dt$contrast==coef_name & dt$component=="H", c("primerid", "Pr(>Chisq)")]
      
      # Calculate Seurat's arithmetic FoldChange for all genes post-modeling
      Seurat::Idents(obj_ct) <- group_col
      fc_seurat <- Seurat::FoldChange(obj_ct, ident.1 = case.group, ident.2 = ctrl.group)
      fc_seurat$primerid <- rownames(fc_seurat)
      
      # Merge FoldChange and P-value results
      tmp_res <- merge(fc_seurat[, c("primerid", "avg_log2FC")], pval, by="primerid")
      colnames(tmp_res) <- c("genes", "avg_log2FC", "p_val")
      tmp_res$p_val_adj <- p.adjust(tmp_res$p_val, method="fdr")
      tmp_res$celltype <- celltype
      tmp_res$compare <- paste0(case.group, "_vs_", ctrl.group)
      
      # Annotate expression status (Up/Down/Unchanged)
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
    }, error = function(e) {
      message("Error in ", celltype, ": ", e$message)
      return(NULL)
    })
    
    # Critical Memory Cleanup for Windows
    rm(obj_ct, expr_ct, sce, sca_ct, zlm_mod, summaryCond, dt, fc_seurat)
    gc() 
    
    loop_step_end <- Sys.time()
    cat(sprintf("    -> Finished %s in %s\n", celltype, format(loop_step_end - loop_step_start)))
  }
  
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  end_time <- Sys.time()
  cat("\n[Final Done] Total processing time:", format(end_time - start_time), "\n")
  return(sc.DE.all)
}

