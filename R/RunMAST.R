RunMAST <- function(seurat_obj, 
                    celltype_col="celltypes",
                    group_col="condition",
                    case.group="WT",
                    ctrl.group="KO",
                    preserved.genes=c("Igsf11", "IGSF11", "Vsir", "VSIR"), 
                    min_cells=3, 
                    min_expr=0.001, 
                    logfc_cutoff=0.25, 
                    cut.off.method="p_val_adj", 
                    significant.cut.off=0.05, 
                    filter_rows=FALSE, ncores=8) {
  
  start_time <- Sys.time()
  cat("[0] Initializing Windows Parallel Backend (SnowParam)...\n")
  bp_param <- BiocParallel::SnowParam(workers = ncores, progressbar = FALSE)
  BiocParallel::register(bp_param) 
  
  cat("[1] Merging Seurat layers and extracting compact dgCMatrix...\n")
  flush.console()
  seurat_obj <- Seurat::JoinLayers(seurat_obj)
  expr <- as(Seurat::GetAssayData(seurat_obj, assay="RNA", layer="data"), "dgCMatrix")
  meta <- seurat_obj@meta.data
  
  cat("[2] Applying strict noise filtering ( >", min_expr, ") and protecting targets...\n")
  flush.console()
  keep_genes_logical <- Matrix::rowSums(expr > min_expr) >= min_cells
  keep_genes_logical[rownames(expr) %in% preserved.genes] <- TRUE 
  expr <- expr[keep_genes_logical, ]
  cat("    -> Genes remaining after filtering:", nrow(expr), "\n")
  
  cat("[3] Calculating Cellular Detection Rate (CDR)...\n")
  meta$CDR <- scale(Matrix::colSums(expr > min_expr))
  
  cat("[4] Creating SingleCellAssay object (Memory Safe Bridge)...\n")
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=expr), colData=meta)
  SummarizedExperiment::rowData(sce)$primerid <- rownames(expr)
  sca <- MAST::SceToSingleCellAssay(sce)
  
  celltypes_raw <- unique(as.character(meta[[celltype_col]]))
  valid_celltypes <- c()
  for (ct in celltypes_raw) {
    meta_ct <- meta[meta[[celltype_col]] == ct, ]
    if (case.group %in% meta_ct[[group_col]] && ctrl.group %in% meta_ct[[group_col]]) {
      valid_celltypes <- c(valid_celltypes, ct)
    }
  }
  
  total_tasks <- length(valid_celltypes)
  DEG.list <- list()
  coef_name <- paste0(group_col, case.group)
  
  cat("[5] Running Windows-optimized MAST analysis on", total_tasks, "celltypes...\n")
  flush.console()
  loop_start_time <- Sys.time()
  current_task <- 0
  
  for (celltype in valid_celltypes) {
    current_task <- current_task + 1
    current_time <- Sys.time()
    cat("\n--------------------------------------------------\n")
    cat("[", format(current_time, "%H:%M:%S"), "] Processing:", celltype, "(", current_task, "/", total_tasks, ")\n")
    
    if (current_task > 1) {
      elapsed_secs <- as.numeric(difftime(current_time, loop_start_time, units="secs"))
      avg_time <- elapsed_secs / (current_task - 1)
      rem_secs <- avg_time * (total_tasks - current_task + 1)
      cat("    -> ETA: ~", floor(rem_secs/60), "mins", round(rem_secs%%60), "secs remaining\n")
    } else {
      cat("    -> ETA: Calculating based on first run...\n")
    }
    flush.console()
    
    idx_keep <- meta[[celltype_col]] == celltype & (meta[[group_col]] %in% c(case.group, ctrl.group))
    sca_ct <- sca[, idx_keep]
    SummarizedExperiment::colData(sca_ct)[[group_col]] <- factor(SummarizedExperiment::colData(sca_ct)[[group_col]], levels=c(ctrl.group, case.group))
    SummarizedExperiment::colData(sca_ct)$CDR <- scale(SummarizedExperiment::colData(sca_ct)$CDR)
    
    tryCatch({
      fmla <- as.formula(paste("~", group_col, "+ CDR"))
      
      zlm_mod <- MAST::zlm(fmla, sca_ct, method="glm", ebayes=TRUE, parallel=TRUE)
      summaryCond <- MAST::summary(zlm_mod, doLRT=coef_name)
      dt <- summaryCond$datatable
      
      # data.table syntax
      fc <- dt[dt$contrast==coef_name & dt$component=="logFC", c("primerid", "coef")]
      pval <- dt[dt$contrast==coef_name & dt$component=="H", c("primerid", "Pr(>Chisq)")]
      res <- merge(fc, pval, by="primerid")
      colnames(res) <- c("genes", "avg_log2FC", "p_val")
      res$p_val_adj <- p.adjust(res$p_val, method="fdr")
      res$celltype <- celltype
      res$compare <- paste0(case.group, "_vs_", ctrl.group)
      
      res <- dplyr::mutate(res,
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
      
      if(filter_rows) res <- dplyr::filter(res, .data[[cut.off.method]] < significant.cut.off)
      
      DEG.list[[celltype]] <- res
      gc()
    }, error=function(e){
      cat("Error in", celltype, ":", e$message, "\n")
    })
  }
  
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  end_time <- Sys.time()
  cat("\n[Done] Total processing time:", format(end_time - start_time), "\n")
  return(sc.DE.all)
}