sc.linux.RunMAST = function(seurat_obj,
                            celltype_col="celltypes",
                            group_col="condition",
                            case.group="WT",
                            ctrl.group="KO",
                            preserved.genes=c("LAIR1", "VSIR", "Lair1", "Vsir"),
                            min_cells=3,
                            min_expr=0.001,
                            logfc_cutoff=0.25,
                            cut.off.method="p_val_adj",
                            significant.cut.off=0.05,
                            filter_rows=FALSE,
                            ncores=16) {
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
  coef_name <- paste0(group_col, case.group)
  cat("[2] Starting Sequential Analysis (Gene-level Parallelism enabled)...\n")
  for (i in seq_along(valid_celltypes)) {
    celltype <- valid_celltypes[i]
    cat(sprintf("\n[%d/%d] Processing: %s (Time: %s)\n", i, total_tasks, celltype, format(Sys.time(), "%H:%M:%S")))
    obj_ct <- subset(seurat_obj, cells = rownames(meta_all[meta_all[[celltype_col]] == celltype, ]))
    expr_ct <- Seurat::GetAssayData(obj_ct, assay="RNA", layer="data")
    keep_genes <- Matrix::rowSums(expr_ct > min_expr) >= min_cells
    keep_genes[rownames(expr_ct) %in% preserved.genes] <- TRUE
    expr_ct <- expr_ct[keep_genes, ]
    obj_ct$CDR <- scale(Matrix::colSums(expr_ct > 0))
    sce <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=as.matrix(expr_ct)),
                                                      colData=obj_ct@meta.data)
    SummarizedExperiment::rowData(sce)$primerid <- rownames(expr_ct)
    sca_ct <- MAST::SceToSingleCellAssay(sce)
    SummarizedExperiment::colData(sca_ct)[[group_col]] <- factor(
      SummarizedExperiment::colData(sca_ct)[[group_col]], levels=c(ctrl.group, case.group)
    )
    res <- tryCatch({
      fmla <- as.formula(paste("~", group_col, "+ CDR"))
      zlm_mod <- MAST::zlm(fmla, sca_ct, method="glm", ebayes=TRUE, parallel=ncores)
      summaryCond <- MAST::summary(zlm_mod, doLRT=coef_name)
      dt <- summaryCond$datatable
      fc <- dt[dt$contrast==coef_name & dt$component=="logFC", c("primerid", "coef")]
      pval <- dt[dt$contrast==coef_name & dt$component=="H", c("primerid", "Pr(>Chisq)")]
      tmp_res <- merge(fc, pval, by="primerid")
      colnames(tmp_res) <- c("genes", "avg_log2FC", "p_val")
      tmp_res$p_val_adj <- p.adjust(tmp_res$p_val, method="fdr")
      tmp_res$celltype <- celltype
      tmp_res$compare <- paste0(case.group, "_vs_", ctrl.group)
      tmp_res <- dplyr::mutate(tmp_res,
                               Expression = dplyr::case_when(
                                 avg_log2FC > logfc_cutoff & p_val_adj < significant.cut.off ~ "Up",
                                 avg_log2FC < -logfc_cutoff & p_val_adj < significant.cut.off ~ "Down",
                                 TRUE ~ "Unchanged"
                               ),
                               color = dplyr::case_when(
                                 Expression == "Up" ~ "#f2210a",
                                 Expression == "Down" ~ "#3158e8",
                                 TRUE ~ "#afb0b3"
                               )
      )
      if(filter_rows) tmp_res <- dplyr::filter(tmp_res, p_val_adj < significant.cut.off)
      DEG.list[[celltype]] <- tmp_res
    }, error=function(e){
      message("Error in ", celltype, ": ", e$message)
    })
    rm(obj_ct, expr_ct, sce, sca_ct, zlm_mod, summaryCond, dt); gc()
  }
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  cat("\n[Done] Total processing time:", format(Sys.time() - start_time), "\n")
  return(sc.DE.all)
}