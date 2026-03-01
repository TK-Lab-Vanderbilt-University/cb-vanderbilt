#' @title Robust Pseudo-bulk Differential Expression Analysis via edgeR
#'
#' @description 
#' Performs pseudo-bulk DE analysis by aggregating raw counts across biological replicates. 
#' This function implements the edgeR Likelihood Ratio Test (LRT) framework, which provides 
#' a balanced and scientifically defensible sensitivity (less conservative than QLFit) 
#' for clinical single-cell datasets.
#'
#' @details 
#' **Key Improvements:**
#' \itemize{
#'   \item **Integer Summation:** Replaces "mean" with "sum" for count aggregation, 
#'   adhering to edgeR's Negative Binomial requirements and boosting statistical power.
#'   \item **FDR Rescue:** Uses `edgeR::filterByExpr` to remove low-count noise genes, 
#'   drastically improving the adjusted p-value (FDR) for true biological targets.
#'   \item **Memory Efficiency:** Avoids heavy Intermediate objects (like SCE) and 
#'   utilizes sequential processing with explicit garbage collection (`gc()`) to 
#'   run safely on machines with limited RAM.
#'   \item **Scientific Defensibility:** Uses the Likelihood Ratio Test (LRT), a 
#'   gold-standard in bulk RNA-seq, which is less prone to false negatives in 
#'   clinical cohorts compared to the overly strict QLF test.
#' }
#'
#' @param seurat_obj A Seurat object containing raw 'counts'.
#' @param celltype_col Character. Metadata column for cell type annotations. Default: "detailed.celltypes".
#' @param sample_col Character. Metadata column for biological replicates (e.g., "orig.ident"). Default: "orig.ident".
#' @param group_col Character. Metadata column for experimental groups (e.g., "condition"). Default: "condition".
#' @param case.group Character. Name of the test group (e.g., "MDS"). Default: "MDS".
#' @param ctrl.group Character. Name of the control group (e.g., "Healthy.donor"). Default: "Healthy.donor".
#' @param logfc_cutoff Numeric. Threshold for log2 fold change categorization. Default: 0.25.
#' @param significant.cut.off Numeric. P-value threshold for significance labeling. Default: 0.05.
#' @param filter_rows Logical. If TRUE, returns only significant results. Default: FALSE.
#'
#' @return A tidy data frame with pseudo-bulk DE results.
#' @export

sc.pseudobulk.DE = function(seurat_obj,
                            celltype_col="detailed.celltypes",
                            sample_col="orig.ident",
                            group_col="condition",
                            case.group="MDS",
                            ctrl.group="Healthy.donor",
                            logfc_cutoff=0.25,
                            significant.cut.off=0.05,
                            filter_rows=FALSE) {
  
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install 'edgeR'.")
  if (!requireNamespace("Matrix.utils", quietly = TRUE)) stop("Please install 'Matrix.utils'.")
  
  start_time <- Sys.time()
  
  # [1] Data Preparation
  cat("[1] Extracting raw counts and building metadata...\n")
  counts <- Seurat::GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
  meta <- seurat_obj@meta.data
  
  # Create a robust grouping string using a unique delimiter
  # This prevents errors if celltype names contain underscores
  group_strings <- paste(meta[[celltype_col]], meta[[sample_col]], sep = "|||")
  
  # [2] Aggregation (CRITICAL: sum for integer counts)
  cat("[2] Summing counts into pseudo-bulk profiles...\n")
  pb_mat <- Matrix.utils::aggregate.Matrix(t(counts), groupings = group_strings, fun = "sum")
  pb_mat <- t(pb_mat) # Genes x (Celltype|||Sample)
  
  # Parse back the aggregated names
  pb_colnames <- colnames(pb_mat)
  parsed_ct <- sapply(strsplit(pb_colnames, "|||", fixed = TRUE), `[`, 1)
  parsed_sample <- sapply(strsplit(pb_colnames, "|||", fixed = TRUE), `[`, 2)
  
  # Consolidate sample-level metadata
  ei <- unique(meta[, c(sample_col, group_col)])
  colnames(ei) <- c("sample_id", "condition")
  
  all_clusters <- unique(parsed_ct)
  DE_list <- list()
  
  cat("[3] Starting Sequential edgeR LRT Analysis...\n")
  for (cluster in all_clusters) {
    cat(sprintf(" -> Analyzing Cluster: %s\n", cluster))
    
    # Subset matrix for the current cell type
    idx <- which(parsed_ct == cluster)
    y_mat <- pb_mat[, idx, drop = FALSE]
    colnames(y_mat) <- parsed_sample[idx]
    
    # Align metadata
    sub_ei <- ei[ei$sample_id %in% colnames(y_mat), , drop = FALSE]
    sub_ei <- sub_ei[match(colnames(y_mat), sub_ei$sample_id), ]
    
    # Define factor levels (Control vs Case)
    sub_ei$condition <- factor(sub_ei$condition, levels = c(ctrl.group, case.group))
    
    # Check for sufficient replicates (at least 2 per group for dispersion estimation)
    if (length(levels(sub_ei$condition)) < 2 || any(table(sub_ei$condition) < 2)) {
      cat("    Skip: Insufficient biological replicates.\n")
      next
    }
    
    tryCatch({
      # Initialize edgeR object
      y <- edgeR::DGEList(counts = y_mat, samples = sub_ei)
      
      # [CRITICAL] filterByExpr: Removes low-count genes to boost FDR power
      design <- model.matrix(~ condition, data = y$samples)
      keep <- edgeR::filterByExpr(y, design = design)
      y <- y[keep, , keep.lib.sizes = FALSE]
      
      # Normalization and Dispersion Estimation
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, design)
      
      # [CRITICAL] glmLRT: More sensitive (less strict) than QLFit for clinical data
      fit <- edgeR::glmFit(y, design)
      lrt <- edgeR::glmLRT(fit, coef = 2) 
      
      # Formatting Results
      res <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
      res <- tibble::rownames_to_column(res, var = "gene") %>%
        dplyr::mutate(cluster_id = cluster,
                      compare = paste0(case.group, "_vs_", ctrl.group)) %>%
        dplyr::rename(p_val = PValue, p_adj = FDR, log2FC = logFC)
      
      DE_list[[cluster]] <- res
    }, error = function(e) {
      cat("    Error in", cluster, ":", e$message, "\n")
    })
    
    # Explicit memory cleanup
    rm(y, fit, lrt); gc()
  }
  
  # [4] Post-processing and Labeling
  cat("[4] Merging results and assigning regulation status...\n")
  final_res <- dplyr::bind_rows(DE_list) %>%
    dplyr::mutate(
      Regulation = dplyr::case_when(
        log2FC > logfc_cutoff & p_val < significant.cut.off ~ "Up",
        log2FC < -logfc_cutoff & p_val < significant.cut.off ~ "Down",
        TRUE ~ "Unchanged"
      ),
      color = dplyr::case_when(
        Regulation == "Up" ~ "#f2210a",
        Regulation == "Down" ~ "#3158e8",
        TRUE ~ "#afb0b3"
      )
    ) %>%
    dplyr::filter(!is.na(p_val))
  
  if(filter_rows) {
    final_res <- final_res %>% dplyr::filter(p_val < significant.cut.off)
  }
  
  cat("\n[Final Done] Total processing time:", format(Sys.time() - start_time), "\n")
  return(final_res)
}