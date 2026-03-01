#' @title Memory-Efficient Seurat Wilcoxon Analysis (sc.RunWilcox)
#'
#' @description 
#' A highly optimized, cross-platform wrapper for Seurat's native `FindMarkers` 
#' (Wilcoxon Rank Sum test), designed for large-scale single-cell datasets. 
#' It strictly preserves standard Seurat analytical logic while dramatically 
#' improving memory efficiency for clinical cohorts. By isolating cell types 
#' sequentially and rigorously forcing garbage collection (`gc()`), it prevents 
#' the severe RAM fragmentation and Out-of-Memory (OOM) errors typically 
#' encountered when running `FindMarkers` across dozens of samples.
#'
#' @details 
#' **Performance & Integrity:**
#' \itemize{
#'   \item **Seurat Native Engine:** Executes `Seurat::FindMarkers(test.use = "wilcox")` 
#'   under the hood, ensuring 100% identical outputs to standard Seurat workflows.
#'   \item **Dynamic Aggregation:** Automatically generates aggregate metadata columns 
#'   (e.g., `celltype_condition`) to accurately subset and contrast groups.
#'   \item **Memory Safety:** Clears large `Seurat` sub-objects and temporary data frames 
#'   immediately after each cell type's iteration, ensuring stable execution on any OS.
#' }
#'
#' @param seurat_obj A Seurat object. Ensure `JoinLayers()` has been run for V5 objects.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param case.group Character. Name of the test group (e.g., "sAML"). Default: "sAML".
#' @param ctrl.group Character. Name of the reference group (e.g., "Healthy.donor"). Default: "Healthy.donor".
#' @param logfc_threshold Numeric. Limit testing to genes showing, on average, at least this fold change. Default: 0.25.
#' @param min_pct Numeric. Only test genes detected in at least this fraction of cells. Default: 0.
#' @param cut.off.method Character. Significance metric to filter by ("p_val_adj" or "p_val"). Default: "p_val".
#' @param significant.cut.off Numeric. Alpha threshold for statistical significance. Default: 0.05.
#' @param filter_rows Logical. If TRUE, the output table only retains significant DEGs based on cut-offs. Default: TRUE.
#'
#' @return A tidy data frame with Wilcoxon DE results across all specified cell types.
#' @export
#'
#' @examples
#' # [Standard Clinical Cohort Analysis]
#' # wilcox_results <- sc.RunWilcox(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "Healthy.donor",
#' #   filter_rows = TRUE
#' # )

sc.RunWilcox = function(seurat_obj,
                        celltype_col="detailed.celltypes",
                        group_col="condition",
                        case.group="sAML",
                        ctrl.group="Healthy.donor",
                        logfc_threshold=0.25,
                        min_pct=0,
                        cut.off.method="p_val",
                        significant.cut.off=0.05,
                        filter_rows=TRUE) {
  
  start_time <- Sys.time()
  cat("[1] Joining layers and preparing meta...\n")
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj)
  
  # 2. Re-arrangement meta data (Exactly as provided)
  aggregate_col <- paste0(celltype_col, ".aggregate")
  seurat_obj@meta.data[[aggregate_col]] <- paste(seurat_obj@meta.data[[celltype_col]], 
                                                 seurat_obj@meta.data[[group_col]], sep = "_")
  
  valid_celltypes <- unique(as.character(seurat_obj@meta.data[[celltype_col]]))
  groups <- unique(seurat_obj@meta.data[[aggregate_col]])
  
  total_tasks <- length(valid_celltypes)
  DEG.list <- list()
  
  cat("[3] Starting Sequential Analysis (Memory-Optimized Seurat Wilcoxon)...\n")
  for (i in seq_along(valid_celltypes)) {
    celltype <- valid_celltypes[i]
    group1 <- paste0(celltype, "_", case.group)
    group2 <- paste0(celltype, "_", ctrl.group)
    
    # Check if both groups exist in the dataset for this cell type
    if (group1 %in% groups && group2 %in% groups) {
      loop_step_start <- Sys.time()
      cat(sprintf("\n[%d/%d] Processing Wilcoxon for: %s (Time: %s)\n", i, total_tasks, celltype, format(loop_step_start, "%H:%M:%S")))
      
      # Subset object
      sub_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[celltype_col]] == celltype, ]))
      
      tryCatch({
        # Run exact user-provided FindMarkers logic
        deg <- Seurat::FindMarkers(
          object = sub_obj,
          ident.1 = group1,
          ident.2 = group2,
          group.by = aggregate_col,
          test.use = "wilcox",
          slot = "data",  
          logfc.threshold = logfc_threshold,
          min.pct = min_pct,
          min.cells.feature = 0,
          min.cells.group = 0,
          only.pos = FALSE,
          densify = TRUE
        ) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "genes") %>%
          dplyr::mutate(celltype = celltype,
                        compare = paste0(case.group, "_vs_", ctrl.group)) %>%
          dplyr::mutate(Expression = dplyr::case_when(
            avg_log2FC > logfc_threshold ~ "Up",
            avg_log2FC < -logfc_threshold ~ "Down",
            TRUE ~ "Unchanged"
          ),
          color = dplyr::case_when(
            Expression == "Up" ~ "#f2210a",
            Expression == "Down" ~ "#3158e8",
            TRUE ~ "#afb0b3"
          ))
        
        # Apply user's filtering logic dynamically
        if(filter_rows) {
          deg <- dplyr::filter(deg, .data[[cut.off.method]] < significant.cut.off)
        }
        
        DEG.list[[celltype]] <- deg
        
      }, error = function(e){
        cat("  -> Error in", celltype, ":", e$message, "\n")
      })
      
      # CRITICAL: Force garbage collection to prevent memory explosion
      rm(sub_obj, deg); gc()
      
      loop_step_end <- Sys.time()
      cat(sprintf("    -> Finished %s in %s\n", celltype, format(loop_step_end - loop_step_start)))
      
    } else {
      cat(sprintf("\n[%d/%d] Skipping: %s - Missing groups for comparison.\n", i, total_tasks, celltype))
    }
  }
  
  sc.DE.all <- dplyr::bind_rows(DEG.list)
  cat("\n[Final Done] Total processing time:", format(Sys.time() - start_time), "\n")
  return(sc.DE.all)
}