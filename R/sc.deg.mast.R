#' @title Systematic Cohort-Level Differential Expression Analysis
#'
#' @description
#' Performs Differential Expression (DE) analysis across all cell types in a Seurat object,
#' comparing two specific conditions. Optimized for memory-intensive tasks with real-time progress logging.
#'
#' @param seurat_obj A Seurat object. (v5 supported)
#' @param celltype_col Character. Metadata column for cell type annotations. Default: "detailed.celltypes".
#' @param group_col Character. Metadata column for experimental conditions. Default: "condition".
#' @param ident.1 Character. The identity class for the first group (e.g., "sAML").
#' @param ident.2 Character. The identity class for the second group (e.g., "MDS").
#' @param test.use Character. DE test to use. Default is "MAST".
#' @param logfc.threshold Numeric. Limit testing to genes which show at least X-fold difference. Default: 0.
#' @param min.pct Numeric. Only test genes that are detected in a minimum fraction of cells. Default: 0.
#' @param gc_interval Integer. Frequency of calling garbage collection (gc) during the loop. Default: 5.
#'
#' @return A merged data frame containing DE results for all cell types.
#' @author Hyundong Yoon
#' @export
sc.deg.mast <- function(seurat_obj,
                          ident.1,
                          ident.2,
                          celltype_col = "detailed.celltypes",
                          group_col = "condition",
                          test.use = "MAST",
                          logfc.threshold = 0,
                          min.pct = 0,
                          gc_interval = 5) {
  
  # Required libraries
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tibble)
    library(MAST)
  })
  
  # 1. Preparation: Create an aggregate identifier
  cat("[1/3] Preparing aggregate metadata...\n")
  seurat_obj$tmp_aggregate <- paste(seurat_obj[[celltype_col, drop = TRUE]], 
                                    seurat_obj[[group_col, drop = TRUE]], sep = "_")
  
  all_celltypes <- unique(as.character(seurat_obj[[celltype_col, drop = TRUE]]))
  total_types <- length(all_celltypes)
  valid_groups <- unique(seurat_obj$tmp_aggregate)
  
  DEG_list <- list()
  start_time <- Sys.time()
  
  cat(sprintf("[2/3] Starting DE Analysis for %s cell types (Test: %s)\n", total_types, test.use))
  cat(sprintf("Comparison: %s vs %s\n", ident.1, ident.2))
  
  # 2. Main Loop
  for (i in seq_along(all_celltypes)) {
    celltype <- all_celltypes[i]
    group1_id <- paste0(celltype, "_", ident.1)
    group2_id <- paste0(celltype, "_", ident.2)
    
    cat(sprintf("\n[%d/%d] (%s) Processing: %s ... ", i, total_types, format(Sys.time(), "%H:%M:%S"), celltype))
    
    # Check if both comparison groups exist for this cell type
    if (group1_id %in% valid_groups && group2_id %in% valid_groups) {
      
      # Subset to minimize memory footprint during FindMarkers
      sub_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj[[celltype_col, drop = TRUE]] == celltype])
      
      tryCatch({
        deg <- FindMarkers(
          object = sub_obj,
          ident.1 = group1_id,
          ident.2 = group2_id,
          group.by = "tmp_aggregate",
          test.use = test.use,
          slot = "data",
          logfc.threshold = logfc.threshold,
          min.pct = min.pct,
          only.pos = FALSE,
          densify = TRUE
        ) %>%
          as.data.frame() %>%
          rownames_to_column(var = "genes") %>%
          mutate(
            celltype = celltype,
            comparison = paste0(ident.1, "_vs_", ident.2),
            Expression = case_when(
              avg_log2FC > 0.25 & p_val < 0.05 ~ "Up",
              avg_log2FC < -0.25 & p_val < 0.05 ~ "Down",
              TRUE ~ "Unchanged"
            ),
            color = case_when(
              Expression == "Up" ~ "#f2210a",
              Expression == "Down" ~ "#3158e8",
              TRUE ~ "#afb0b3"
            )
          )
        
        DEG_list[[celltype]] <- deg
        cat("Success!")
        
      }, error = function(e) {
        cat(sprintf("Failed: %s", e$message))
      })
      
      # Clean up sub-object to free memory
      rm(sub_obj)
      if (i %% gc_interval == 0) gc(verbose = FALSE)
      
    } else {
      cat("Skipped (One or both groups missing)")
    }
  }
  
  # 3. Consolidation
  cat("\n\n[3/3] Merging results into final dataframe...\n")
  DEG_final_df <- dplyr::bind_rows(DEG_list)
  
  end_time <- Sys.time()
  processing_duration <- end_time - start_time
  
  cat("========================================================\n")
  cat(sprintf("Analysis Completed in: %s\n", format(processing_duration)))
  cat("========================================================\n")
  
  return(DEG_final_df)
}