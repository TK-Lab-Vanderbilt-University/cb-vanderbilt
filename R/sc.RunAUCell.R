#' @title Fast Parallel AUCell Scoring with MSigDB (sc.RunAUCell)
#'
#' @description 
#' This function seamlessly fetches molecular signatures from MSigDB (e.g., Hallmark, C2) 
#' and calculates AUCell enrichment scores for single cells. It is heavily optimized 
#' for large datasets by utilizing `BiocParallel` for multicore processing and keeping 
#' matrices in a sparse format (`dgCMatrix`).
#'
#' @details 
#' **Optimization & Parallelization:**
#' \itemize{
#'   \item **Multicore Execution:** Automatically detects the OS. Uses `MulticoreParam` 
#'   on Linux/macOS (forking) and `SnowParam` on Windows to prevent crashes while 
#'   maximizing CPU utilization.
#'   \item **Memory Efficiency:** Temporarily disables `plotStats` to prevent RAM bottlenecks.
#' }
#' 
#' **MSigDB Integration:**
#' Supports dynamic fetching of multiple MSigDB categories (e.g., `category = c("H", "C2")`).
#'
#' @param seurat_obj A Seurat object.
#' @param species Character. Species name for MSigDB (default: "Homo sapiens").
#' @param categories Character vector. MSigDB categories to fetch (default: c("H", "C2")).
#' @param assay Character. Assay to use for expression data (default: "RNA").
#' @param n_cores Integer. Number of CPU cores to use for parallel processing (default: 4).
#'
#' @return A Seurat object with calculated AUC scores appended to the metadata.
#' @export
#'
#' @import Seurat AUCell msigdbr dplyr BiocParallel
#'
#' @examples
#' # [Scoring Example]
#' # Integ_object <- sc.RunAUCell(seurat_obj = Integ_object, n_cores = 20)
sc.RunAUCell <- function(seurat_obj, 
                         species = "Homo sapiens", 
                         categories = c("H", "C2"), 
                         assay = "RNA", 
                         n_cores = 4) {
  
  # 1. Setup Parallel Backend (OS-specific to prevent Windows errors)
  if (.Platform$OS.type == "windows") {
    bp_param <- BiocParallel::SnowParam(workers = n_cores, progressbar = TRUE)
  } else {
    bp_param <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = TRUE)
  }
  
  # 2. Fetch MSigDB Pathways
  cat(sprintf("[1] Fetching MSigDB pathways for %s (Categories: %s)...\n", species, paste(categories, collapse=", ")))
  msigdb_combined <- data.frame()
  for (cat in categories) {
    tmp_df <- msigdbr::msigdbr(species = species, category = cat)
    msigdb_combined <- dplyr::bind_rows(msigdb_combined, tmp_df)
  }
  msigdb_list <- split(x = msigdb_combined$gene_symbol, f = msigdb_combined$gs_name)
  cat(sprintf(" -> Total %d Pathway lists are ready.\n", length(msigdb_list)))
  
  # 3. Extract Expression Matrix Safely
  expr_mat <- seurat_obj[[assay]]$data
  
  # 4. Build Rankings
  cat(sprintf("[2] Building gene expression rankings per cell (using %d cores)...\n", n_cores))
  cells_rankings <- AUCell::AUCell_buildRankings(
    exprMat = expr_mat, 
    plotStats = FALSE,   
    splitByBlocks = TRUE,
    BPPARAM = bp_param, 
    verbose = TRUE
  )
  
  # 5. Calculate AUC
  cat("[3] Calculating AUC for pathways...\n")
  cells_AUC <- AUCell::AUCell_calcAUC(
    geneSets = msigdb_list, 
    rankings = cells_rankings, 
    nCores = n_cores,
    verbose = TRUE
  )
  
  # 6. Add to Seurat Metadata
  cat("[4] Adding AUC scores to Seurat metadata...\n")
  auc_matrix <- t(AUCell::getAUC(cells_AUC)) 
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = as.data.frame(auc_matrix))
  
  cat("🎉 All Scoring Completed Successfully!\n")
  return(seurat_obj)
}