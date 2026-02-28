#' @title Memory-Efficient Seurat Object Optimization (sc.CompactSeurat)
#'
#' @description 
#' This function optimizes and "downsizes" Seurat objects to minimize memory (RAM) 
#' usage and storage footprint. It is particularly useful for large-scale datasets 
#' (e.g., cohorts with 40+ samples) by stripping redundant data, merging V5 layers, 
#' and ensuring all matrices are stored in a highly compressed sparse format (dgCMatrix).
#'
#' @details 
#' **Optimization Steps:**
#' \itemize{
#'   \item **Layer Integration:** Automatically detects and merges Seurat V5 layers 
#'   using `JoinLayers()` to ensure a unified data structure.
#'   \item **Dietary Reduction:** Uses `DietSeurat()` to remove unnecessary graphs, 
#'   neighbor objects, and extra assays, keeping only the essential "RNA" assay.
#'   \item **Dimensionality Reduction Preservation:** Retains critical reductions: 
#'   PCA, Harmony, and UMAP.
#'   \item **Sparse Matrix Conversion:** Forces "counts" and "data" slots into 
#'   `dgCMatrix` format, which can reduce memory usage by up to 60-80% compared 
#'   to standard dense matrices.
#' }
#' 
#' **COMPATIBILITY NOTE:**
#' This function is specifically designed for **Seurat V5**. It may **not be 
#' compatible with Seurat V4** or earlier versions due to significant changes 
#' in how assays and layers are handled (Assay5 vs. Assay).
#'
#' @param seurat_obj A Seurat object (V5 recommended).
#'
#' @return A "compacted" Seurat object with optimized memory footprint and 
#' re-normalized RNA data.
#' @export
#'
#' @examples
#' # [Optimization Example]
#' # optimized_obj <- sc.CompactSeurat(seurat_obj = large_integrated_obj)
#' # 
#' # # Verify memory reduction
#' # object.size(large_integrated_obj)
#' # object.size(optimized_obj)

sc.CompactSeurat = function(seurat_obj) {
DefaultAssay(seurat_obj) <- "RNA"
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
seurat_obj = SeuratObject::JoinLayers(seurat_obj)
}
seurat_obj = DietSeurat(
object = seurat_obj,
assays = "RNA",
layers = c("counts", "data"),
dimreducs = c("pca", "harmony", "umap")
)
if (is(seurat_obj[["RNA"]]$counts, "matrix")) {
seurat_obj[["RNA"]]$counts <- as(seurat_obj[["RNA"]]$counts, "dgCMatrix")
}
if (is(seurat_obj[["RNA"]]$data, "matrix")) {
seurat_obj[["RNA"]]$data <- as(seurat_obj[["RNA"]]$data, "dgCMatrix")
}
seurat_obj = NormalizeData(seurat_obj)
return(seurat_obj)
}
