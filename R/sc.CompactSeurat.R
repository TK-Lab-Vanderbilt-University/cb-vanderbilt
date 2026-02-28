#' Compact Seurat Object by Removing Unnecessary Data and Coercing to Sparse Matrix
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A compacted Seurat object with normalized RNA assay.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- sc.CompactSeurat(seurat_obj)
#' }
CompactSeurat = function(seurat_obj) {
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
