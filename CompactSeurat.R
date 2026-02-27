CompactSeurat <- function(seurat_obj) {
  # Load library within function to be safe
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat package is required.")
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  
  # Join layers if it's Seurat v5
  if (inherits(seurat_obj[["RNA"]], "Assay5")) {
    seurat_obj <- Seurat::JoinLayers(seurat_obj)
  }
  
  # DietSeurat to keep only necessary data
  seurat_obj <- Seurat::DietSeurat(
    object = seurat_obj,
    assays = "RNA",
    layers = c("counts", "data"),
    dimreducs = c("pca", "harmony", "umap")
  )
  
  # Final normalization
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  
  return(seurat_obj)
}
