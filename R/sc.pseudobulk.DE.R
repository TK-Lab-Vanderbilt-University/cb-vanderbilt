#' Pseudobulk Differential Expression Analysis for Single-Cell Data
#'
#' @description
#' This function performs a cluster-wise pseudobulk DE analysis using a Seurat object. 
#' It aggregates counts per sample and cell type, then utilizes the edgeR 
#' Quasi-Likelihood framework combined with limma contrasts for robust statistical testing.
#'
#' @param obj A Seurat object containing the single-cell data.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "detailed.celltypes".
#' @param condition_col Character. Metadata column name for experimental conditions. Default: "condition".
#' @param sample_col Character. Metadata column name for sample/source identifiers. Default: "orig.ident".
#' @param case Character. The condition level representing the test group (e.g., "MDS").
#' @param control Character. The condition level representing the reference group (e.g., "Healthy.donor").
#'
#' @return A consolidated data frame of significant DEGs (p < 0.05) across all clusters.
#' @export
#'
#' @import Seurat
#' @import SingleCellExperiment
#' @import edgeR
#' @import limma
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import Matrix.utils
#' @import magrittr
sc.pseudobulk.DE <- function(obj, 
                             celltype_col = "detailed.celltypes", 
                             condition_col = "condition", 
                             sample_col = "orig.ident",
                             case = "MDS",
                             control = "Healthy.donor") {

# 1. Metadata Preparation
# ----------------------------------------------------------------------------
message("Step 1: Preparing metadata and identifiers...")

# Create combined ID for pseudobulk grouping
obj$pseudobulk_id_aggregate <- paste(
  obj[[celltype_col, drop = TRUE]], 
  obj[[condition_col, drop = TRUE]], 
  sep = "_"
)

# Set identity to handle Seurat v5 vector requirements
obj <- Seurat::SetIdent(obj, value = obj[["pseudobulk_id_aggregate", drop = TRUE]])

# Extract core data components
counts <- Seurat::GetAssayData(object = obj, slot = "counts", assay = "RNA")
metadata <- obj@meta.data
metadata$cluster_id <- factor(Seurat::Idents(obj))

# 2. Pseudobulk Aggregation
# ----------------------------------------------------------------------------
message("Step 2: Aggregating single-cell counts to pseudobulk profiles...")

# Convert to SingleCellExperiment object
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts), 
  colData = metadata
)

# Prepare experimental info (ei) mapping
n_cells <- table(sce[[sample_col]]) %>% as.vector()
names(n_cells) <- names(table(sce[[sample_col]]))
m <- match(names(n_cells), sce[[sample_col]])

ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>%
  dplyr::select(dplyr::all_of(c(sample_col, condition_col, "n_cells")))

# Set up grouping variables for matrix summation
groups <- colData(sce)[, c("cluster_id", sample_col)]
groups[[sample_col]] <- factor(groups[[sample_col]])

# Aggregate counts using matrix multiplication for efficiency
pb <- Matrix.utils::aggregate.Matrix(t(SingleCellExperiment::counts(sce)), 
                                     groupings = groups, fun = "sum")

# Parse rownames to extract cluster identities
splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_", n = 2), `[`, 1)

# Split aggregated matrix into a cluster-wise list
pb_list <- split.data.frame(pb, factor(splitf)) %>%
  lapply(function(u) {
    new_colnames <- gsub(".*_", "", rownames(u))
    t(u) %>% magrittr::set_colnames(new_colnames)
  })

# 3. Cluster-wise DE Analysis Loop (edgeR + limma)
# ----------------------------------------------------------------------------
message("Step 3: Executing differential expression analysis...")

all_clusters <- names(pb_list)
DE_list <- list()

for (cluster in all_clusters) {
  cat("  Processing:", cluster, "... ")
  
  y_raw <- pb_list[[cluster]]
  sample_ids <- colnames(y_raw)
  
  # Filter metadata for samples active in current cluster
  sub_ei <- ei %>% dplyr::filter(!!rlang::sym(sample_col) %in% sample_ids)
  sub_ei <- sub_ei[match(sample_ids, sub_ei[[sample_col]]), , drop = FALSE]
  
  # Check if both Case and Control groups exist in the cluster
  if (length(unique(sub_ei[[condition_col]])) < 2) {
    cat("Skipped: Missing comparison group.\n")
    next
  }
  
  tryCatch({
    # Initialize edgeR DGEList
    y <- edgeR::DGEList(y_raw, remove.zeros = TRUE)
    y <- edgeR::calcNormFactors(y)
    
    # Construct Design Matrix
    sub_ei[[condition_col]] <- factor(sub_ei[[condition_col]])
    design_formula <- as.formula(paste("~ 0 +", condition_col))
    design <- model.matrix(design_formula, data = sub_ei)
    colnames(design) <- levels(sub_ei[[condition_col]])
    
    # Contrast Definition (limma required)
    contrast_logic <- paste0(case, " - ", control)
    contrast <- limma::makeContrasts(contrasts = contrast_logic, levels = design)
    
    # Model Fitting
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
    
    # Results Processing
    res_table <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table %>%
      dplyr::mutate(gene = rownames(.), cluster_id = cluster) %>%
      dplyr::rename(p_val = PValue, p_adj = FDR)
    
    DE_list[[cluster]] <- res_table
    cat("Success!\n")
    
  }, error = function(e) {
    cat("Error in cluster", cluster, ":", e$message, "\n")
  })
}

# 4. Final Formatting and Consolidation
# ----------------------------------------------------------------------------
message("Step 4: Finalizing and merging results...")

if (length(DE_list) == 0) {
  warning("No significant DEGs found or no clusters were successfully processed.")
  return(data.frame())
}

# Apply significance-based coloring and regulation labels
DE_final_list <- lapply(DE_list, function(DE) {
  DE %>%
    dplyr::mutate(
      color = dplyr::case_when(
        logFC > 0.25 & p_val < 0.05 ~ '#f2210a',
        logFC < -0.25 & p_val < 0.05 ~ '#3158e8',
        TRUE ~ "#afb0b3"
      ),
      Regulation = dplyr::case_when(
        logFC > 0.5 ~ "Up",
        logFC < -0.5 ~ "Down",
        TRUE ~ "Unchanged"
      )
    ) %>%
    # Retain only significant entries (p < 0.05)
    dplyr::filter(p_val < 0.05 & !is.na(p_val))
})

# Merge list into a single long-format data frame
DE.all.df <- do.call(rbind, DE_final_list)

message("[COMPLETE] All clusters processed successfully.")
return(DE.all.df)
}
