#' Pseudobulk Differential Expression Analysis for Seurat Objects
#'
#' @description
#' This function takes a Seurat object, aggregates counts into pseudobulk profiles 
#' based on samples and clusters, and performs DE analysis using edgeR. 
#' It dynamically filters samples for each cluster to ensure statistical validity.
#'
#' @param obj A Seurat object.
#' @param celltype_col Character. Metadata column name for cell types. Default: "detailed.celltypes".
#' @param condition_col Character. Metadata column name for conditions. Default: "condition".
#' @param sample_col Character. Metadata column name for sample IDs. Default: "orig.ident".
#' @param case Character. The test group condition (e.g., "MDS").
#' @param control Character. The reference/control group condition (e.g., "Healthy.donor").
#'
#' @return A consolidated data frame of significant DEGs across all clusters.
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

  # ----------------------------------------------------------------------------
  # 1. Metadata & Pseudobulk Preparation
  # ----------------------------------------------------------------------------
  message("Step 1: Preparing metadata and aggregating counts...")
  
  # Create aggregate ID for grouping
  obj$pseudobulk_id_aggregate <- paste(
    obj[[celltype_col, drop = TRUE]], 
    obj[[condition_col, drop = TRUE]], 
    sep = "_"
  )

  # Set identity and extract counts
  obj <- Seurat::SetIdent(obj, value = obj[["pseudobulk_id_aggregate", drop = TRUE]])
  counts <- Seurat::GetAssayData(object = obj, slot = "counts", assay = "RNA")
  metadata <- obj@meta.data
  metadata$cluster_id <- factor(Seurat::Idents(obj))

  # Convert to SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts), 
    colData = metadata
  )

  # Prepare Experimental Info (ei)
  n_cells <- table(sce[[sample_col]]) %>% as.vector()
  names(n_cells) <- names(table(sce[[sample_col]]))
  m <- match(names(n_cells), sce[[sample_col]])

  ei <- data.frame(colData(sce)[m, ], n_cells, row.names = NULL) %>%
    dplyr::select(dplyr::all_of(c(sample_col, condition_col, "n_cells")))

  # Aggregate counts
  groups <- colData(sce)[, c("cluster_id", sample_col)]
  groups[[sample_col]] <- factor(groups[[sample_col]])
  
  pb_matrix <- Matrix.utils::aggregate.Matrix(t(SingleCellExperiment::counts(sce)), 
                                              groupings = groups, fun = "sum")
  
  # Split aggregated matrix into a cluster-wise list (pb)
  splitf <- sapply(stringr::str_split(rownames(pb_matrix), pattern = "_", n = 2), `[`, 1)
  pb <- split.data.frame(pb_matrix, factor(splitf)) %>%
    lapply(function(u) {
      new_colnames <- gsub(".*_", "", rownames(u))
      t(u) %>% magrittr::set_colnames(new_colnames)
    })

  # ----------------------------------------------------------------------------
  # 2. Cluster-wise DE Analysis Loop (User Provided Logic)
  # ----------------------------------------------------------------------------
  message("Step 2: Running differential expression analysis per cluster...")
  
  all_clusters <- names(pb)
  DE_list <- list()

  for (cluster in all_clusters) {
    cat("Processing:", cluster, "\n")
    
    # 1. Extract count matrix
    y_raw <- pb[[cluster]]
    
    # 2. Identify samples in this cluster
    sample_ids <- colnames(y_raw)
    
    # 3. Subset and align metadata
    sub_ei <- ei %>% dplyr::filter(!!rlang::sym(sample_col) %in% sample_ids)
    sub_ei <- sub_ei[match(sample_ids, sub_ei[[sample_col]]), , drop = FALSE]
    
    # Skip if conditions are missing
    if (length(unique(sub_ei[[condition_col]])) < 2) {
      cat("  - Skipped: Not enough conditions present in this cluster.\n")
      next
    }
    
    tryCatch({
      # edgeR Analysis
      y <- edgeR::DGEList(y_raw, remove.zeros = TRUE)
      y <- edgeR::calcNormFactors(y)
      
      # Design matrix
      design <- model.matrix(as.formula(paste("~ 0 + factor(", condition_col, ")")), data = sub_ei)
      colnames(design) <- levels(factor(sub_ei[[condition_col]]))
      rownames(design) <- sub_ei[[sample_col]]
      
      # Validate contrast groups
      if (!all(c(case, control) %in% colnames(design))) {
        cat("  - Skipped: Missing conditions (", case, " or ", control, ").\n", sep="")
        next
      }
      
      # Contrast & Fit
      contrast_logic <- paste0(case, " - ", control)
      contrast <- limma::makeContrasts(contrasts = contrast_logic, levels = design)
      
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmQLFit(y, design)
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
      
      # Results formatting
      res_table <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table %>%
        dplyr::mutate(
          gene = rownames(.), 
          cluster_id = cluster,
          # Regulation logic based on 0.25/-0.25 threshold
          Regulation = dplyr::case_when(
            logFC > 0.25 ~ "Up",
            logFC < -0.25 ~ "Down",
            TRUE ~ "Unchanged"
          )
        ) %>%
        dplyr::rename(p_val = PValue, p_adj = FDR)
      
      DE_list[[cluster]] <- res_table
      cat("  - Success!\n")
      
    }, error = function(e) {
      cat("  - Error in", cluster, ":", e$message, "\n")
    })
  }

  # ----------------------------------------------------------------------------
  # 3. Consolidate Results
  # ----------------------------------------------------------------------------
  if (length(DE_list) == 0) {
    message("No DEGs found across clusters.")
    return(data.frame())
  }
  
  final_pb_DE <- dplyr::bind_rows(DE_list)
  message("\n[DONE] Pseudobulk DE analysis complete.")
  
  return(final_pb_DE)
}
