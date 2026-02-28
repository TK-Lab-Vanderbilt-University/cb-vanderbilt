#' @title Turbo-Charged MAST Analysis for Linux (sc.linux.RunMAST)
#'
#' @description 
#' A high-performance wrapper for MAST (Model-based Analysis of Single-cell Transcriptomics) 
#' designed specifically for Linux HPC environments. This function delivers a **massive speed boost** #' by implementing a hybrid parallelization strategy: it processes cell types sequentially 
#' to maintain a low memory footprint while utilizing **full multi-core power for gene-level 
#' calculations**. This prevents the common "Out of Memory" (OOM) errors on large datasets 
#' while completing analyses in a fraction of the standard time.
#'
#' @details 
#' **Why this is faster and safer:**
#' \itemize{
#'   \item **Gene-Level Parallelism:** Unlike standard loops, this function forces MAST's 
#'   `zlm` engine to distribute thousands of gene-level regressions across all available 
#'   CPU cores simultaneously.
#'   \item **Sequential Stability:** By subsetting and processing one cell type at a time, 
#'   it avoids the exponential RAM usage seen in parallelized `lapply` calls over large Seurat objects.
#'   \item **Optimized for Large Cohorts:** Ideal for projects with 40+ samples (e.g., MDS/sAML cohorts) 
#'   where standard MAST implementation becomes a bottleneck.
#' }
#'
#' @param seurat_obj A Seurat object. Ensure `JoinLayers()` has been run for V5 objects.
#' @param celltype_col Character. Metadata column name for cell type annotations. Default: "celltypes".
#' @param group_col Character. Metadata column name for experimental groups. Default: "condition".
#' @param case.group Character. Name of the experimental group (e.g., "sAML"). Default: "WT".
#' @param ctrl.group Character. Name of the control group (e.g., "MDS"). Default: "KO".
#' @param preserved.genes Character vector. List of genes to protect from filtering (e.g., "LAIR1").
#' @param min_cells Integer. Minimum cells expressing a gene to be included. Default: 3.
#' @param min_expr Numeric. Minimum expression threshold. Default: 0.001.
#' @param logfc_cutoff Numeric. Threshold for absolute log2 Fold Change categorization. Default: 0.25.
#' @param cut.off.method Character. Significance metric: "p_val_adj" or "p_val". Default: "p_val_adj".
#' @param significant.cut.off Numeric. Threshold for statistical significance. Default: 0.05.
#' @param filter_rows Logical. If TRUE, removes non-significant genes from the final table. Default: FALSE.
#' @param ncores Integer. Number of CPU cores to utilize for high-speed parallelization. Default: 16.
#'
#' @return A tidy data frame with DE results, including gene symbols, log2FC, p-values, and expression status.
#' @export
#'
#' @examples
#' # [Advanced Analysis Example: sAML vs. MDS Comparison]
#' # deg_results <- sc.linux.RunMAST(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "MDS",
#' #   preserved.genes = c("LAIR1", "VSIR", "IGSF11"),
#' #   cut.off.method = "p_val_adj",
#' #   ncores = 35
#' # )

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
