#' @title High-Performance MAST Analysis for Windows (sc.windows.RunMAST)
#'
#' @description 
#' A specialized wrapper for MAST (Model-based Analysis of Single-cell Transcriptomics) 
#' optimized for Windows OS. Unlike the Linux version, this function utilizes a 
#' **Socket-based parallel backend (SnowParam)**. It processes cell types 
#' sequentially to maintain memory stability while distributing gene-level 
#' regressions across multiple CPU cores for a **significant speed boost**.
#'
#' @details 
#' **Windows-Specific Optimization:**
#' \itemize{
#'   \item **SnowParam Integration:** Explicitly manages Windows-compatible parallel 
#'   clusters, ensuring that MAST's `zlm` engine utilizes multi-threading without 
#'   the "fork" limitations of Windows.
#'   \item **Memory Safety:** By isolating one cell type at a time, it keeps the 
#'   active RAM usage low, which is critical on Windows machines where memory 
#'   fragmentation is common.
#'   \item **Efficiency:** Designed to handle complex single-cell datasets on 
#'   workstations without requiring a Linux server environment.
#' }
#' 
#' **CRITICAL WARNING:**
#' Windows handles memory differently than Linux. If `ncores` is set too high 
#' relative to your available RAM, the R session or the entire system may 
#' become unresponsive or "down." It is recommended to monitor Task Manager 
#' during your first run to find the optimal core count for your hardware.
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
#' @param ncores Integer. Number of workers for the SnowParam cluster. Default: 8.
#'
#' @return A tidy data frame with DE results, including gene symbols, log2FC, p-values, and expression status.
#' @export
#'
#' @examples
#' # [Windows Workstation Example: sAML vs. Healthy.donor]
#' # deg_results <- sc.windows.RunMAST(
#' #   seurat_obj = Integ_object,
#' #   celltype_col = "detailed.celltypes",
#' #   group_col = "condition",
#' #   case.group = "sAML",
#' #   ctrl.group = "Healthy.donor",
#' #   preserved.genes = c("LAIR1", "VSIR"),
#' #   ncores = 8
#' # )

sc.windows.RunMAST <- function(seurat_obj,
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
ncores=8) {
start_time <- Sys.time()
cat(sprintf("[0] Initializing Windows Parallel Backend (SnowParam) with %d workers...\n", ncores))
bp_param <- BiocParallel::SnowParam(workers = ncores, progressbar = FALSE)
cat("[1] Joining Seurat layers for unified analysis...\n")
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
cat("[2] Starting Sequential Cell-type Processing with Gene-level Parallelism...\n")
for (i in seq_along(valid_celltypes)) {
celltype <- valid_celltypes[i]
loop_step_start <- Sys.time()
cat(sprintf("\n[%d/%d] Current Target: %s | Start: %s\n",
i, total_tasks, celltype, format(loop_step_start, "%H:%M:%S")))
cells_to_keep <- rownames(meta_all[meta_all[[celltype_col]] == celltype, ])
obj_ct <- subset(seurat_obj, cells = cells_to_keep)
expr_ct <- Seurat::GetAssayData(obj_ct, assay="RNA", layer="data")
keep_genes <- Matrix::rowSums(expr_ct > min_expr) >= min_cells
keep_genes[rownames(expr_ct) %in% preserved.genes] <- TRUE
expr_ct <- expr_ct[keep_genes, ]
obj_ct$CDR <- scale(Matrix::colSums(expr_ct > 0))
sce <- SingleCellExperiment::SingleCellExperiment(
assays = list(logcounts = as.matrix(expr_ct)),
colData = obj_ct@meta.data
)
SummarizedExperiment::rowData(sce)$primerid <- rownames(expr_ct)
sca_ct <- MAST::SceToSingleCellAssay(sce)
SummarizedExperiment::colData(sca_ct)[[group_col]] <- factor(
SummarizedExperiment::colData(sca_ct)[[group_col]], levels = c(ctrl.group, case.group)
)
res <- tryCatch({
fmla <- as.formula(paste("~", group_col, "+ CDR"))
zlm_mod <- MAST::zlm(fmla, sca_ct, method="glm", ebayes=TRUE, parallel = bp_param)
summaryCond <- MAST::summary(zlm_mod, doLRT = coef_name)
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
}, error = function(e) {
message("Error in ", celltype, ": ", e$message)
return(NULL)
})
rm(obj_ct, expr_ct, sce, sca_ct, zlm_mod, summaryCond, dt)
gc()
loop_step_end <- Sys.time()
cat(sprintf("    -> Finished %s in %s\n", celltype, format(loop_step_end - loop_step_start)))
}
sc.DE.all <- dplyr::bind_rows(DEG.list)
end_time <- Sys.time()
cat("\n[Final Done] Total processing time:", format(end_time - start_time), "\n")
return(sc.DE.all)
}
