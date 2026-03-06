#' Execute the Full MultiNicheNet Analysis Pipeline with Enhanced Parameters
#'
#' @description
#' A comprehensive wrapper function that executes the entire MultiNicheNetR workflow. 
#' It supports specific group selection, contrast definitions, and rigorous filtering thresholds.
#'
#' @param sce A preprocessed SingleCellExperiment object.
#' @param lr_db_list A list containing the `lr_network` and `ligand_target_matrix`.
#' @param group_oi Character string indicating the group of interest (e.g., "MDS").
#' @param group_id Metadata column for condition/grouping (e.g., "condition").
#' @param celltype_id Metadata column for cell type annotations.
#' @param conditions_keep Character vector of conditions to retain for analysis.
#' @param contrasts_oi Character string defining the DE contrasts.
#' @param contrast_tbl Tibble defining the mapping between contrasts and groups.
#' @param sample_id Metadata column for sample identifiers.
#' @param covariates Optional metadata columns for DE model adjustment. Default is NA.
#' @param batches Optional metadata columns for batch effect adjustment. Default is NA.
#' @param min_cells Minimum number of cells per cell type. Default is 0.
#' @param min_sample_prop Minimum proportion of samples that must meet min_cells. Default is 0.10.
#' @param fraction_cutoff Minimum fraction of cells expressing a gene. Default is 0.01.
#' @param logFC_threshold Log-fold change threshold for DE genes. Default is 0.25.
#' @param p_val_threshold P-value threshold for DE genes. Default is 0.05.
#' @param p_val_adj Logical; whether to use adjusted p-values. Default is FALSE.
#' @param empirical_pval Logical; whether to use empirical p-values. Default is FALSE.
#' @param top_n_target Number of top targets to consider for ligand activity. Default is 250.
#' @param ligand_activity_down Logical; whether to consider downregulated ligand activity. Default is FALSE.
#' @param n.cores Number of CPU cores for parallel processing. Default is 30.
#' @param verbose Logical; whether to print progress messages. Default is TRUE.
#'
#' @return A list containing prioritized interaction tables and metadata.
#' 
#' @export
#' @import MultiNicheNet
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble tibble
#' @importFrom dplyr mutate group_by summarise filter count pull setdiff intersect distinct select rename
run_multinichenet_pipeline <- function(sce, 
                                       lr_db_list, 
                                       group_oi = "MDS",
                                       group_id = "condition", 
                                       celltype_id = "detailed.celltypes", 
                                       conditions_keep = c("MDS", "Healthy.donor"),
                                       contrasts_oi = c("'MDS-Healthy.donor','Healthy.donor-MDS'"),
                                       contrast_tbl = tibble::tibble(contrast = c("MDS-Healthy.donor", "Healthy.donor-MDS"), group = c("MDS", "Healthy.donor")),
                                       sample_id = "orig.ident",
                                       covariates = NA, 
                                       batches = NA,
                                       min_cells = 0, 
                                       min_sample_prop = 0.10, 
                                       fraction_cutoff = 0.01,
                                       logFC_threshold = 0.25, 
                                       p_val_threshold = 0.05, 
                                       p_val_adj = FALSE, 
                                       empirical_pval = FALSE,
                                       top_n_target = 250, 
                                       ligand_activity_down = FALSE, 
                                       n.cores = 30,
                                       verbose = TRUE) {
  
  lr_network <- lr_db_list$lr_network
  ligand_target_matrix <- lr_db_list$ligand_target_matrix
  
  # 0. Subset SCE based on conditions_keep
  if(verbose) message("Subsetting SCE for conditions: ", paste(conditions_keep, collapse = ", "))
  sce <- sce[, SummarizedExperiment::colData(sce)[,group_id] %in% conditions_keep]
  
  # 1. Define Senders and Receivers (Auto-detection from the subset)
  senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() %>% as.character()
  receivers_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() %>% as.character()
  
  # 2. Get Abundance & Cell type Filtering
  abundance_info <- get_abundance_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches
  )
  
  # Filter out cell types present in fewer than 2 samples per condition (MultiNicheNet requirement)
  abundance_df_summarized <- abundance_info$abundance_data %>% 
    dplyr::mutate(keep = as.logical(keep)) %>% 
    dplyr::group_by(group_id, celltype_id) %>% 
    dplyr::summarise(samples_present = sum((keep)), .groups = "drop")
  
  total_nr_conditions <- length(unique(SummarizedExperiment::colData(sce)[,group_id]))
  absent_celltypes <- abundance_df_summarized %>% 
    dplyr::filter(samples_present < 2) %>% 
    dplyr::group_by(celltype_id) %>% 
    dplyr::count() %>% 
    dplyr::filter(n == total_nr_conditions) %>% 
    dplyr::pull(celltype_id)
  
  senders_oi <- setdiff(senders_oi, absent_celltypes)
  receivers_oi <- setdiff(receivers_oi, absent_celltypes)
  sce <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  
  # 3. Expression Fraction & Gene Filtering
  frq_list <- get_frac_exprs(
    sce = sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id, 
    batches = batches, min_cells = min_cells, fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop
  )
  
  genes_oi <- frq_list$expressed_df %>% dplyr::filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
  sce <- sce[genes_oi, ]
  
  # 4. Processing Abundance & Expression Info
  abundance_expression_info <- process_abundance_expression_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, 
    lr_network = lr_network, batches = batches, frq_list = frq_list, abundance_info = abundance_info
  )
  
  # 5. Differential Expression (DE) Analysis
  if(verbose) message("Performing Differential Expression Analysis...")
  DE_info <- get_DE_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, 
    min_cells = min_cells, expressed_df = frq_list$expressed_df
  )
  
  if(empirical_pval) {
    DE_info_emp <- get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
    celltype_de <- DE_info_emp$de_output_tidy_emp %>% 
      dplyr::select(-p_val, -p_adj) %>% 
      dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
  } else {
    celltype_de <- DE_info$celltype_de$de_output_tidy
  } 
  
  sender_receiver_de <- combine_sender_receiver_de(
    sender_de = celltype_de, receiver_de = celltype_de, senders_oi = senders_oi, 
    receivers_oi = receivers_oi, lr_network = lr_network
  )
  
  # 6. Ligand Activity Prediction
  if(verbose) message("Predicting Ligand Activities...")
  ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
    get_ligand_activities_targets_DEgenes(
      receiver_de = celltype_de, 
      receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
      ligand_target_matrix = ligand_target_matrix, 
      logFC_threshold = logFC_threshold,
      p_val_threshold = p_val_threshold, 
      p_val_adj = p_val_adj, 
      top_n_target = top_n_target,
      verbose = verbose, 
      n.cores = n.cores
    )
  ))
  
  # 7. Metadata Grouping & Prioritization
  sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
  metadata_combined <- SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  
  grouping_tbl <- metadata_combined[, c(sample_id, group_id, if(!is.na(batches)) batches else NULL)] %>% 
    dplyr::distinct()
  colnames(grouping_tbl)[1:2] <- c("sample", "group")
  
  if(verbose) message("Generating Prioritization Tables...")
  prioritization_tables <- suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de, 
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl, 
    sender_receiver_tbl = sender_receiver_tbl, 
    grouping_tbl = grouping_tbl,
    scenario = "regular", 
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender, 
    ligand_activity_down = ligand_activity_down
  ))
  
  # 8. Inference of Ligand-Target Correlation
  lr_target_prior_cor <- lr_target_prior_cor_inference(
    receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
    abundance_expression_info = abundance_expression_info, 
    celltype_de = celltype_de, 
    grouping_tbl = grouping_tbl, 
    prioritization_tables = prioritization_tables, 
    ligand_target_matrix = ligand_target_matrix, 
    logFC_threshold = logFC_threshold, 
    p_val_threshold = p_val_threshold, 
    p_val_adj = p_val_adj
  )
  
  # Consolidate Output
  multinichenet_output <- list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
  
  if(verbose) message("Analysis Complete.")
  return(make_lite_output(multinichenet_output))
}
