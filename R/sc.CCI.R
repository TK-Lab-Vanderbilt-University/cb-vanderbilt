#' Build and Extend Ligand-Receptor Database for MultiNicheNet
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate distinct bind_rows filter
#' @importFrom MultiNicheNet convert_alias_to_symbols convert_human_to_mouse_symbols
build_custom_lr_db <- function(organism = "human", 
                               custom_ligands = c("SELPLG", "LGALS9", "MMP13", "LRIG1", "SDC2", "VSIG8", "IGSF11"), 
                               custom_receptor = "VSIR") {
  
  options(timeout = 120)
  
  if (organism == "human") {
    lr_network_all <- readRDS(url("https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds")) %>% 
      mutate(ligand = convert_alias_to_symbols(ligand, organism = organism), 
             receptor = convert_alias_to_symbols(receptor, organism = organism))
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    
    cur_ligands <- make.names(custom_ligands)
    cur_receptor <- make.names(custom_receptor)
  } else if (organism == "mouse") {
    lr_network_all <- readRDS(url("https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds")) %>% 
      mutate(ligand = convert_alias_to_symbols(ligand, organism = organism), 
             receptor = convert_alias_to_symbols(receptor, organism = organism))
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    
    cur_ligands <- convert_human_to_mouse_symbols(custom_ligands) %>% make.names()
    cur_receptor <- convert_human_to_mouse_symbols(custom_receptor) %>% make.names()
  }
  
  lr_network_all <- lr_network_all %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  lr_network <- lr_network_all %>% distinct(ligand, receptor)
  
  custom_lr_pairs <- data.frame(ligand = cur_ligands, receptor = cur_receptor)
  lr_network <- bind_rows(lr_network, custom_lr_pairs) %>% distinct()
  
  colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  missing_ligands <- setdiff(cur_ligands, colnames(ligand_target_matrix))
  if (length(missing_ligands) > 0) {
    for (l in missing_ligands) {
      zero_col <- matrix(0, nrow = nrow(ligand_target_matrix), ncol = 1, dimnames = list(rownames(ligand_target_matrix), l))
      ligand_target_matrix <- cbind(ligand_target_matrix, zero_col)
    }
  }
  
  lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]
  
  message("--- Custom Ligand-Receptor Verification ---")
  check_lr <- lr_network %>% filter(receptor == cur_receptor & ligand %in% cur_ligands)
  print(check_lr)
  
  return(list(lr_network = lr_network, ligand_target_matrix = ligand_target_matrix))
}

#' Preprocess Seurat Object for MultiNicheNet Analysis
#'
#' @export
#' @importFrom Seurat sc.CompactSeurat
#' @importFrom SingleCellExperiment as.SingleCellExperiment
#' @importFrom MultiNicheNet alias_to_symbol_SCE makenames_SCE
prepare_sce_object <- function(seurat_obj, organism = "human") {
  
  seurat_obj <- sc.CompactSeurat(seurat_obj)
  seurat_obj[["RNA"]]$counts <- as(seurat_obj[["RNA"]]$counts, "dgCMatrix")
  
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
  rownames(sce) <- make.names(rownames(sce))
  
  sce <- alias_to_symbol_SCE(sce, organism) %>% makenames_SCE()
  
  return(sce)
}

#' Execute the Full MultiNicheNet Analysis Pipeline with Enhanced Parameters
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
  
  if(verbose) message("Subsetting SCE for conditions: ", paste(conditions_keep, collapse = ", "))
  sce <- sce[, SummarizedExperiment::colData(sce)[,group_id] %in% conditions_keep]
  
  senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() %>% as.character()
  receivers_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique() %>% as.character()
  
  abundance_info <- get_abundance_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches
  )
  
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
  
  frq_list <- get_frac_exprs(
    sce = sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id, 
    batches = batches, min_cells = min_cells, fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop
  )
  
  genes_oi <- frq_list$expressed_df %>% dplyr::filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
  sce <- sce[genes_oi, ]
  
  abundance_expression_info <- process_abundance_expression_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, 
    lr_network = lr_network, batches = batches, frq_list = frq_list, abundance_info = abundance_info
  )
  
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

