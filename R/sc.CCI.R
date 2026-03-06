#' Build and Extend Ligand-Receptor Database for MultiNicheNet
#'
#' @description
#' This function retrieves the standard NicheNet databases (LR network and ligand-target matrix) 
#' and integrates user-defined ligand-receptor pairs. This is particularly useful for 
#' specialized research targets (e.g., VISTA/VSIR) that may not be fully annotated in 
#' general databases.
#'
#' @param organism A character string indicating the organism: "human" (default) or "mouse".
#' @param custom_ligands A character vector of custom ligand symbols to be added.
#' @param custom_receptor A character string of the custom receptor symbol.
#'
#' @return A list containing the integrated `lr_network` (data frame) and 
#' the `ligand_target_matrix` (numeric matrix).
#' 
#' @export
#' @import tidyverse MultiNicheNet
build_custom_lr_db <- function(organism = "human", 
                               custom_ligands = c("SELPLG", "LGALS9", "MMP13", "LRIG1", "SDC2", "VSIG8", "IGSF11"), 
                               custom_receptor = "VSIR") {
  
  options(timeout = 120)
  
  # 1. Download base databases from Zenodo
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
  
  # 2. Harmonize network format
  lr_network_all <- lr_network_all %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  lr_network <- lr_network_all %>% distinct(ligand, receptor)
  
  # 3. Integrate custom LR pairs
  custom_lr_pairs <- data.frame(ligand = cur_ligands, receptor = cur_receptor)
  lr_network <- bind_rows(lr_network, custom_lr_pairs) %>% distinct()
  
  # 4. Standardize Target Matrix names
  colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>% 
    convert_alias_to_symbols(organism = organism) %>% make.names()
  
  # 5. Handle missing custom ligands in the target matrix by adding zero-columns
  missing_ligands <- setdiff(cur_ligands, colnames(ligand_target_matrix))
  if (length(missing_ligands) > 0) {
    for (l in missing_ligands) {
      zero_col <- matrix(0, nrow = nrow(ligand_target_matrix), ncol = 1, dimnames = list(rownames(ligand_target_matrix), l))
      ligand_target_matrix <- cbind(ligand_target_matrix, zero_col)
    }
  }
  
  # 6. Filter and Subset
  lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
  ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]
  
  message("--- Custom Ligand-Receptor Verification ---")
  check_lr <- lr_network %>% filter(receptor == cur_receptor & ligand %in% cur_ligands)
  print(check_lr)
  
  return(list(lr_network = lr_network, ligand_target_matrix = ligand_target_matrix))
}


#' Preprocess Seurat Object for MultiNicheNet Analysis
#'
#' @description
#' Converts a Seurat object into a SingleCellExperiment (SCE) object, 
#' ensures the matrix is in sparse format, and standardizes gene symbols.
#'
#' @param seurat_obj A Seurat object (e.g., integrated object).
#' @param organism A character string: "human" or "mouse".
#'
#' @return A processed SingleCellExperiment object.
#' 
#' @export
#' @import Seurat SingleCellExperiment
prepare_sce_object <- function(seurat_obj, organism = "human") {
  
  # Optimize object size
  seurat_obj <- sc.CompactSeurat(seurat_obj)
  seurat_obj[["RNA"]]$counts <- as(seurat_obj[["RNA"]]$counts, "dgCMatrix")
  
  # Convert to SCE
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
  rownames(sce) <- make.names(rownames(sce))
  
  # Standardize gene symbols
  sce <- alias_to_symbol_SCE(sce, organism) %>% makenames_SCE()
  
  return(sce)
}


#' Execute the Full MultiNicheNet Analysis Pipeline
#'
#' @description
#' A comprehensive wrapper function that executes the entire MultiNicheNetR workflow. 
#' It covers abundance processing, differential expression (DE), ligand activity prediction, 
#' and interaction prioritization across experimental groups.
#'
#' @param sce A preprocessed SingleCellExperiment object.
#' @param lr_db_list A list containing the `lr_network` and `ligand_target_matrix` (from build_custom_lr_db).
#' @param sample_id Metadata column for sample identifiers.
#' @param group_id Metadata column for condition/grouping (e.g., "condition").
#' @param celltype_id Metadata column for cell type annotations.
#' @param contrasts_oi Character string defining the DE contrasts.
#' @param contrast_tbl Tibble defining the mapping between contrasts and groups.
#' @param covariates Optional metadata columns for DE model adjustment.
#' @param batches Optional metadata columns for batch effect adjustment.
#' @param min_cells Minimum number of cells per cell type to include.
#' @param fraction_cutoff Minimum fraction of cells expressing a gene.
#' @param logFC_threshold Log-fold change threshold for DE genes.
#' @param n.cores Number of CPU cores for parallel processing.
#'
#' @return A list containing prioritized interaction tables and metadata (Lite output).
#' 
#' @export
#' @import MultiNicheNet SummarizedExperiment
run_multinichenet_pipeline <- function(sce, 
                                       lr_db_list, 
                                       sample_id = "orig.ident", 
                                       group_id = "condition", 
                                       celltype_id = "detailed.celltypes", 
                                       contrasts_oi = c("'MDS-Healthy.donor','Healthy.donor-MDS'"),
                                       contrast_tbl = tibble(contrast = c("MDS-Healthy.donor", "Healthy.donor-MDS"), group = c("MDS", "Healthy.donor")),
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
                                       n.cores = 30) {
  
  lr_network <- lr_db_list$lr_network
  ligand_target_matrix <- lr_db_list$ligand_target_matrix
  
  # 1. Define Senders and Receivers
  senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
  receivers_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
  sce <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  
  # 2. Get Abundance & Filtering
  abundance_info <- get_abundance_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, batches = batches
  )
  
  abundance_df_summarized <- abundance_info$abundance_data %>% 
    mutate(keep = as.logical(keep)) %>% 
    group_by(group_id, celltype_id) %>% 
    summarise(samples_present = sum((keep)), .groups = "drop")
  
  total_nr_conditions <- length(unique(SummarizedExperiment::colData(sce)[,group_id]))
  absent_celltypes <- abundance_df_summarized %>% 
    filter(samples_present < 2) %>% 
    group_by(celltype_id) %>% 
    count() %>% 
    filter(n == total_nr_conditions) %>% 
    pull(celltype_id)
  
  senders_oi <- setdiff(senders_oi, absent_celltypes)
  receivers_oi <- setdiff(receivers_oi, absent_celltypes)
  sce <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  
  # 3. Fraction Expression & Gene Filtering
  frq_list <- get_frac_exprs(
    sce = sce, sample_id = sample_id, celltype_id = celltype_id, group_id = group_id, 
    batches = batches, min_cells = min_cells, fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop
  )
  
  genes_oi <- frq_list$expressed_df %>% filter(expressed == TRUE) %>% pull(gene) %>% unique() 
  sce <- sce[genes_oi, ]
  
  # 4. Processing Info
  abundance_expression_info <- process_abundance_expression_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, 
    lr_network = lr_network, batches = batches, frq_list = frq_list, abundance_info = abundance_info
  )
  
  # 5. Differential Expression Analysis
  DE_info <- get_DE_info(
    sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
    batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, 
    min_cells = min_cells, expressed_df = frq_list$expressed_df
  )
  
  if(empirical_pval) {
    DE_info_emp <- get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
    celltype_de <- DE_info_emp$de_output_tidy_emp %>% 
      select(-p_val, -p_adj) %>% 
      rename(p_val = p_emp, p_adj = p_adj_emp)
  } else {
    celltype_de <- DE_info$celltype_de$de_output_tidy
  } 
  
  sender_receiver_de <- combine_sender_receiver_de(
    sender_de = celltype_de, receiver_de = celltype_de, senders_oi = senders_oi, 
    receivers_oi = receivers_oi, lr_network = lr_network
  )
  
  # 6. Ligand Activity Prediction
  ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
    get_ligand_activities_targets_DEgenes(
      receiver_de = celltype_de, 
      receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
      ligand_target_matrix = ligand_target_matrix, 
      logFC_threshold = logFC_threshold,
      p_val_threshold = p_val_threshold, 
      p_val_adj = p_val_adj, 
      top_n_target = top_n_target,
      verbose = TRUE, 
      n.cores = n.cores
    )
  ))
  
  # 7. Metadata Grouping & Prioritization
  sender_receiver_tbl <- sender_receiver_de %>% distinct(sender, receiver)
  metadata_combined <- SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
  
  grouping_tbl <- metadata_combined[, c(sample_id, group_id, if(!is.na(batches)) batches else NULL)] %>% 
    distinct()
  colnames(grouping_tbl)[1:2] <- c("sample", "group")
  
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
  
  # 8. Inference Correlation
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
  
  # Final Consolidated Output
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
  
  return(make_lite_output(multinichenet_output))
}