library(tidyverse)

#####################################
# Create a full, multimodal data matrix from all of the individual modalities.
# Then, create masked data set in mosaic patch setting.
# Inputs:
## view_df: Data frame, showing available modalities for each observation.
## data_dir: location of modality-concatenated data files.
## cancer_folds: vector of cancer types to be included in analysis.
### Needs to match cancer column in full_sample_sheet/view_df
## full_sample_sheet: data frame of dictionary reference of all files downloaded.
## n: maximum number of observations; will subsample if less than number of total subjects.
### Make number small for testing.
## p_max: maximum number of features per modality; ; will subsample if less than number of total columns.
### Make number small for testing.
## v: how many modalities (or vector of indices of which modalities, in order) to keep in simulated data. 
### Make number smaller than full data for testing.
## b: number of different observation patches
## h: views seen per patch, same for all.
## iters: number of simulation iterations.
## wd: directory where simulated mosaic patch observed data should be saved.
mvp_create_tcga <- function(view_df, data_dir, cancer_folds, full_sample_sheet,
                            n, p_max, v, b, h, iters, wd){
  
  # Load individual modalities
  dna_meth_betas <- read.csv(paste0(data_dir, "/all_dna_betas.csv"), row.names = 1)
  mirna_quant <- read.csv(paste0(data_dir, "/all_trans_mirna.csv"), row.names = 1)
  rppa_quant <- read.csv(paste0(data_dir, "/all_RPPA.csv"), row.names = 1)
  genex_quant <- read.csv(paste0(data_dir, "/all_trans_genex.csv"), row.names = 1)
  cnv_mask <- read.csv(paste0(data_dir, "/all_cnv_masked.csv"), row.names = 1)
  cnv_allelic <- read.csv(paste0(data_dir, "/all_cnv_allelic.csv"), row.names = 1)
  cnv_genelevel <- read.csv(paste0(data_dir, "/all_cnv_genelevel.csv"), row.names = 1)
  cnv_seg <- read.csv(paste0(data_dir, "/all_cnv_seg.csv"), row.names = 1)
  
  cnv_full <- rbind(cnv_mask, cnv_allelic, cnv_genelevel, cnv_seg)
  
  p_full <- c(nrow(rppa_quant), nrow(mirna_quant), nrow(genex_quant), nrow(dna_meth_betas),
              nrow(cnv_full))
  
  if(length(p_max) == 1){p_max <- rep(p_max, length(p_full))}
  p <- pmin(p_full, p_max)
  p_mat <- matrix(c(p_full, p), ncol = 2)
  
  setwd(wd)
  
  if(length(v) == 1){
    v_vec <- 1:v
  } else {
    v_vec <- v
  }
  v_l <- length(v_vec)
  
  dir_name <- paste0("./clusts_", paste(cancer_folds[1:clusts], collapse = ""), "_n_", n, "_p_", sum(p_full), "_k_", b, "_views_", paste(v, collapse = ""),  
                     "_vpb_", h)
  
  dir.create(dir_name)
  setwd(dir_name)
  
  n_old <- min(n, nrow(view_df))
  
  # Create information list for return
  return_info <- list()
  clusts <- length(unique(cancer_folds))
  # Iterate
  for(iter in 1:iters){
    n <- n_old
    print(paste0("Iter: ", iter))
    sub_dir_name <-  paste0("./clusts_", paste(cancer_folds[1:clusts], collapse = ""), "_n_", n_old,  "_p_", sum(p[v_vec]), "_k_", b, "_views_", 
                            paste(v, collapse = ""),  "_vpb_", h, "_", iter)
    dir.create(sub_dir_name)
    setwd(sub_dir_name)
    
    # Build list of column names
    p_list <- apply(p_mat, 1, function(x){return(sample(x[1], x[2]))})
    df_column_names <- list()
    df_column_names[[1]] <- c(sapply(isoforms_quant[p_list[[1]], 1], 
                                     function(x){return(paste(x, c("reads_per_million_miRNA_mapped"), sep = "_"))}))
    df_column_names[[2]] <- c(sapply(mirna_quant[p_list[[2]], 1], 
                                     function(x){return(paste(x, c("reads_per_million_miRNA_mapped"), sep = "_"))}))
    df_column_names[[3]] <- c(sapply(genex_quant[p_list[[3]], 1], 
                                     function(x){return(paste(x, c("uunstranded", "stranded_first", "stranded_second",
                                                                   "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"), sep = "_"))}))
    df_column_names[[4]] <- c(sapply(dna_meth_betas[p_list[[4]], 1], 
                                     function(x){return(paste(x, c("measure"), sep = "_"))}))
    df_column_names[[5]] <- c(sapply(cnv_genelevel[p_list[[5]], 1], 
                                     function(x){return(paste(x, c("copy_number"), sep = "_"))}),
                              sapply(cnv_allelic[p_list[[9]], 1], 
                                     function(x){return(paste(x, c("Copy_Number"), sep = "_"))}),
                              sapply(cnv_seg[p_list[[10]], 1], 
                                     function(x){return(paste(x, c("Segment_Mean"), sep = "_"))}),
                              sapply(cnv_mask[p_list[[11]], 1], 
                                     function(x){return(paste(x, c("Segment_Mean"), sep = "_"))})
                              )
    
    load_df <- view_df[, -c(1, 2)]
    df_cnames <- c(colnames(view_df[, -c(1, 2)]))
    df_pnames <- as.list(df_cnames[1:4])
    df_pnames[[5]] <- c(df_cnames[5:8])

    ### Create fully loaded data set, across modalities and observations
    n_samp <- c()
    clust_ass_vec <- c()
    for(cc in 1:clusts){
      n_samp <- c(n_samp, sample(which(view_df$clust == cancer_folds[cc])), floor(n_old / clusts) + floor(cc / clusts))
      clust_ass_vec <- c(clust_ass_vec, rep(cancer_folds[cc], floor(n_old / clusts) + floor(cc / clusts)))
    }
    order_randomizer <- sample(n_old)
    n_samp <- n_samp[order_randomizer]
    clust_ass_vec <- as.numeric(factor(clust_ass_vec[order_randomizer]))
    
    ## Create data frame
    data_mat <- matrix(0, n_old, bp)
    rownames(data_mat) <- load_df$samples[n_samp]
    colnames(data_mat) <- unlist(df_column_names[v_vec])
    names(colnames(data_mat)) <- NULL
    
    ## Build combined dataset
    for(nn in 1:length(n_samp)){
      if(nn %% 1 == 0){
        print(paste0("Sample: ", nn, " / ", length(n_samp)))
      }
      for(vv in 1:length(v_vec)){
        row_ind <- which(full_sample_sheet$data_type %in% df_pnames[[v_vec[vv]]] & full_sample_sheet$Case.ID == view_df$samples[n_samp[nn]])
        if(length(row_ind) > 1){
          row_ind <- which(full_sample_sheet$data_type %in% df_pnames[[v_vec[vv]]] & full_sample_sheet$Case.ID == view_df$samples[n_samp[nn]] & 
                             full_sample_sheet$Sample.Type == "Primary Tumor")
          if(length(row_ind) != 1){
            row_ind <- which(full_sample_sheet$data_type %in% df_pnames[[v_vec[vv]]] & full_sample_sheet$Case.ID == view_df$samples[n_samp[nn]])[1]
          }
        }
        file_name <- paste(data_dir, full_sample_sheet$cancer[row_ind],  full_sample_sheet$dvfold[row_ind], full_sample_sheet$File.ID[row_ind], 
                           full_sample_sheet$File.Name[row_ind] , sep = "/")
        # list.files(data_dir, pattern = full_sample_sheet$File.Name[row_ind], full.names = TRUE, include.dirs = TRUE, recursive = TRUE)
        if(grepl("augmented_star_gene_counts.tsv", file_name)) {
          dff <- read.delim(file_name, comment.char="#")[-c(1:4), ]
          dff <- dff[which(dff$gene_id %in% df_sub_names[[v_vec[vv]]]), ]
        } else if(grepl("level3betas.txt", file_name)) {
          dff <- read.delim(file_name, header = FALSE)
          dff[which(is.na(dff[, 2])), 2] <- 0
          colnames(dff) <- c("V1", "measure")
          dff <- dff[which(dff$V1 %in% df_sub_names[[4]]), ]
        } else if(grepl("mirnas.quantification.txt", file_name)) {
          dff <- read.delim(file_name)
          dff <- dff[which(dff$miRNA_ID %in% df_sub_names[[v_vec[vv]]]), ]
        } else if(grepl("isoforms.quantification.txt", file_name)) {
          dff <- read.delim(file_name)
          dff <- dff[which(paste(dff$miRNA_ID, dff$isoform_coords, sep = "_") %in% df_sub_names[[v_vec[vv]]]), ]
        } else if(grepl("nocnv_grch38.seg.v2", file_name)) {
          dff <- read.delim(file_name)
          dff <- dff[which(paste(dff$Chromosome, dff$Start, sep = "_") %in% df_sub_names[[11]]), ]
        } else if(grepl("\\.grch38.seg.v2", file_name)) {
          dff <- read.delim(file_name)
          dff <- dff[which(paste(dff$Chromosome, dff$Start, sep = "_") %in% df_sub_names[[10]]), ]
        } else if(grepl("allelic_specific.seg", file_name)) {
          dff <- read.delim(file_name, comment.char="#")
          dff <- dff[which(paste(dff$Chromosome, dff$Start, sep = "_") %in% df_sub_names[[9]]), ]
        } else if(grepl("gene_level_copy_number", file_name)) {
          dff <- read.delim(file_name)
          dff <- dff[which(paste(dff$gene_id, dff$gene_name, sep = "_") %in% df_sub_names[[v_vec[vv]]]), ]
          dff$copy_number[which(is.na(dff$copy_number))] <- 0
        } else {
          next
        }
        for(pp in 1:length(df_column_names[[v_vec[vv]]])){
          placement_col <- which(colnames(data_mat) == df_column_names[[v_vec[vv]]][pp])
          if(grepl("augmented_star_gene_counts.tsv", file_name)) {
            full_data_row <- which(sapply(dff$gene_id, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- c("fpkm_unstranded")[
              which(sapply(c("fpkm_unstranded"), 
                           function(x){return(grepl(x, df_column_names[[v_vec[vv]]][pp]))}))]
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("level3betas.txt", file_name)) {
            full_data_row <- which(sapply(dff$V1, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "measure"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("mirnas.quantification.txt", file_name)) {
            full_data_row <- which(sapply(dff$miRNA_ID, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "reads_per_million_miRNA_mapped"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("isoforms.quantification.txt", file_name)) {
            full_data_row <- which(sapply(dff$miRNA_ID, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}) & 
                                     sapply(dff$isoform_coords, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "reads_per_million_miRNA_mapped"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("nocnv_grch38.seg.v2", file_name)) {
            full_data_row <- which(sapply(dff$Chromosome, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}) & 
                                     sapply(dff$Start, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "Segment_Mean"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("\\.grch38.seg.v2", file_name)) {
            full_data_row <- which(sapply(dff$Chromosome, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}) & 
                                     sapply(dff$Start, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "Segment_Mean"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("allelic_specific.seg", file_name)) {
            full_data_row <- which(sapply(dff$Chromosome, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}) & 
                                     sapply(dff$Start, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "Segment_Mean"
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else if(grepl("gene_level_copy_number", file_name)) {
            full_data_row <- which(sapply(dff$gene_name, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}) & 
                                     sapply(dff$gene_id, function(x){return(grepl(paste0(x, "_"), df_column_names[[v_vec[vv]]][pp], fixed = TRUE))}))
            full_data_colname <- "copy_number"
            if(is.na(dff[full_data_row, full_data_colname])){
              print(view_df$samples[n_samp[nn]])
              print(file_name)
              return(dff)
            }
            ff <- try(data_mat[nn, placement_col] <- dff[full_data_row, full_data_colname])
          } else {
            return(dff)
          }
          if(class(ff) == "try-error"){
            return(list(dff, df_column_names[[v_vec[vv]]][pp], vv, pp, file_name))
          }
        }
      }
    }
    
    # Remove 0 variance rows
    z_var_rows <- which(apply(as.matrix(data_mat), 1, var) == 0)
    if(length(z_var_rows) > 0){
      print("0 var row")
      data_mat <- data_mat[-z_var_rows, ]
    }
    
    blocks_n_list <- list()
    overlap_size <- 0
    n <- n - length(z_var_rows)
    n_per <- floor(n / b)
    
    for(kk in 1:b){
      if(kk == b){
        blocks_n_list[[kk]] <- (((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (n))
      } else {
        blocks_n_list[[kk]] <- ((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (kk * n_per - overlap_size * (kk - 1))
      }
    }
    
    ## Remove 0 variance columns
    z_var_cols <- c()
    z_var_colnames <- c()
    for(n1 in 1:length(blocks_n_list)){
      sub_mat <- data_mat[blocks_n_list[[n1]], ]
      z_var_cols  <- c(z_var_cols, which(apply(as.matrix(sub_mat), 2, var) == 0))
      z_var_colnames  <- c(z_var_colnames, colnames(data_mat)[which(apply(as.matrix(sub_mat), 2, var) == 0)])
    }
    
    if(length(z_var_colnames) > 0){
      print("0 var col")
      print(z_var_cols)
      data_mat <- data_mat[, -z_var_cols]
      for(vv in 1:length(v_vec)){
        df_column_names[[v_vec[vv]]] <- setdiff(df_column_names[[v_vec[vv]]], z_var_colnames)
      }
    }
    
    
    ## Create feature set, view list
    blocks_p_list <- list()
    bp <- 0
    for(vv in 1:v_l){
      blocks_p_list[[vv]] <- bp + (1:length(df_column_names[[v_vec[vv]]]))
      bp <- bp + length(df_column_names[[v_vec[vv]]])
    }
    
    ### Create block information
    possible_view_comb <- combn(v_l, h)
    if(b > ncol(possible_view_comb)){
      possible_view_comb <- combn(v_l + 1, h) - 1
      while(TRUE){
        shuf_pan <- sample(1:ncol(possible_view_comb), b)
        if((length(unique(c(possible_view_comb[, shuf_pan]))) == v_l + 1) && 
           (max(table(c(possible_view_comb[, shuf_pan]))) - min(table(c(possible_view_comb[, shuf_pan]))) <= 1)){
          break
        }
      }
      
    } else {
      possible_view_comb <- combn(v_l, h)
      while(TRUE){
        shuf_pan <- sample(1:ncol(possible_view_comb), b)
        if((length(unique(c(possible_view_comb[, shuf_pan]))) == v_l) && 
           (max(table(c(possible_view_comb[, shuf_pan]))) - min(table(c(possible_view_comb[, shuf_pan]))) <= 1)){
          break
        }
      }
      
    }
    
    panels_list <- list()
    obs_panel_list <- list()
    panel_obs_block_list <- list()
    panel_obs_list <- list()
    
    for(bb in 1:b){
      panels_list[[bb]] <- possible_view_comb[, shuf_pan][, bb]
      panels_list[[bb]] <- panels_list[[bb]][panels_list[[bb]] != 0]
      obs_panel_list[[bb]] <- unlist(blocks_p_list[c(panels_list[[bb]])])
    }
    
    for(vv in 1:length(v_vec)){
      panel_obs_block_list[[vv]] <- unique(c(which(possible_view_comb[, shuf_pan] == vv, arr.ind = TRUE)[, 2]))
      panel_obs_list[[vv]] <- unlist(blocks_n_list[c(panel_obs_block_list[[vv]])])
    }
    
    ## Create masked data frame
    masked_dat <- matrix(0, n, bp)
    rownames(masked_dat) <- load_df$samples[n_samp]
    colnames(masked_dat) <- unlist(df_column_names[v_vec])
    names(colnames(masked_dat)) <- NULL
    
    for(bb in 1:b){
      tt <- try(masked_dat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])] <- 
                  data_mat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])])
      if(class(tt) == "try-error"){
        return(list(data_mat, blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])))
      }
      nan_dat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])] <- 
        data_mat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])]
    }
    
    ## Save the info necessary to run comparison methods in code.
    write.csv(masked_dat, "masked_dat.csv", row.names = FALSE)
    write.csv(data_mat, "full_dat.csv", row.names = FALSE)
    
    cfun <- function(L) {
      pad.na <- function(x,len) {
        c(x, rep(NaN,len - length(x)))
      }
      maxlen <- max(sapply(L, length))
      do.call(data.frame, lapply(L, pad.na, len = maxlen))
    }
    
    panel_obs_mat <- matrix(0, n, v_l)
    for(vv in 1:v_l){
      panel_obs_mat[panel_obs_list[[vv]], vv] <- 1
    }
    write.csv(panel_obs_mat, "panel_obs_mat.csv", row.names = FALSE)
    
    panel_times_mat <- matrix(0, ncol(data_mat), v_l)
    for(vv in 1:v_l){
      hh <- try(panel_times_mat[blocks_p_list[[vv]], vv] <- 1)
      if(class(hh) == "try-error"){
        return(list(blocks_p_list, panel_times_mat))
      }
    }
    write.csv(panel_times_mat, "panel_times_mat.csv", row.names = FALSE)
    
    ss2 <- cfun(blocks_n_list)
    colnames(ss2) <- NULL
    write.csv(ss2, "sobs.csv", row.names = FALSE)
    
    ss4 <- cfun(obs_panel_list)
    colnames(ss4) <- NULL
    write.csv(ss4, "stimes.csv", row.names = FALSE)
    
    ss3 <- which(masked_dat != 0)
    write.csv(ss3, "omega.csv", row.names = FALSE)
    
    if(length(z_var_rows > 0)) {
      clust_ass_vec <- clust_ass_vec[-z_var_rows]
    }
    write.csv(clust_ass_vec, "clustass.csv", row.names = FALSE)
    
    # Save info as an .Rdata file for ease of use
    save(clust_ass_vec, blocks_p_list, blocks_n_list, obs_panel_list,
         panel_obs_block_list, panel_obs_list, clusts, df_column_names,
         df_sub_names, panels_list, v_vec, n_samp, view_df, load_df,
         file = "trueinfo.RData")
    
    # Saved output information:
    ## masked_dat: data with patch missingness applied
    ## full_dat: full dat matrix
    ## clust_assignments: cluster label matrix. 
    ## view_features: same as patch_features, kept here for parallel structure with simulations
    ## patch_obs: observations within each patch
    ## patch_features: features within each patch
    ## view_obs:  same as patch_obs, kept here for parallel structure with simulations
    ## omega: vector of unmasked entries
    ## patch_obs_mat: observations within each patch in matrix form, savable as csv. 
    ## patch_features_mat: features within each patch in matrix form, savable as csv. 
    ## view_features_mat: features within each individual patch in matrix form, savable as csv.  
    ## view_obs_mat: observations for which each view is unmasked in matrix form, savable as csv. 
    ## num_clusts: number of clusters, i.e. cancer types
    ## rank: not used; exists for parallel return structure with simulated data
    ## centroid_mats: not used; exists for parallel return structure with simulated data
    ## sigma_mats: not used; exists for parallel return structure with simulated data
    
    return_info[[iter]] <- 
      list(masked_dat = t(masked_dat),
         full_dat = t(data_mat),
         clust_assignments = clust_ass_vec, 
         view_features = blocks_p_list, 
         patch_obs = blocks_n_list, 
         patch_features = blocks_p_list, 
         view_obs = blocks_n_list, 
         omega = ss3,
         view_features_mat = panel_times_mat, 
         patch_obs_mat = ss2, 
         patch_features_mat = ss4, 
         view_obs_mat = panel_obs_mat,
         num_clusts = length(cancer_folds), 
         rank = NULL, 
         centroid_mats = NULL,
         sigma_mats = NULL)
    
    setwd("../")
  }
}

#####################################
#####################################
#####################################

# Sample of how to create synthetic masked data.
# data_dir <- "./output_dir/" #Change to directory where data is saved locally
# save_dir <- "./results_dir/" #Change to directory where new data csv should be saved
# settings_dir <- "./settings_dir/" #Change to directory where settings file is saved locally
# cancer_folds <- c("Lung", "Breast", "Kidney")
# data_view_folds <- c("RPPA", "Transcriptome", "Gene", "Methylation", "CNV")
# view_df <- read.csv(paste0(data_dir, "/view_miss.csv"), row.names = 1) # File created in fuse_data.R.
# full_sample_sheet <- read.csv(paste0(settings_dir, "/full_sample_sheet.csv"), row.names = 1) # File created in fuse_data.R.
# 
# 
# base_n <- 500 # For testing, can set maximum observations to small number to reduce runtime.
# base_p <- rep(100000, 8) # For testing, can set maximum features per modality to small number to reduce runtime.
# base_v <- 5
# base_b <- 4
# base_h <- 3
# iters <- 1 # Use to run multiple simulations at once.
# 
# 
# ttt <- mvp_create_tcga(view_df = view_df, data_dir = data_dir, cancer_folds = cancer_folds,
#                        settings_dir = settings_dir, full_sample_sheet = full_sample_sheet,
#                        n = base_n, p_max = base_p, v = base_v, b = base_b, h = base_h,
#                       iters = iters, wd = save_dir)
