source("./SpectralClustering.R")
source("./HeterogeneousSNRPatchOrdering.R")
source("./HomogenousSNRPatchOrdering.R")
source("./ClusterQuitling.R")
source("./PredictionValidation.R")

## Data driven selection of rank and number of clusters for MICRONS data.
## Inputs:
# data_dir: Overall directory of where data and info are currently stored.
# Should contain .Rdata created by the make_full.R script in this directory.
# est_dir: Directory where results across different methods & HP settings are stored.
# Results for different HP settings is in a separate folder for each method.
# rank_solve: the ranks to be tested in hyperparameter selection for Cluster Quillting
# clust_solve: the number of clusters in hyperparameter selection for Cluster Quillting
# Will test all valid combinations of the two (rank <= clusters).
# homo_SNR: homogeneous SNR assumption for patch ordering?

compile_ari_microns <- function(data_dir, est_dir, rank_solve = NULL, clust_solve = NULL,
                           homo_SNR = FALSE){
    
  # Compile evaluation metrics (ARI) across iterations of simulation setting
    
  setwd(data_dir)
  masked_dat <- read.csv(paste0("masked_dat_fl.csv"), row.names = NULL)
  load("./trueinfo.RData") # Reload of info from simulation code
  
  # Derive combinations of rank and clusters to try.
  param_df <- expand.grid(x = rank_solve, y = clust_solve)
  param_df <- param_df[param_df$x <= param_df$y, ]
  if(nrow(param_df) == 0){
    stop("Must be at least one rank lower than number of clusters.")
  }
  rank_est <- param_df[, 1]
  clust_est <- param_df[, 2]
  
  ## Assume each method has its own subfolder of results, one file for each HP setting.
  ## Assume that directory is named after the method.
  ## Change to correct location otherwise.
  
  setwd(est_dir)
  
  # IMSC_AGL
  ls3 <- list.files("./IMSCAGL/", full.names = TRUE)
  best_ii <- 0
  best_ari <- 0
  for(ii in 1:length(ls3)){
    IMSCAGL_res <- read.table(ls3[ii])
    dd_ari <- pred_val(dataset, c(IMSCAGL_res$V1), 0.2)
    if(dd_ari > best_ari){
      best_ari <- dd_ari
      best_ii <- ii
    }
  }
  IMSCAGL_res <- read.table(paste0(ls3[best_ii]))
  write.csv(IMSCAGL_res, paste0(est_dir, "/IMSCAGL_res.csv"), row.names = FALSE)
  
  # IMG
  ls3 <- list.files("./IMG/", full.names = TRUE)
  best_ii <- 0
  best_ari <- 0
  for(ii in 1:length(ls3)){
    IMG_res <- read.table(ls3[ii])
    dd_ari <- pred_val(dataset, c(IMG_res$V1), 0.2)
    if(dd_ari > best_ari){
      best_ari <- dd_ari
      best_ii <- ii
    }
  }
  IMG_res <- read.table(paste0(ls3[best_ii]))
  write.csv(IMG_res, paste0(est_dir, "/IMG_res.csv"), row.names = FALSE)
  
  # DAIMC
  ls3 <- list.files("./DAIMC/", full.names = TRUE)
  best_ii <- 0
  best_ari <- 0
  for(ii in 1:length(ls3)){
    DAIMC_res <- read.table(ls3[ii])
    dd_ari <- pred_val(dataset, c(DAIMC_res$V1), 0.2)
    if(dd_ari > best_ari){
      best_ari <- dd_ari
      best_ii <- ii
    }
  }
  DAIMC_res <- read.table(paste0(ls3[best_ii]))
  write.csv(DAIMC_res, paste0(est_dir, "/DAIMC_res.csv"), row.names = FALSE)
  
  # OPIMC
  ls3 <- list.files("./OPIMC/", full.names = TRUE)
  best_ii <- 0
  best_ari <- 0
  for(ii in 1:length(ls3)){
    OPIMC_res <- read.table(ls3[ii])
    dd_ari <- pred_val(dataset, c(OPIMC_res$V1), 0.2)
    if(dd_ari > best_ari){
      best_ari <- dd_ari
      best_ii <- ii
    }
  }
  OPIMC_res <- read.table(paste0(ls3[best_ii]))
  write.csv(OPIMC_res, paste0(est_dir, "/OPIMC_res.csv"), row.names = FALSE)
  
  # NN
  ls3 <- list.files("./NN/", full.names = TRUE)
  best_nn <- 0
  best_ii <- 0
  best_ari <- 0
  for(ii in 1:length(ls3)){
    for(nn in 1:length(rank_est)){
      test_datimpute <- read.csv(paste0(ls3[ii]), header=FALSE)
      test_sc_nn <- spec_clust(t(test_datimpute), rank_est[nn], clust_est[nn])
      dd_ari <- pred_val(dataset, c(test_sc_nn), 0.2)
      if(dd_ari > best_ari){
        best_ari <- dd_ari
        best_nn <- nn
        best_ii <- ii
      }
    }
  }
  test_datimpute <- read.csv(paste0(ls3[best_ii]), header=FALSE)
  test_sc_nn <- spec_clust(t(test_datimpute), rank_est[best_nn], clust_est[best_nn])
  write.csv(test_sc_nn, paste0(est_dir, "/NN_res.csv"), row.names = FALSE)
  
  # Cluster Quilting
  best_nn <- 0
  best_ari <- 0
  for(nn in 1:length(rank_est)){
    if(homo_SNR){
      tpp <- traverse_path(masked_dat, panel_obs_list, blocks_p_list)
    } else {
      tpp <- traverse_path_het(masked_dat, panel_obs_list, blocks_p_list, rank_est[nn])
    }
    cq_est <- ClusterQuilting(masked_dat, blocks_p_list, panel_obs_list, 
                              rank_est[nn], clust_est[nn], tpp)
    dd_ari <- pred_val(dataset, c(cq_est), 0.2)
    if(dd_ari > best_ari){
      best_ari <- dd_ari
      best_nn <- nn
    }
  }
  
  if(homo_SNR){
    tpp <- traverse_path(masked_dat, panel_obs_list, blocks_p_list)
  } else {
    tpp <- traverse_path_het(masked_dat, panel_obs_list, blocks_p_list, rank_est[best_nn])
  }
  cq_est <- ClusterQuilting(t(masked_dat), blocks_p_list, panel_obs_list, 
                            rank_est[best_nn], clust_est[best_nn], tpp)
  write.csv(cq_est, paste0(est_dir, "/CQ_res.csv"), row.names = FALSE)
  
}

  
