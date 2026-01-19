library(mclust)
library(pracma)
library(maotai)

source("./SpectralClustering.R")
source("./HeterogeneousSNRPatchOrdering.R")
source("./HomogenousSNRPatchOrdering.R")
source("./ClusterQuitling.R")
source("./PredictionValidation.R")

## Compiles mean and std.err of ARI across multiple iterations of a single simulation setting.
## Oracle selection of rank and number of clusters.
## Also runs spectral clustering on full data + NN imputed data to get comparison clustering results.
## Also runs Cluster Quilting to get the cluster estimates and performance results.
## Inputs:
# sim_dir: Overall directory of where data are currently stored.
# Assumes results_dir is an overall directory of a full simulation study.
# Each subfolder contains subsubfolders for all of the runs across an individual simulation setting.
# Each subsubfolder contains a single iteration of a single simulation setting.
# homo_SNR: homogeneous SNR assumption for patch ordering?
compile_ari_oracle <- function(sim_dir, homo_SNR = FALSE){
  # Get list of subfolders, containing all iterations of a particular simulation setting
  lss <- list.dirs(sim_dir, full.names = TRUE)
  for(simfold in 1:length(lss)){
    # Get list of subsubfolders, containing single iteration for a single setting
    lss2 <- list.dirs(lss[simfold], full.names = TRUE)

    # Create save arrays
    qq <- c()
    qq2 <- c()
    qq3 <- c()
    qq4 <- c()
    qq5 <- c()
    qq6 <- c()
    qq8 <- c()
    cqq <- c()
    

    # Compile evaluation metrics (ARI) across iterations of simulation setting
    for (subfold in 1:length(lss2)){
      setwd(lss2[subfold])
      data_mat <- read.csv("./full_dat.csv")
      masked_dat <- read.csv("./masked_dat.csv")
      load("./trueinfo.RData") # Reload of info from simulation code
      
      # Defaults for oracle.
      rank_est <- rank_sim
      clust_est <- clusts
      
      # Spectral clustering on unmasked data
      test_sc <- spec_clust(t(data_mat), rank_sim, TRUE)
      qq <- c(qq, adjustedRandIndex(clust_ass_vec, test_sc))
      
      # IMSC_AGL
      IMSCAGL_res <- read.table("./IMSCAGL_res.csv")
      qq2 <- c(qq2, adjustedRandIndex(clust_ass_vec, c(IMSCAGL_res$V1)))
      
      # IMG
      IMG_res <- read.table("./IMG_res.csv")
      qq3 <- c(qq3, adjustedRandIndex(clust_ass_vec, c(IMG_res$V1)))
      
      # DAIMC
      DAIMC_res <- read.table("./DAIMC_res.csv")
      qq5 <- c(qq5, adjustedRandIndex(clust_ass_vec, c(DAIMC_res$V1)))
      
      # OPIMC
      OPIMC_res <- read.table("./OPIMC_res.csv")
      qq6 <- c(qq6, adjustedRandIndex(clust_ass_vec, c(OPIMC_res$V1)))
      
      # NN
      test_datimpute <- read.csv("./nn_impute.csv", header=FALSE)
      test_sc_nn <- spec_clust(t(test_datimpute), rank_est, clust_est)
      qq8 <- c(qq8, adjustedRandIndex(clust_ass_vec, test_sc_nn))
      
      # Cluster Quilting
      if(homo_SNR){
        tpp <- traverse_path(masked_dat, panel_obs_list, blocks_p_list)
      } else {
        tpp <- traverse_path_het(masked_dat, panel_obs_list, blocks_p_list, rank_est)
      }
      cq_est <- ClusterQuilting(t(masked_dat), blocks_p_list, panel_obs_list, 
                                rank_est, clust_est, tpp)
      write.csv(cq_est, "./CQ_res.csv", row.names = NULL)
      cqq <- c(cqq, adjustedRandIndex(clust_ass_vec, cq_est))
    }
    mod_names <- c("Full Data", "IMSC_AGL", "IMG", "DAIMC", "OPIMC", "NN Data", "ClustQuilt")
    
    # Calculate average, std.err. of ARI for all methods across iterations of sim parameter settings.
    ari_res <- c(mean(qq, na.rm = TRUE), mean(qq2, na.rm = TRUE),
                 mean(qq3, na.rm = TRUE),
                 mean(qq5, na.rm = TRUE), mean(qq6, na.rm = TRUE),
                 mean(qq8, na.rm = TRUE), mean(cqq, na.rm = TRUE))
    sd_res <- c(sd(qq, na.rm = TRUE), sd(qq2, na.rm = TRUE),
                sd(qq3, na.rm = TRUE), 
                sd(qq5, na.rm = TRUE), sd(qq6, na.rm = TRUE),
                sd(qq8, na.rm = TRUE), sd(cqq, na.rm = TRUE))
    sd_res <- sd_res / sqrt(length(lss2))
    full_res <- as.data.frame(cbind(qq, qq2, qq3, qq5, qq6, qq8, cqq))
    colnames(full_res) <- mod_names
    save_df <- data.frame(mods = mod_names,
                          ari = ari_res,
                          sd = sd_res)
    
    # Save results for simulation setting in subfolder (not subsubfolder or folder)
    write.csv(save_df, paste0(lss[simfold], "/ari_results.csv"), row.names = FALSE)
    print(save_df)
    write.csv(full_res, paste0(lss[simfold], "/full_results.csv"), row.names = FALSE)
  }
  
}

################################################################################

## Compiles mean and std.err of ARI across multiple iterations of a single simulation setting.
## Data driven selection of rank and number of clusters.
## Also runs spectral clustering on full data + NN imputed data to get comparison clustering results.
## Also runs Cluster Quilting to get the cluster estimates and performance results.
## Inputs:
# sim_dir: Overall directory of where data are currently stored.
# Assumes results_dir is an overall directory of a full simulation study.
# Each subfolder contains subsubfolders for all of the runs across an individual simulation setting.
# Each subsubfolder contains a single iteration of a single simulation setting.
# Results for different HP settings is in a separate folder for each method.
# rank_solve: the ranks to be tested in hyperparameter selection for Cluster Quillting
# clust_solve: the number of clusters in hyperparameter selection for Cluster Quillting
# Will test all valid combinations of the two (rank <= clusters).
# homo_SNR: homogeneous SNR assumption for patch ordering?
compile_ari_dd <- function(sim_dir, rank_solve = NULL, clust_solve = NULL,
                               homo_SNR = FALSE){
  # Get list of subfolders, containing all iterations of a particular simulation setting
  lss <- list.dirs(sim_dir, full.names = TRUE)
  for(simfold in 1:length(lss)){
    # Get list of subsubfolders, containing single iteration for a single setting
    lss2 <- list.dirs(lss[simfold], full.names = TRUE)
    
    # Create save arrays
    qq <- c()
    qq2 <- c()
    qq3 <- c()
    qq4 <- c()
    qq5 <- c()
    qq6 <- c()
    qq8 <- c()
    cqq <- c()
    
    
    # Compile evaluation metrics (ARI) across iterations of simulation setting
    for (subfold in 1:length(lss2)){
      setwd(lss2[subfold])
      data_mat <- read.csv("./full_dat.csv")
      masked_dat <- read.csv("./masked_dat.csv")
      load("./trueinfo.RData") # Reload of info from simulation code
      
      # Defaults to oracle if unset, but should be user-specified.
      if(is.null(rank_solve)){
        rank_solve <- rank_sim
      } 
      
      if(is.null(clust_solve)){
        clust_solve <- clusts
      } 
      
      # Derive combinations of rank and clusters to try.
      param_df <- expand.grid(x = rank_solve, y = clust_solve)
      param_df <- param_df[param_df$x <= param_df$y, ]
      if(nrow(param_df) == 0){
        stop("Must be at least one rank lower than number of clusters.")
      }
      rank_est <- param_df[, 1]
      clust_est <- param_df[, 2]
      
      
      # Spectral clustering on unmasked data
      test_sc <- spec_clust(t(data_mat), rank_sim, TRUE)
      qq <- c(qq, adjustedRandIndex(clust_ass_vec, test_sc))
      
      ## Assume each method has its own subfolder of results, one file for each HP setting.
      ## Assume that directory is named after the method.
      ## Change to correct location otherwise.
      
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
      qq2 <- c(qq2, adjustedRandIndex(clust_ass_vec, c(IMSCAGL_res$V1)))
      
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
      qq3 <- c(qq3, adjustedRandIndex(clust_ass_vec, IMG_res$V1))
      
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
      qq5 <- c(qq5, adjustedRandIndex(clust_ass_vec, c(DAIMC_res$V1)))
      
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
      qq6 <- c(qq6, adjustedRandIndex(clust_ass_vec, c(OPIMC_res$V1)))
      
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
      qq8 <- c(qq8, adjustedRandIndex(clust_ass_vec, test_sc_nn))
      
      # Cluster Quilting
      best_nn <- 0
      best_ari <- 0
      for(nn in 1:length(rank_est)){
        if(homo_SNR){
          tpp <- traverse_path(masked_dat, panel_obs_list, blocks_p_list, blocks_n_list)
        } else {
          tpp <- traverse_path_het(masked_dat, panel_obs_list, blocks_p_list, blocks_n_list, rank_est[nn])
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
      write.csv(cq_est, "./CQ_res.csv", row.names = NULL)
      cqq <- c(cqq, adjustedRandIndex(clust_ass_vec, cq_est))
    }
    mod_names <- c("Full Data", "IMSC_AGL", "IMG", "DAIMC", "OPIMC", "NN Data", "ClustQuilt")
    
    # Calculate average, std.err. of ARI for all methods across iterations of sim parameter settings.
    ari_res <- c(mean(qq, na.rm = TRUE), mean(qq2, na.rm = TRUE),
                 mean(qq3, na.rm = TRUE),
                 mean(qq5, na.rm = TRUE), mean(qq6, na.rm = TRUE),
                 mean(qq8, na.rm = TRUE), mean(cqq, na.rm = TRUE))
    sd_res <- c(sd(qq, na.rm = TRUE), sd(qq2, na.rm = TRUE),
                sd(qq3, na.rm = TRUE), 
                sd(qq5, na.rm = TRUE), sd(qq6, na.rm = TRUE),
                sd(qq8, na.rm = TRUE), sd(cqq, na.rm = TRUE))
    sd_res <- sd_res / sqrt(length(lss2))
    full_res <- as.data.frame(cbind(qq, qq2, qq3, qq5, qq6, qq8, cqq))
    colnames(full_res) <- mod_names
    save_df <- data.frame(mods = mod_names,
                          ari = ari_res,
                          sd = sd_res)
    
    # Save results for simulation setting in subfolder (not subsubfolder or folder)
    write.csv(save_df, paste0(lss[simfold], "/ari_results_dd.csv"), row.names = FALSE)
    print(save_df)
    write.csv(full_res, paste0(lss[simfold], "/full_results_dd.csv"), row.names = FALSE)
  }
  
}

