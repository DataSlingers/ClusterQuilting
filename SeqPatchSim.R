#############################
#############################
#############################

library(pracma)
library(mvtnorm)

# n: total observations
# p: total features
# n_per: number of observations per block.
# b: number of different patches
# h: views seen per patch; must be 1 for this simulation
# clusts: number of clusters
# rank_sim: simulation of rank
# means_cons: average distance between clusters.
# sd_cons: standard deviation of noise.

seq_create <- function(n, p, n_per, b, clusts, rank_sim, mean_cons, sd_cons){
  
  
  # Make true underlying cluster assignments
  clust_ass <- matrix(sample(n), ncol = clusts)
  clust_ass_vec <- rep(NA, n)
  for(ii in 1:clusts){
    clust_ass_vec[clust_ass[, ii]] <- ii
  }
  
  clust_means <- list()
  clust_sigmas <- list()

  
  # Create data matrix
  data_mat <- matrix(NA, n, p)
  
  mc <- c(mean_cons, 0, -mean_cons)
  create_matrix <- paste0("clust_lr2 <- as.matrix(expand.grid(", paste(rep("mc", rank_sim), collapse = ", "), "))")
  eval(parse(text = create_matrix))
  iter <- 0
  print(dim(clust_lr2))
  while(iter < 100){
    iter <- iter + 1
    clust_lr <- clust_lr2[sample(nrow(clust_lr2), clusts), , drop = FALSE] 
    if(svd(clust_lr)$d[length(svd(clust_lr)$d)] > 10^-10){
      break
    }
  }
  clust_lr = clust_lr + rnorm(nrow(clust_lr) * ncol(clust_lr), 0, sd_cons / 4)
  
  clust_means <- list()
  clust_sigmas <- list()
  feature_lr <- t(t(randortho(p))[, 1:rank_sim])
  clust_means_mat <- clust_lr %*% feature_lr
  for(cc in 1:clusts){
    clust_means[[cc]] <- clust_means_mat[cc, ]
    clust_sigmas[[cc]] <- diag(p) * sd_cons
  }
  

  data_mat <- matrix(NA, n, p)
  for(ii in 1:clusts){
    for(jj in 1:length(clust_ass[, ii])){
      data_mat[clust_ass[jj, ii], ] <- rmvnorm(1, clust_means[[ii]], clust_sigmas[[ii]])
    }
  }
  
  ### Data masking
  # Selected views for each patch
  blocks_n_list <- list()
  blocks_p_list <- list()
  overlap_size <- ceiling((n_per * b - n) / (b - 1))
  for(kk in 1:b){
    if(kk == b){
      blocks_n_list[[kk]] <- (((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (n))
      blocks_p_list[[kk]] <-  (((kk - 1) * floor(p / b) + 1) : (p))
    } else {
      blocks_n_list[[kk]] <- ((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (kk * n_per - overlap_size * (kk - 1))
      blocks_p_list[[kk]] <-  (((kk - 1) * floor(p / b) + 1) : ((kk) * floor(p / b)))
    }
  }
  
  # List of obs blocks
  obs_panel_list <- list()
  panel_obs_block_list <- list()
  panel_obs_list <- list()
  
  masked_dat <- matrix(0, n, p)
  for(bb in 1:b){
    obs_panel_list[[bb]] <- blocks_p_list[[bb]]
    masked_dat[blocks_n_list[[bb]], blocks_p_list[[bb]]] <- data_mat[blocks_n_list[[bb]], blocks_p_list[[bb]]]
  }

  for(bb in 1:b){
    panel_obs_block_list[[bb]] <- bb
    panel_obs_list[[bb]] <- blocks_n_list[[bb]]
  }
  # Create versions of matrices compatible with saving as csvs.
  ## Used for fitting comparison methods in Matlab.
  cfun <- function(L) {
    pad.na <- function(x,len) {
      c(x,rep(NaN,len-length(x)))
    }
    maxlen <- max(sapply(L,length))
    do.call(data.frame,lapply(L,pad.na,len=maxlen))
  }
  
  # 
  panel_obs_mat <- matrix(0, n, b)
  for(bb in 1:b){
    panel_obs_mat[panel_obs_list[[bb]], bb] <- 1
  }
  # write.csv(panel_obs_mat, "panel_obs_mat.csv", row.names = FALSE)
  
  panel_times_mat <- matrix(0, p, b)
  for(bb in 1:b){
    panel_times_mat[blocks_p_list[[bb]], bb] <- 1
  }
  # write.csv(panel_times_mat, "panel_times_mat.csv", row.names = FALSE)
  
  ss2 <- cfun(blocks_n_list)
  colnames(ss2) <- NULL
  # write.csv(ss2, "sobs.csv", row.names = FALSE)
  
  ss4 <- cfun(obs_panel_list)
  colnames(ss4) <- NULL
  # write.csv(ss4, "stimes.csv", row.names = FALSE)
  
  ss3 <- which(masked_dat != 0)
  # write.csv(ss3, "omega.csv", row.names = FALSE)
  
  # Outputs:
  ## masked_dat: data with patch missingness applied
  ## full_dat: full dat matrix
  ## clust_assignments: cluster label matrix. 
  ## patch_obs: observations within each patch
  ## patch_features: features within each patch
  ## view_obs: observations for which each view is unmasked
  ## omega: vector of unmasked entries
  ## patch_obs_mat: observations within each patch in matrix form, savable as csv. 
  ## patch_features_mat: features within each patch in matrix form, savable as csv. 
  ## view_obs_mat: observations for which each view is unmasked in matrix form, savable as csv. 
  ## num_clusts: number of clisters 
  ## rank: rank of simulation
  ## centroid_mats: centroid matrices
  ## sigma_mats: covariance matrices
  return(list(masked_dat = t(masked_dat),
              full_dat = t(data_mat),
              clust_assignments = clust_ass_vec, 
              patch_obs = blocks_n_list, 
              patch_features = blocks_p_list, 
              view_obs = panel_obs_list, 
              omega = ss3,
              patch_obs_mat = ss2, 
              patch_features_mat = ss4, 
              view_obs_mat = panel_obs_mat,
              num_clusts = clusts, 
              rank = rank_sim, 
              centroid_mats = clust_means,
              sigma_mats = clust_sigmas))
  
}

# set.seed(44)
# p_base <- 100
# n_base <- 710
# n_per_base <- 210
# b_base <- 4
# mean_cons_base <- 4.5
# sd_cons_base <- 1
# base_clusts <- 3
# base_rank <- 2

# test_output <- seq_create(n_base, p_base, n_per_base, b_base, base_clusts, base_rank,
#                           mean_cons_base, sd_cons_base)