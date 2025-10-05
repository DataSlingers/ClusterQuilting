#############################
#############################
#############################
library(pracma)

# n: total observations
# p: features per view. vector or single number if all same
# v: number of views
# b: number of different patches
# h: views seen per patch, same for all.
# rank_sim: simulation of rank
# means_cons: average distance between clusters.
# sd_cons: standard deviation of noise.
mos_create <- function(n, p, v, b, h, clusts, rank_sim, mean_cons, sd_cons){
  
  
  # Make true underlying cluster assignments
  clust_ass <- matrix(sample(n), ncol = clusts)
  clust_ass_vec <- rep(NA, n)
  for(ii in 1:clusts){
    clust_ass_vec[clust_ass[, ii]] <- ii
  }
  
  clust_means <- list()
  clust_sigmas <- list()
  
  # List of views
  blocks_p_list <- list()
  if(length(p) == 1){
    for(vv in 1:v){
      blocks_p_list[[vv]] <- c(((vv - 1) * p + 1): (vv * p))
    }
  } else {
    for(vv in 1:v){
      blocks_p_list[[vv]] <- c((sum(p[0:(vv-1)]) + 1):(sum(p[0:(vv)])))
    }
  }
  
  # Create data matrix
  data_mat <- matrix(NA, n, sum(p))
  
  # Create centroid matrices
  mc <- c(mean_cons, 0, -mean_cons)
  create_matrix <- paste0("clust_lr2 <- as.matrix(expand.grid(", paste(rep("mc", rank_sim), collapse = ", "), "))")
  eval(parse(text = create_matrix))
  iter <- 0
  while(iter < 100){
    iter <- iter + 1
    clust_lr <- clust_lr2[sample(nrow(clust_lr2), clusts), , drop = FALSE] 
    if(svd(clust_lr)$d[length(svd(clust_lr)$d)] > 10^-10){
      break
    }
  }
  clust_lr = clust_lr + rnorm(nrow(clust_lr) * ncol(clust_lr), 0, sd_cons / 2)
  
  ## Create centroid matrix for individual view
  for(vv in 1:v){
    if(length(p) == 1){
      p_size <- p
    } else {
      p_size <- p[vv]
    }
    clust_means[[vv]] <- list()
    clust_sigmas[[vv]] <- list()
    feature_lr <- t(t(randortho(p_size))[, 1:rank_sim])
    clust_means_mat <- clust_lr %*% feature_lr
    for(cc in 1:clusts){
      clust_means[[vv]][[cc]] <- clust_means_mat[cc, ]
      clust_sigmas[[vv]][[cc]] <- diag((p_size))

    }
    
    # Simulate data for individual view
    # Change distribution as necessary
    for(cc in 1:clusts){
      for(dd in 1:length(clust_ass[, cc])){
        row_ass <- clust_ass[, cc]
        data_mat[row_ass[dd], blocks_p_list[[vv]]] <- 
          rnorm(length(row_ass[dd]) * length(blocks_p_list[[vv]]), clust_means[[vv]][[cc]], sd_cons)
      }
    }
  }
  # Selected views for each patch
  possible_view_comb <- combn(v, h)
  if(b > ncol(possible_view_comb)){
    possible_view_comb <- combn(v + 1, h) - 1
    while(TRUE){
      shuf_pan <- sample(1:ncol(possible_view_comb), b)
      if((length(unique(c(possible_view_comb[, shuf_pan]))) == v + 1) && 
         (max(table(c(possible_view_comb[, shuf_pan]))) - min(table(c(possible_view_comb[, shuf_pan]))) <= 1)){
        break
      }
    }
    
  } else {
    possible_view_comb <- combn(v, h)
    while(TRUE){
      shuf_pan <- sample(1:ncol(possible_view_comb), b)
      if((length(unique(c(possible_view_comb[, shuf_pan]))) == v) && 
         (max(table(c(possible_view_comb[, shuf_pan]))) - min(table(c(possible_view_comb[, shuf_pan]))) <= 1)){
        break
      }
    }
    
  }
  
  
  # List of obs blocks
  blocks_n_list <- list()
  panels_list <- list()
  obs_panel_list <- list()
  panel_obs_block_list <- list()
  panel_obs_list <- list()
  overlap_size <- 0
  n_per <- floor(n / b)
  for(kk in 1:b){
    if(kk == b){
      blocks_n_list[[kk]] <- (((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (n))
    } else {
      blocks_n_list[[kk]] <- ((kk - 1) * n_per + 1 - overlap_size * (kk - 1)) : (kk * n_per - overlap_size * (kk - 1))
    }
    panels_list[[kk]] <- possible_view_comb[, shuf_pan][, kk]
    panels_list[[kk]] <- panels_list[[kk]][panels_list[[kk]] != 0]
  }
  
  masked_dat <- matrix(0, n, sum(p))
  for(bb in 1:b){
    obs_panel_list[[bb]] <- unlist(blocks_p_list[c(panels_list[[bb]])])
    masked_dat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])] <- 
      data_mat[blocks_n_list[[bb]], unlist(blocks_p_list[c(panels_list[[bb]])])]
  }
  
  for(vv in 1:v){
    panel_obs_block_list[[vv]] <- unique(c(which(possible_view_comb[, shuf_pan] == vv, arr.ind = TRUE)[, 2]))
    panel_obs_list[[vv]] <- unlist(blocks_n_list[c(panel_obs_block_list[[vv]])])
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
  panel_obs_mat <- matrix(0, n, v)
  for(vv in 1:v){
    panel_obs_mat[panel_obs_list[[vv]], vv] <- 1
  }
  # write.csv(panel_obs_mat, "panel_obs_mat.csv", row.names = FALSE)
  
  panel_times_mat <- matrix(0, sum(p), v)
  for(vv in 1:v){
    panel_times_mat[blocks_p_list[[vv]], vv] <- 1
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
  ## view_features: features within each individual view. 
  ## patch_obs: observations within each patch
  ## patch_features: features within each patch
  ## patch_views: views within each patch
  ## view_obs: observations for which each view is unmasked
  ## omega: vector of unmasked entries
  ## view_features_mat: features within each individual view in matrix form, savable as csv.  
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
              view_features = blocks_p_list, 
              patch_obs = blocks_n_list, 
              patch_features = obs_panel_list, 
              patch_views = panel_obs_block_list, 
              view_obs = panel_obs_list, 
              omega = ss3,
              view_features_mat = panel_times_mat, 
              patch_obs_mat = ss2, 
              patch_features_mat = ss4, 
              view_obs_mat = panel_obs_mat,
              num_clusts = clusts, 
              rank = rank_sim, 
              centroid_mats = clust_means,
              sigma_mats = clust_sigmas))
    
}

set.seed(42)
base_n <- 840
base_v <- 12
base_p <- rep(50, base_v)
base_b <- 4
base_h <- 6
base_clusts <- 3
base_rank <- 2
base_mean <- 4.5
sdd <- 1

# test_output <- mos_create(n = base_n, p = base_p, v = base_v, b = base_b, h = base_h, clusts = base_clusts, rank_sim = base_rank, mean_cons = base_mean, sdd)
