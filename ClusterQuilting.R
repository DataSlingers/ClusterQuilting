
library(RMTstat)

#### Cluster Quilting function
## dat: raw data matrix with features as columns and observations as rows.
## p_obs: list of arrays of features for each patch.
## n_obs: list of arrays of observations for each patch.
## traverse_order: order in which patches should be merged.
cq_r <- function(data_matrix, p_obs, n_obs, rr, traverse_order){
  
  
  # Matrix imputation
  HH <- matrix(0, nrow(data_matrix), rr)
  VV <- matrix(0, ncol(data_matrix), rr)
  p_block <- p_obs[[traverse_order[1]]]
  n_block <- n_obs[[traverse_order[1]]]
  block_dat <- data_matrix[p_block, n_block]
  svd_block <- svd(block_dat, rr, rr)
  prev_n <- n_obs[[traverse_order[1]]]
  VV[n_block, ] <- svd_block$v
  HH[p_block, ] <- svd_block$u %*% diag(svd_block$d[1:rr], nrow = rr, ncol = rr)
  
  for(bb in 2:length(traverse_order)){
    p_block <- p_obs[[traverse_order[bb]]]
    n_block <- n_obs[[traverse_order[bb]]]
    block_dat <- data_matrix[p_block, n_block]
    svd_block <- svd(block_dat, rr, rr)
    overlaps <- prev_n[which(prev_n %in% n_obs[[traverse_order[bb]]])]
    overlaps2 <- which(n_obs[[traverse_order[bb]]] %in% prev_n)
    setdiff1 <- setdiff(n_obs[[traverse_order[bb]]], prev_n)
    setdiff2 <- which(!(n_obs[[traverse_order[bb]]] %in% prev_n))
    G_m <- (pinv(t(svd_block$v[overlaps2, ]) %*% (svd_block$v[overlaps2, ])) %*% (t(svd_block$v[overlaps2, ]) %*% VV[overlaps, ])) 
    VV[setdiff1, ] <- svd_block$v[setdiff2, ] %*% G_m
    HH[p_block, ] <- svd_block$u %*% diag(svd_block$d[1:rr], nrow = rr, ncol = rr) %*% solve(t(G_m))
    
    prev_n <- unique(c(prev_n, n_obs[[traverse_order[bb]]]))
  }
  
  # Do SVD
  post_svd <- svd(HH %*% t(VV), rr, rr)
  post_proc_spec <- diag(post_svd$d[1:rr], nrow = rr, ncol = rr) %*% t(post_svd$v)
  
  # Fit k-means
  sc_cq <- kmeans(t(post_proc_spec), clusts)
  return(sc_cq)
}
