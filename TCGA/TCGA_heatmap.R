# This script includes functions to produce the heatmap of Figure 5 in Zheng, Chang, and Allen 2025. 
library(tidyverse)

#### Cluster Quilting function that only returns imputed matrix from first step.
## dat: raw data matrix with features as columns and observations as rows.
## p_obs: list of vectors of features for each patch.
## n_obs: list of vectors of observations for each patch.
## rr: ranks
## clusts: number of clusters.
## traverse_order: order in which patches should be merged.
ClusterQuilting_est <- function(data_matrix, p_obs, n_obs, rr, clusts, traverse_order){
  
  
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
  
  hv <- HH %*% t(VV)
  return(hv)
}

########################################################################

## Produces a heatmap of TCGA data.
## Inputs:
# data_dir: location of folder containing a run of masked data.
# save_dir: directory where plot should be saved.
# homo_SNR: homogeneous SNR assumption for patch ordering?
produce_tcga_heatmap <- function(data_dir, save_dir, homo_SNR = FALSE){
  data_mat <- read.csv("./full_dat.csv")
  masked_dat <- read.csv("./masked_dat.csv")
  plot_dat <- masked_dat
  
  # Truncate for plotting
  plot_dat[which(masked_dat  > 3.5, arr.ind = TRUE)] <- 3.5
  plot_dat[which(masked_dat  < -3.5, arr.ind = TRUE)] <- -3.5
  
  load("./trueinfo.RData")
  
  if(homo_SNR){
    tpp <- traverse_path(masked_dat, panel_obs_list, blocks_p_list)
  } else {
    tpp <- traverse_path_het(masked_dat, panel_obs_list, blocks_p_list, rank_est)
  }
  cq_val <- ClusterQuilting_est(t(masked_dat), blocks_p_list, panel_obs_list, 
                            rank_sim, clusts, tpp)
  fff2 <- t(cq_val)
  cq_est <- ClusterQuilting(t(masked_dat), blocks_p_list, panel_obs_list, 
                                rank_sim, clusts, tpp)
  
  row_order <- order(cq_est)
  cbp <- c("#004D40", "#1E88E5", "#D81B60")

  # Make true cluster bar (image editor to place next to each heatmap)
  real_ass <- data.frame(assignment = factor(rep(clust_ass_vec[row_order], each = 2), 
                                             labels = c("Kidney", "Breast", "Lung")),
                         x = rep(1:nrow(masked_dat), each = 2),
                         y = rep(-1:0, times = nrow(masked_dat)))
  
  p <- ggplot() +
    geom_tile(aes(x = x, y = y, fill = assignment), color = "white", data = real_ass, show.legend = FALSE) +
    scale_fill_manual(values = cbp) +
    theme_void() +
    coord_flip() +
    labs(fill = "True\nLabels") +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          legend.title = element_text(size = 32),
          legend.text = element_text(size = 24))
  
  ggsave(paste0(save_dir, "cq_ass.png"), p, dpi = 300, height = 7, width = 7, units = "in")
  
  
  # Create heatmaps
  
  for(ii in 1:length(blocks_p_list)){
    dd <- dist(t(fff2[row_order, blocks_p_list[[ii]]]))
    hh <- hclust(dd)
    plot_df1 <- data.frame(Observations = rep(1:nrow(masked_dat), times = rep(length(blocks_p_list[[ii]]))),
                           Features = rep(1:length(blocks_p_list[[ii]]), each = nrow(masked_dat)),
                           Value = c(as.matrix(fff2[row_order, blocks_p_list[[ii]]][, hh$order])))
    
    p <-ggplot() +
      geom_tile(aes(x = Observations, y = Features, fill = Value), color = "grey",
                data = plot_df1, show.legend = FALSE) +
      scale_fill_gradient2(high = "#fde725", low = "#440154", mid = "white") +
      theme_void() +
      coord_flip() + 
      theme(plot.margin = margin(0,0,0,0),
            plot.title = element_text(size = 32, margin=margin(5,0,-15,0)))
    
    ggsave(paste0(save_dir, "heatmap_block", ii, ".png"), p, dpi = 300)
  }
    
}
