### Code for analyzing simulation dataset: 
## Creating simulated GMM data.
## Getting method estimates.
## Creating results plots. 

setwd("~/ClusterQuilting")  # Change to point to top of repository directory

################################################################################
################################################################################
################################################################################

# Mosaic patch GMM simulations; example Fig 2
source("./MosaicPatchSim.R")
source("./HeterogeneousSNRPatchOrdering.R")
source("./HomogenousSNRPatchOrdering.R")
source("./ClusterQuilting.R")
source("./CompileARI.R")
source("./MakeResultsPlots.R")

################################################################################

## Set simulation parameters
### What these parameters are can be found in MosaicPatchSim.R
base_n <- 840
base_v <- 12
base_p <- rep(50, base_v)
base_b <- 4
base_h <- 6
base_clusts <- 3
base_rank <- 2
base_mean <- 4.5
sdd <- 1
iters <- 50
save_dir <- "~/IWantToSaveHere/" # Change this to desired save location

################################################################################

### Create simulation data
# Changing views per block
test_h <- c(4, 5, 6, 7, 8)
save_list <- list()
for(ii in test_h){
  save_info <- mos_create(n = base_n, p = base_p, v = base_v, b = base_b, h = ii, clusts = base_clusts,
                    rank_sim = base_rank, mean_cons = base_mean, sd_cons = sdd, iters = iters, 
                    wd = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "h.Rdata"))

# Changing number of blocks
test_b <- c(3, 4, 5, 6, 7)
save_list <- list()
for(ii in test_b){
  save_info <- mos_create(n = base_n, p = base_p, v = base_v, b = ii, h = base_h, clusts = base_clusts,
                    rank_sim = base_rank, mean_cons = base_mean, sd_cons = sdd, iters = iters,
                    wd = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "b.Rdata"))

# Changing number of underlying clusters
test_clust <- c(2, 3, 4, 5, 6)
save_list <- list()
for(ii in test_clust){
  save_info <- mos_create(n = base_n, p = base_p, v = base_v, b = base_b, h = base_h, clusts = ii,
                    rank_sim = base_rank, mean_cons = base_mean, sd_cons = sdd, iters = iters, 
                    wd = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "clust.Rdata"))

# Changing cluster distance
test_mean <- c(0.5, 1.5, 2.5, 3.5, 4.5)
save_list <- list()
for(ii in test_mean){
  save_info <- mvp_create(n = base_n, p = base_p, v = base_v, b = base_b, h = base_h, clusts = base_clusts,
                    rank_sim = base_rank, mean_cons = ii, sd_cons = sdd, iters = iters, 
                    wd = save_dir)
}
save(save_list, file = paste0(save_dir, "dist.Rdata"))

setwd("~/ClusterQuilting") 

################################################################################

### Run comprunner.m in Matlab with corresponding save_dir pointer
### in order to get estimates of clusters from comparison methods used in paper.

# Get CQ estimate, compile evaluation metric results for each simulation setting.
## For Supp. Figure 5 with homogeneous SNR: change homo_SNR = TRUE
compile_ari_oracle(save_dir, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

## For data driven tuning results:
# ranks <- c(2, 3, 4, 5, 6)
# clusts <- c(2, 3, 4, 5, 6)
# 
# compile_ari_dd(save_dir, ranks, clusts, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

################################################################################

### Make subfigures
# Changing views per block
folds <- c("clusts_3_n_840_p_600_k_4_views_12_vpb_4_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_5_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_7_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_8_means_4.5_rank_2_off_0")
plot_func(paste0(save_dir, folds), c("4", "5", "6", "7", "8"),
          "Views Per Block", "vpb.png")

# Changing number of blocks
folds <- c("clusts_3_n_840_p_600_k_3_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_5_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_6_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_7_views_12_vpb_6_means_4.5_rank_2_off_0")
plot_func(paste0(save_dir, folds), c("3", "4", "5", "6", "7"), 
          "Number of Blocks", "numblocks.png")

# Changing number of underlying clusters
folds <- c("clusts_2_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_4_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_5_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_6_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0")
plot_func(paste0(save_dir, folds), c("2", "3", "4", "5", "6"),
          "Number of Clusters", "clusts.png")

# Changing cluster distance
folds <- c("clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_0.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_1.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_2.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_3.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0")
plot_func(paste0(save_dir, folds),  c("0.5", "1.5", "2.5", "3.5", "4.5"),
          "Cluster Distance", "dist.png")

setwd("~/ClusterQuilting") 

################################################################################
################################################################################
################################################################################

# Sequential patch GMM simulations; example Supp Fig 3
p_base <- 100
n_base <- 710
n_per_base <- 210
num_blocks_base <- 4
mean_cons_base <- 4.5
sd_cons_base <- 1
iters <- 50
save_dir <- "~/IWantToSaveHere/" # Change this to desired save location

################################################################################

### Create simulation data
# Changing block sizes
n_per_s <- c(180, 190, 200, 210, 220)
save_list <- list()
for(ii in n_per_s){
  save_info <- seq_create(n = n_base, p = p_base, n_per = ii, b = b_base, 
                          clusts = base_clusts, rank_sim = base_rank,
                          mean_cons = mean_cons_base, sd_cons = sd_cons_base,
                          iters = iters, save_dir = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "n_per.Rdata"))

# Changing number of features
test_p <- c(25, 50, 75, 100, 125)
save_list <- list()
for(ii in test_p){
  save_info <- seq_create(n = n_base, p = ii, n_per = n_per_base, b = b_base, 
                          clusts = base_clusts, rank_sim = base_rank,
                          mean_cons = mean_cons_base, sd_cons = sd_cons_base,
                          iters = iters, save_dir = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "p.Rdata"))

# Changing number of clusters
num_clusts_s <- c(2, 3, 4, 5, 6)
save_list <- list()
for(ii in num_clusts_s){
  save_info <- seq_create(n = n_base, p = test_p, n_per = n_per_base, b = b_base, 
                          clusts = base_clusts, rank_sim = base_rank,
                          mean_cons = ii, sd_cons = sd_cons_base,
                          iters = iters, save_dir = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "clusts.Rdata"))

# Changing distance between clusters
means_cons_s <- c(2, 3, 4, 5, 6) + 0.5
save_list <- list()
for(ii in means_cons_s){
  save_info <- seq_create(n = n_base, p = test_p, n_per = n_per_base, b = b_base, 
                          clusts = ii, rank_sim = base_rank,
                          mean_cons = mean_cons_base, sd_cons = sd_cons_base,
                          iters = iters, save_dir = save_dir)
  save_list[[as.character(ii)]] <- save_info
}
save(save_list, file = paste0(save_dir, "dists.Rdata"))

setwd("~/ClusterQuilting") 

################################################################################

### Run comprunner.m in Matlab with corresponding save_dir pointer
### in order to get estimates of clusters from comparison methods used in paper.

# Get CQ estimate, compile evaluation metric results for each simulation setting.
## For Supp. Figure 4 with homogeneous SNR: change homo_SNR = TRUE
compile_ari_oracle(save_dir, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

## For data driven tuning results:
# ranks <- c(2, 3, 4, 5, 6)
# clusts <- c(2, 3, 4, 5, 6)
# 
# compile_ari_dd(save_dir, ranks, clusts, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

################################################################################

### Make subfigures
# Changing block sizes
folds <- c("clusts_3_n_710_p_100_o_180_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_190_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_200_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_220_k_4_means_4.5_sd_1")
plot_func(paste0(save_dir, folds), c("180", "190", "200", "210", "220"), 
          "Blocks Size", "blocksize.png")

# Changing number of blocks
folds <- c("clusts_3_n_710_p_100_o_280_k_3_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_168_k_5_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_140_k_6_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_120_k_7_means_4.5_sd_1")
plot_func(paste0(save_dir, folds), c("3", "4", "5", "6", "7"), 
          "Number of Blocks", "numblocks.png")

# Changing number of underlying clusters
folds <- c("clusts_2_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_3_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_4_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_5_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0",
           "clusts_6_n_840_p_600_k_4_views_12_vpb_6_means_4.5_rank_2_off_0")
plot_func(paste0(save_dir, folds), c("2", "3", "4", "5", "6"),
          "Number of Clusters", "clusts.png")

# Changing cluster distance
folds <- c("clusts_3_n_710_p_100_o_210_k_4_means_2.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_3.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_4.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_5.5_sd_1",
           "clusts_3_n_710_p_100_o_210_k_4_means_6.5_sd_1")
plot_func(paste0(save_dir, folds),  c("2.5", "3.5", "4.5", "5.5", "6.5"),
          "Cluster Distance", "dists.png")

setwd("~/ClusterQuilting") 
