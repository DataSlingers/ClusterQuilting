setwd("../") # Make sure this points to top of repository directory

################################################################################

# Mosaic patch simulation data
source("./MosaicPatchSim.R")
source("./HeterogeneousSNRPatchOrdering.R")
source("./HomogenousSNRPatchOrdering.R")
source("./ClusterQuilting.R")
source("./CompileARI.R")

## Set simulation parameters
### What these parameters are can be found in MosaicPatchSim.R
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
iters <- 50
save_dir <- "~/IWantToSaveHere/" # Change this to desired save location

## Create simulated data + mosaic patch missingness
test_output <- mos_create(n = base_n, p = base_p, v = base_v, b = base_b, 
                          h = base_h, clusts = base_clusts, rank_sim = base_rank, 
                          mean_cons = base_mean, sd_cons = sdd, iters = iters, 
                          save_dir = save_dir)

### Run comprunner.m in Matlab with corresponding save_dir pointer
### in order to get estimates of clusters from comparison methods used in paper.

## Create patch ordering
### Homogeneous SNR patch ordering.
tp <- traverse_path(test_output$masked_dat, test_output$view_obs, test_output$view_features)
### Hetergeneous SNR patch ordering. Requires a rank.
# tp <- traverse_path_het(test_output$masked_dat, test_output$view_obs, test_output$view_features, 2)

### Fit CQ, under assumption of rank 2 with 3 clusters. Save results.
cq_est <- ClusterQuilting(test_output$masked_dat, test_output$view_features,
                          test_output$view_obs, 2, 3, tp)
write.csv(cq_est, paste0(save_dir, "CQ_res.csv"), row.names = FALSE)


################################################################################
################################################################################
################################################################################

# Sequential patch simulation data

## Load scripts
source("./SeqPatchSim.R")
source("./HeterogeneousSNRPatchOrdering.R")
source("./HomogenousSNRPatchOrdering.R")
source("./ClusterQuilting.R")
source("./CompileARI.R")

## Set simulation parameters
### What these parameters are can be found in SeqPatchSim.R
set.seed(44)
p_base <- 100
n_base <- 710
n_per_base <- 210
b_base <- 4
mean_cons_base <- 4.5
sd_cons_base <- 1
base_clusts <- 3
base_rank <- 2
iters <- 50
save_dir <- "~/IWantToSaveHere/" # Change this to desired save location

## Create simulated data + sequential patch missingness
test_output <- seq_create(n = n_base, p = p_base, n_per = n_per_base, b = b_base, 
                          clusts = base_clusts, rank_sim = base_rank,
                          mean_cons = mean_cons_base, sd_cons = sd_cons_base,
                          iters = iters, save_dir = save_dir)

### Run comprunner.m in Matlab with corresponding save_dir pointer
### in order to get estimates of clusters from comparison methods used in paper.

## Create patch ordering
### Homogeneous SNR patch ordering.
tp <- traverse_path(test_output$masked_dat, test_output$view_obs, test_output$view_features)
### Hetergeneous SNR patch ordering. Requires a rank.
# tp <- traverse_path_het(test_output$masked_dat, test_output$view_obs, test_output$view_features, 2)

## Fit CQ, under assumption of rank 2 with 3 clusters
cq_est <- ClusterQuilting(test_output$masked_dat, test_output$patch_features,
                          test_output$patch_obs, 2, 3, tp)
write.csv(cq_est, paste0(save_dir, "CQ_res.csv"), row.names = FALSE)
