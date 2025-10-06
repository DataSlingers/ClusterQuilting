
# Mosaic patch simulation data
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
test_output <- mos_create(n = base_n, p = base_p, v = base_v, b = base_b, h = base_h, clusts = base_clusts, rank_sim = base_rank, mean_cons = base_mean, sdd)


### Homogeneous patch ordering.
tp <- traverse_path_het(test_output$masked_dat, test_output$patch_obs, test_output$patch_features)
### Hetergeneous patch ordering.
# tp <- traverse_path_het(test_output$masked_dat, test_output$view_obs, test_output$view_features, 2)

### Fit CQ
cq_est <- ClusterQuilting(test_output$masked_dat, test_output$view_features,
                          test_output$view_obs, 2, 3, tp)


################################################################################
################################################################################
################################################################################

# Sequential patch simulation data
set.seed(44)
p_base <- 100
n_base <- 710
n_per_base <- 210
b_base <- 4
mean_cons_base <- 4.5
sd_cons_base <- 1
base_clusts <- 3
base_rank <- 2
test_output <- seq_create(n_base, p_base, n_per_base, b_base, base_clusts, base_rank,
                          mean_cons_base, sd_cons_base)

### Homogeneous patch ordering.
tp <- traverse_path_het(test_output$masked_dat, test_output$patch_obs, test_output$patch_features)
### Hetergeneous patch ordering.
# tp <- traverse_path_het(test_output$masked_dat, test_output$patch_obs, test_output$patch_features, 2)

### Fit CQ
cq_est <- ClusterQuilting(test_output$masked_dat, test_output$patch_features,
                          test_output$patch_obs, 2, 3, tp)