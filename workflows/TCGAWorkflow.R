### Code for analyzing TCGA dataset: 
## Creating data frame from individual files, and metadata info.
## Getting method estimates.
## Creating results plots. 

setwd("~/ClusterQuilting") # Change to point to top of code directory

################################################################################

# Instructions for downloading individual files can be found in the TCGA Readme.

################################################################################

# Create observation information; runs as a script
source("./TCGA/fuse_data.R")
setwd("~/ClusterQuilting") 

################################################################################

# Functions for creating full data frame, and simulated mosaic patch missingness
source("./TCGA/TCGA_datamake.R")

# Create synthetic masked data.
data_dir <- "./output_dir/" #Change to directory where data is saved locally
save_dir <- "./results_dir/" #Change to directory where new data csv should be saved
settings_dir <- "./settings_dir/" #Change to directory where settings file is saved locally
cancer_folds <- c("Lung", "Breast", "Kidney")
data_view_folds <- c("RPPA", "Transcriptome", "Gene", "Methylation", "CNV")
view_df <- read.csv(paste0(data_dir, "/view_miss.csv"), row.names = 1) # File created in fuse_data.R.
full_sample_sheet <- read.csv(paste0(settings_dir, "/full_sample_sheet.csv"), row.names = 1) # File created in fuse_data.R.


base_n <- 10000 # For testing, can set maximum observations to small number to reduce runtime.
base_p <- rep(100000, 8) # For testing, can set maximum features per modality to small number to reduce runtime.
base_v <- 5
base_b <- c(2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6)
base_h <- c(3, 4, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4)
iters <- 50 # Use to run multiple simulations at once.

save_list <- list()
for(ii in 1:length(base_b)){
  ttt <- mvp_create_tcga(view_df = view_df, data_dir = data_dir, cancer_folds = cancer_folds,
                         settings_dir = settings_dir, full_sample_sheet = full_sample_sheet,
                         n = base_n, p_max = base_p, v = base_v, b = base_b[ii], h = base_h[ii],
                         iters = iters, wd = save_dir)
  save_list[[ii]] <- list()
  save_list[[ii]][["info"]] <- save_info
  save_list[[ii]][["blocks"]] <- base_b[ii]
  save_list[[ii]][["view_per_blocks"]] <- base_h[ii]
}
save(save_list, file = paste0(save_dir, "info.Rdata"))

setwd("~/ClusterQuilting") 

################################################################################

### Run comprunner.m in Matlab with corresponding save_dir pointer
### in order to get estimates of clusters from comparison methods used in paper.

# Get CQ estimate, compile evaluation metric results for each simulation setting.
source("./CompileARI.R")

## Oracle tuning
compile_ari_oracle(sim_dir, rank_solve = 2, clust_solve = 3,
                  homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

## Data driven tuning
ranks <- c(2, 3, 4, 5, 6)
clusts <- c(2, 3, 4, 5, 6)
compile_ari_dd(save_dir, ranks, clusts, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

################################################################################

source("./MakeResultsPlots.R")
# Make ARI results figure (Figure 3)
folds <- c("clusts_BrainBreastKidney_k_2_views_5_vpb_3",
           "clusts_BrainBreastKidney_k_2_views_5_vpb_4",
           "clusts_BrainBreastKidney_k_3_views_5_vpb_3",
           "clusts_BrainBreastKidney_k_3_views_5_vpb_4",
           "clusts_BrainBreastKidney_k_4_views_5_vpb_2",
           "clusts_BrainBreastKidney_k_4_views_5_vpb_3",
           "clusts_BrainBreastKidney_k_4_views_5_vpb_4",
           "clusts_BrainBreastKidney_k_5_views_5_vpb_2",
           "clusts_BrainBreastKidney_k_5_views_5_vpb_3",
           "clusts_BrainBreastKidney_k_5_views_5_vpb_4",
           "clusts_BrainBreastKidney_k_6_views_5_vpb_2",
           "clusts_BrainBreastKidney_k_6_views_5_vpb_3",
           "clusts_BrainBreastKidney_k_6_views_5_vpb_4")

labels <- paste0(paste0(base_b, " Blocks"), ", ",  paste0(base_h, " Views"))
plot_func(paste0(save_dir, folds), labels,
          "Number of Blocks, Views Per Block", "tcga.png")

setwd("~/ClusterQuilting") 

################################################################################

# Make heatmap (Figure 4)
source("./TCGA/TCGA_heatmap.R")

produce_tcga_heatmap(paste0(save_dir, "clusts_BrainBreastKidney_k_2_views_5_vpb_4/", 
                            "clusts_BrainBreastKidney_k_2_views_5_vpb_4_1/"),
                     "./new_dir")

setwd("~/ClusterQuilting") 