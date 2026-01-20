### Code for analyzing MICrONs dataset.
### Also makes Figures 5 and 6 in Zheng, Chang, and Allen (2025).
## Creating data frame of functional data from database
## Getting method estimates.
## Creating results plots. 

setwd("~/ClusterQuilting")  # Change to point to top of code directory

################################################################################

# Instructions for downloading database can be found in the MICrONs Readme.
# Then, run trace_extract.py in a Python env to get fluorescence traces.

################################################################################

## Put data in a form that we can use in our current CQ pipeline; runs as a script.
source("./MICRONS/make_full.R")

setwd("~/ClusterQuilting") 

################################################################################

## Data driven tuning of best hyperparameters for each method;
## Then save corresponding cluster estimates
source("./MICRONS/CompileARIMicrons.R")
ranks <- c(2, 3, 4, 5, 6)
clusts <- c(2, 3, 4, 5, 6)
data_dir <- "~/MICRONS" # Change this to where MICRONS data is saved.
est_dir <- "/MICRONS_est" # Change this to where estimates for comparison methods are saved
compile_ari_microns(data_dir, est_dir, ranks, clusts, homo_SNR = FALSE)

setwd("~/ClusterQuilting") 

################################################################################

## Create figures 
source("./MICRONS/MICRONSfigures.R")
data_dir <- "~/MICRONS" # Change this to where MICRONS data is saved.
est_dir <- "/MICRONS_est" # Change this to where estimates for comparison methods are saved
save_dir <- "/MICRONS_plot" # Change this to where plots should be saved

figures_microns(data_dir, est_dir, save_dir)

setwd("~/ClusterQuilting") 