# Script for turning fluorescence trace df from 2 seperate scans into a single scan
# Leaves 
## Names are specific to the scans we use in the paper analysis.
library(tidyverse)

setwd("/") # CHANGE THIS to correct directory with existing data files
save_dir <- "./output_dir" # CHANGE THIS to desired output directory

# Load data; CHANGE NAMES TO YOUR OWN FILE NAMES
location <- read.csv("./location.csv", row.names=1)
ori <- read.csv("./orientation.csv", row.names=1)
block_s4s7 <- read.csv("./session4_scan7.csv")
block_s8s5 <- read.csv("./session8_scan5.csv")

# Standarize columns
block_s4s7 <- apply(block_s4s7, 2, scale)
block_s8s5 <- apply(block_s8s5, 2, scale)

# Filter to specific session/scans in our analysis
s4s7_ori_df <- ori %>% dplyr::filter((session == 4 & scan_idx == 7))
s8s5_ori_df <- ori %>% dplyr::filter((session == 8 & scan_idx == 5))

# Change uid to the same format as column names in fl trace data
s4s7_ori_df$uid <- paste0("u", s4s7_ori_df$unit_id)
s8s5_ori_df$uid <- paste0("u", s8s5_ori_df$unit_id)

# Keep columns that we have information for
existing_ids_s4s7 <- intersect(colnames(block_s4s7), s4s7_ori_df$uid)
existing_ids_s8s5 <- intersect(colnames(block_s8s5), s8s5_ori_df$uid)

# Find the union so we can align columns in full data.
s4s7_ids_num <- as.numeric(sub("^u", "", existing_ids_s4s7))
s4s7_target_ids <- s4s7_ori_df$target_id[match(existing_ids_s4s7, paste0("u", s4s7_ori_df$unit_id))]
s8s5_ids_num <- as.numeric(sub("^u", "", existing_ids_s8s5))
s8s5_target_ids <- s8s5_ori_df$target_id[match(existing_ids_s8s5, paste0("u", s8s5_ori_df$unit_id))]
full_target_ids <- unique(c(s4s7_target_ids, s8s5_target_ids))

# Create data matrix, observation information needed to run 
full_df <- matrix(0, nrow = nrow(block_s4s7) + nrow(block_s8s5),
                  ncol = length(full_target_ids))
colnames(full_df) <- full_target_ids

# Needed entry information for CQ algorithm.
patch_obs <- list()
patch_features <- list()
# obs_panel_list <- list()
# panel_obs_list <- list()

# First block
block1_cols <- match(full_target_ids, s4s7_target_ids)
full_df[c(1:nrow(block_s4s7)), which(!is.na(block1_cols))] <- 
  as.matrix(block_s4s7[, na.omit(block1_cols)])
patch_obs[[1]] <- which(!is.na(block1_cols))
patch_features[[1]] <- c(1:nrow(block_s4s7))

# Second block
block2_cols <- match(full_target_ids, s8s5_target_ids)
full_df[c((nrow(block_s4s7) + 1):(nrow(block_s4s7) + nrow(block_s8s5))), which(!is.na(block2_cols))] <- 
  as.matrix(block_s8s5[, na.omit(block2_cols)])
patch_obs[[2]] patch_features[[2]] <- which(!is.na(block2_cols))
patch_features[[2]] <- c((nrow(block_s4s7) + 1):(nrow(block_s4s7) + nrow(block_s8s5)))

# Create versions of matrices compatible with saving as csvs.
## Used for fitting comparison methods in Matlab.
cfun <- function(L) {
  pad.na <- function(x,len) {
    c(x,rep(NaN,len-length(x)))
  }
  maxlen <- max(sapply(L,length))
  do.call(data.frame,lapply(L,pad.na,len=maxlen))
}

# Save information necessary for downstream CQ analysis and comparison methods.
write.csv(t(full_df), "masked_dat_fl.csv", row.names = FALSE)

### Save observations for which each view is unmasked
panel_obs_mat <- matrix(0, n, 2)
for(bb in 1:2){
  panel_obs_mat[patch_obs[[bb]], bb] <- 1
}
write.csv(panel_obs_mat, "panel_obs_mat.csv", row.names = FALSE)

### Save features within each individual view
panel_times_mat <- matrix(0, p, 2)
for(bb in 1:2){
  panel_times_mat[patch_features[[bb]], bb] <- 1
}
write.csv(panel_times_mat, "panel_times_mat.csv", row.names = FALSE)

### Save observations within each patch
ss2 <- cfun(patch_obs)
colnames(ss2) <- NULL
write.csv(ss2, "sobs.csv", row.names = FALSE)

### Save features within each patch
ss4 <- cfun(patch_features)
colnames(ss4) <- NULL
write.csv(ss4, "stimes.csv", row.names = FALSE)

### Save non-zero entries
ss3 <- which(masked_dat != 0)
write.csv(ss3, "omega.csv", row.names = FALSE)

clust_ass <- rep(0, ncol(full_df)) # No true clusters, but need to throw something in to some of the comparions method functions
write.csv(clust_ass_vec, "clustass.csv", row.names = FALSE)
