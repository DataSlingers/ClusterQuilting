library(tidyverse)
## Makes figures 5 and 6 from Zheng, Chang, and Allen (2025).
## Inputs:
# data_dir: Overall directory of where data and info are currently stored.
# Should contain .Rdata created by the make_full.R script in this directory.
# est_dir: Directory of where final selected clustering results for different methods are stored.
# Should contain the .csv files created by running compile_ari_microns.
# save_dir: directory where plots should be saved.
figures_microns <- function(data_dir, est_dir, save_dir){
  
  ########################################################################
  
  ### Load the necessary information
  # Load data
  masked_dat <- read.csv(paste0(data_dir, "/masked_dat.csv"), row.names = NULL)
  # Load the CQ pipeline info created in R
  load(paste0(data_dir, "./trueinfo.RData")) # Reload of info 
  
  # Load clustering results
  CQ_res <- read.table(paste0(est_dir,"/CQ_res.csv"), quote="\"", 
                       comment.char="", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  IMG_res <- read.table(paste0(data_dir,"/IMG_res.csv"), quote="\"", comment.char="")
  IMSCAGL_res <- read.table(paste0(data_dir,"/IMSCAGL_res.csv"), quote="\"", comment.char="")
  DAIMC_res <- read.table(paste0(data_dir,"/DAIMC_res.csv"), quote="\"", comment.char="")
  OPIMC_res <- read.table(paste0(data_dir,"/OPIMC_res.csv"), quote="\"", comment.char="")
  NN_res <- read.table(paste0(data_dir,"/NN_res.csv"), quote="\"", comment.char="")
  
  # Process the location data into plottable values.
  location <- read.csv("./location.csv", row.names=1)
  location$x <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(location$pt_position, "[ ]+", 
                                                                           fixed = FALSE), "[[", 1)),"\\["), "[[", 2)))
  location$y <- as.numeric(unlist(lapply(strsplit(location$pt_position, "[ ]+", 
                                                    fixed = FALSE), "[[", 2)))
  location$z <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(location$pt_position, "[ ]+", 
                                                                           fixed = FALSE), "[[", 3)),"\\]"), "[[", 1)))
  
  # Compile clustering results into a wide format data frame.
  results_df <- as.data.frame(CQ_res)
  colnames(results_df) <- "cq"
  results_df$target_id <- unique(location$target_id)
  results_df$nn <- test_sc_kmeans
  results_df$imscagl <- IMSCAGL_res$V1
  results_df$img <- IMG_res$V1
  results_df$daimc <- DAIMC_res$V1
  results_df$opimc <- OPIMC_res$V1
  
  # Combine clustering results with location data
  location <- left_join(location, results_df, by = join_by(target_id))
  location$cq <- factor(location$cq)
  location$nn <- factor(location$nn)
  location$imscagl <- factor(location$imscagl)
  location$daimc <- factor(location$daimc)
  location$opimc <- factor(location$opimc)
  location$img <- factor(location$img)
  
  # Only keep neurons that have estimates (non-missing filter)
  keeprows <- c(1)
  for(ii in 2:nrow(location)){
    if(!(location$target_id[ii] %in% location$target_id[1:(ii-1)])){
      keeprows <- c(keeprows, ii)
    }
  }
  location <- location[keeprows,]
  
  # Calculate (approximate) boundaries for scans
  patch_coord <- location %>% 
    group_by(session, scan_idx) %>% 
    summarise(maxx = max(x), minx = min(x),
              maxy = max(y), miny = min(y),
              maxz = max(z), minz = min(z),
              qmaxz = quantile(z, 0.9), qminz = quantile(z, 0.1))
  
  patch_coord$session <- factor(patch_coord$session)
  
  ########################################################################
  
  ### Code for Figure 5
  # Neurons plotted at spatial location, colored by clustering assignment.
  # (Approximate) patches boundaries colored in background 
  # Show example for CQ cluster estimates
  # Can plot for different methods by changing the color argument in geom_point
  # (e.g. to make corresponding plots in Supplemental Figure 12.)
  p <- ggplot(data = location) +
    geom_rect(aes(xmax = maxx, xmin = minx, ymax = maxy, ymin = miny, fill = session), data = patch_coord, 
              alpha = 0.2, show.legend = F) +
    scale_fill_manual(values = c("brown", "yellow")) +
    geom_point(aes(x = x, y = y, color = cq),
               alpha = 0.4) +
    labs(x = "", y = "", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  p
  ggsave(paste0(save_dir, "/CQ_xy.pdf"))
  
  p <- ggplot(data = location) +
    geom_rect(aes(xmax = maxx, xmin = minx, ymax = maxz, ymin = minz, fill = session), data = patch_coord, 
              alpha = 0.2, show.legend = F) +
    scale_fill_manual(values = c("brown", "yellow")) +
    geom_point(aes(x = x, y = z, color = cq),
               alpha = 0.4) +
    labs(x = "", y = "", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  p
  ggsave(paste0(save_dir, "/CQ_xz.pdf", p))
  
  p <- ggplot(data = location) +
    geom_rect(aes(xmax = maxy, xmin = miny, ymax = maxz, ymin = minz, fill = session), data = patch_coord, 
              alpha = 0.2, show.legend = F) +
    geom_point(aes(x = y, y = z, color = cq),
               alpha = 0.4) +
    scale_fill_manual(values = c("brown", "yellow")) +
    labs(x = "", y = "", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  p
  ggsave(paste0(save_dir, "/CQ_yz.pdf", p))
  
  
  ########################################################################
  
  # Also make average trace plots.
  tt11 <- colMeans(masked_dat[intersect(which(location$cq == 1), blocks_n_list[[1]]),
                                 blocks_p_list[[1]]])
  
  p11 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[1]])), y = tt11)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust1block1.pdf"), p11,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  
  tt12 <- colMeans(masked_dat[intersect(which(location$cq == 2), blocks_n_list[[1]]),
                                 blocks_p_list[[1]]])
  
  p12 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[1]])), y = tt12)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust2block1.pdf"), p12,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  tt13 <- colMeans(masked_dat[intersect(which(location$cq == 3), blocks_n_list[[1]]),
                                 blocks_p_list[[1]]])
  
  p13 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[1]])), y = tt13)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust3block1.pdf"), p13,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  tt21 <- colMeans(masked_dat[intersect(which(location$cq == 1), blocks_n_list[[2]]),
                                 blocks_p_list[[2]]])
  
  p21 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[2]])), y = tt21)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust1block2.pdf"), p21,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  tt22 <- colMeans(masked_dat[intersect(which(location$cq == 2), blocks_n_list[[2]]),
                                 blocks_p_list[[2]]])
  
  p22 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[2]])), y = tt22)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust2block2.pdf"), p22,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  tt23 <- colMeans(masked_dat[intersect(which(location$cq == 3), blocks_n_list[[2]]),
                                 blocks_p_list[[2]]])
  
  p23 <- ggplot() +
    geom_path(aes(x = c(1:length(blocks_p_list[[2]])), y = tt23)) + 
    labs(x = "Time", y = "Average Trace", color = "Cluster") +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  
  ggsave(paste0(save_dir, "/clust3block2.pdf"), p23,
         limitsize = FALSE, width = 6, height = 3, units = "in")
  
  ########################################################################
  
}