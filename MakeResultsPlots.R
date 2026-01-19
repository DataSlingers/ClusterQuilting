library(kableExtra)
library(tidyverse)

## Make line plots that appear in the paper.
## These show how performance of methods (by ARI) change across changing sim parameter value.
## Inputs:
# sim_par_folds: the subfolders of simulation results for different parameter settings;
# this should contain the ari_results.csv file produced by running compile_ari.
# x_params: the x-axis tick labels; used for denoting changing values of a simulation parameter;
# practically, these should correspond to sim_par_folds.
# x_lab: The x-axis label - what simulation parameter is changing along the x-axis?
# save_name: Name of file that plot should be saved to.
plot_func <- function(sim_par_folds, x_params, x_lab,
                      save_name){
  
  cp <- c("black", "#0016ff", "#f5832c", "#954000", "#24f123", "#03ba9b", "#ff0000")
  mod_names <- c("Full Data", "IMSC_AGL", "IMG", "DAIMC", "OPIMC", "NN Data", "ClustQuilt")
  
  mean_tab <- matrix(NA, length(mod_names), 0)
  sd_tab <- matrix(NA, length(mod_names), 0)
  for(spf in 1:length(sim_par_folds)){
    setwd(sim_par_folds[spf])
    # Load results
    fff <- read.csv("ari_results.csv")
    mean_tab <- cbind(mean_tab, fff$ari)
    sd_tab <- cbind(sd_tab, fff$sd)
  }
  
  # Create data frame for ggplot plotting
  plot_df <- data.frame(adjpar = x_params, each = length(mod_names),
                        model = factor(rep(mod_names, times = length(sim_par_folds)),
                                       levels = mod_names))
  
  plot_df$tp <- c(mean_tab)
  plot_df$tp_sd <- c(sd_tab) 
  
  p <- ggplot() +
    geom_path(aes(factor(adjpar), y = tp, group = model, color = model), data = plot_df, alpha = 0.5) +
    geom_point(aes(factor(adjpar), y = tp, color = model), data = plot_df, alpha = 0.5, size = 2.5) +
    geom_errorbar(aes(factor(adjpar), ymax = pmin(tp + tp_sd, 1), ymin = pmax(tp - tp_sd, 0), 
                      group = model, color = model), data = plot_df, width = 0.2, alpha = 0.35) +
    labs(x = x_lab, y = "Adjusted Rand Index", color = "Model") +
    scale_y_continuous(limits = c(-0.1, 1)) + 
    scale_color_manual(values = cp) + 
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 10)
    )
  
  ggsave(save_name, p, dpi = 300)
}
