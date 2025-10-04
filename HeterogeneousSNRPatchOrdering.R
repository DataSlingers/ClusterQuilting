traverse_path_het <- function(dat, pol, bpl, rr){
  
  # Init
  tp_vec <- c()
  MM <- sample(c(1:length(pol)))
  dat <- t(dat)
  sig_hat_patches <- rep(0, length(pol))
  snr_patches <- rep(0, length(pol))

  for(ii in 1:length(pol)){
    ## should set bet = smaller dimension/larger dimension; svr = 1/bet
    bet <- max(length(pol[[ii]]), length(bpl[[ii]])) / min(length(pol[[ii]]), length(bpl[[ii]]))
    median_mp <- qmp(0.5, svr = 1/bet)
    ## Each patch's noise level should be normalized based on this patch's larger dimension (should be bpl[ii]=300 here)
    sig_hat_patches[ii] <- median(svd(dat[pol[[ii]], bpl[[ii]]])$d) / 
      sqrt(min(length(pol[[ii]]), length(bpl[[ii]])) * median_mp)
  
    ## ZZ should be normalized based on estimated sigma_hat
    ZZ <- dat[pol[[ii]], bpl[[ii]]] / sqrt(min(length(pol[[ii]]), length(bpl[[ii]]))) / sig_hat_patches[ii]
    sig_r <- svd(ZZ)$d[rr]
    ## SNR formula
    snr_patches[ii] <- pmax(sqrt(sig_r^2 - bet - 1 + sqrt(pmax((sig_r^2 - bet - 1) ^ 2 - 4 * bet, 0))), 10^-9) / 
      (sqrt(2) * (1 + sqrt(bet)))
  }
  # First iter
  score_mat <- matrix(NA, length(MM), length(MM))
  for(m1 in 1:length(MM)){
    for(m2 in 1:length(MM)){
      if(m1 == m2) {next}
      overlap_rows <- intersect(pol[[MM[m1]]], pol[[MM[m2]]])
      if(length(overlap_rows) == 0) {next}
      sig_r <- svd(dat[unlist(overlap_rows), bpl[[MM[m2]]]])$d[rr]
      score_mat[m1, m2] <- (1.1 * norm(dat[pol[[MM[m2]]], bpl[[MM[m2]]]], "2") / 
                              (sig_r) + 1) * 
       (1 / snr_patches[m1] + 1 / snr_patches[m2])
    }
  }
  kk <- which(score_mat == min(score_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  tp_vec <- c(MM[kk[1]], MM[kk[2]])
  curr_score <- min(score_mat, na.rm = TRUE)
  MM <- MM[-which(MM %in% tp_vec)]

  # Future iter
  while(length(MM) > 0){
    score_vec <- rep(NA, length(MM))
    current_group_row <- unique(unlist(pol[tp_vec]))
    for(m1 in 1:length(MM)){
      overlap_rows <- intersect(pol[[MM[m1]]], current_group_row)
      if(length(overlap_rows) == 0) {next}
      sig_r <- svd(dat[unlist(overlap_rows), bpl[[MM[m1]]]])$d[rr]
      score_vec[m1] <- (1.1 * norm(dat[pol[[MM[m1]]], bpl[[MM[m1]]]], "2") / 
                         (sig_r) + 1) * 
       (1 / snr_patches[m1] + curr_score)
    }
    kk <- which(score_vec == min(score_vec, na.rm = TRUE))[1]
    
    curr_score <- min(score_vec, na.rm = TRUE)
    tp_vec <- c(tp_vec, MM[kk])
    MM <- MM[-kk]
  }
}

# tp <- traverse_path_het(test_output$masked_dat, test_output$view_obs, test_output$view_features, 2)