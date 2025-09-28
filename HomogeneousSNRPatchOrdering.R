#### Find traverse order
### Inputs:
## raw_dat: raw data matrix with features as rows and observations as columns.
## pol: list of observations for each patch to merge.
## bpl: list of features for each patch to merge.
traverse_path <- function(raw_dat, pol, bpl){
  
  # Init
  tp_vec <- c()
  MM <- sample(c(1:length(pol)))
  
  dat <- t(raw_dat)
    
  # First iter
  score_mat <- matrix(0, length(MM), length(MM))
  for(m1 in 1:length(MM)){
    score_denom <- sum(abs(dat[pol[[MM[m1]]], bpl[[MM[m1]]]]))
    score_num <- 0
    for(m2 in 1:length(MM)){
      if(m1 == m2) {next}
      overlap_rows <- intersect(pol[[MM[m1]]], pol[[MM[m2]]])
      if(length(overlap_rows) == 0) {next}
      score_num <- sum(abs(dat[unlist(overlap_rows), bpl[[MM[m1]]]]))
      score_mat[m1, m2] <- score_num / score_denom
    }
  }
  kk <- which(score_mat == max(score_mat), arr.ind = TRUE)[1, ]
  tp_vec <- c(MM[kk[1]], MM[kk[2]])
  MM <- MM[-which(MM %in% tp_vec)]
    
  # Future iter
  while(length(MM) > 1){
    score_vec <- rep(0, length(MM))
    current_group_row <- unique(unlist(pol[tp_vec]))
    for(m1 in 1:length(MM)){
      overlap_rows <- intersect(pol[[MM[m1]]], current_group_row)
      if(length(overlap_rows) == 0) {next}
      score_denom <- sum(abs(dat[pol[[MM[m1]]], bpl[[MM[m1]]]]))
      score_num <- sum(abs(dat[overlap_rows, bpl[[MM[m1]]]]))
      score_vec[m1] <- score_num / score_denom
    }
    kk <- which(score_vec == max(score_vec))[1]
    tp_vec <- c(tp_vec, MM[kk])
    MM <- MM[-kk]
  }
  tp_vec <- c(tp_vec, MM[1])

  
  return(tp_vec)
}

