
#### Normalize matrix 
### Inputs:
## dat: dataset as a matrix.
row_normalize <- function(dat) {
  dat <- dat - rowMeans(dat)
  nrm <- apply(dat, 1, sd)
  nrm[nrm == 0] <- 1
  return(dat / nrm)
}

#### Build similarity matrix
## dat: raw data matrix with features as columns and observations as rows.
similarity_mat <- function(dat) {
  
  # Pairwise squared Euclidean distances
  d_sq <- as.matrix(dist(dat))^2
  
  sigma_vals <- d_sq[upper.tri(d_sq)]
  sigma <- sqrt(median(sigma_vals[sigma_vals > 0], na.rm = TRUE))
  
  # RBF similarity
  W <- exp(-d_sq / (d_sq * sigma^2))
  diag(W) <- 0
  
  return(W)
}

#### Build normalized graph lapacian
### Inputs:
## sim_mat: similarity matrix.
norm_graph_laplacian <- function(sim_mat){
  G <-  colSums(sim_mat) 
  N <- ncol(sim_mat)
  D_half = diag(1 / sqrt(G))
  return(diag(N) - D_half %*% sim_mat %*% D_half)
}

#### Spectral clustering function
### Inputs:
## dat: raw data matrix with features as rows and observations as columns.
## eigs: number of eigenvalues.
## clusts: number of cluster centers
spec_clust <- function(dat, eigs, clusts){
  
  dat <- row_normalize(dat)
  sim_mat <- similarity_mat(t(dat))
  
  graph_lap <- norm_graph_laplacian(sim_mat) 
  eigd <- eigen(graph_lap, symmetric = TRUE)
  rel_cols <- eigd$vectors[, (ncol(graph_lap) - eigs):(ncol(graph_lap) - 1)]
  sc_kmeans <- kmeans(rel_cols, clusts)
  return(sc_kmeans)
}
