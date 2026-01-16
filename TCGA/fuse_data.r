# Raw TCGA files contain data for each individual patient;
# This script fuses each data modality to a single data frame
## (Contains more modalities than what is used in the paper.)
library(tidyverse)

#####################################
#####################################
#####################################

data_dir <- "./" #Change to directory where data is saved locally
cancer_folders <- c("Kidney", "Lung", "Breast") # Assume that each project is saved in separate subfolder
save_dir <- "./output_dir/" #Change to directory where new data csv should be saved

# Create dictionary reference of all files downloaded.
## Used for matching rows in full data matrix further down the pipeline.
sample_list <- list()
sample_unique <- c()
data_types <- c()
cancer_types <- c()
for(d1 in 1:length(cancer_folds)){
  sample_list[[cancer_folds[d1]]] <- list()
  for(d2 in 1:length(data_view_folds)){
    old_len <- length(sample_unique)
    full_dir <- paste(data_dir, cancer_folds[d1], data_view_folds[d2], sep = "/")
    samples <- list.files(full_dir, pattern = "sample_sheet", full.names = TRUE)
    r1 <- read.delim(samples[1])
    cases_id <- unique(unlist(strsplit(r1$Case.ID, ", ")))
    sample_unique <- unique(c(sample_unique, cases_id))
    cancer_types <- c(cancer_types, rep(cancer_folds[d1], length(sample_unique) - old_len))
    sample_list[[cancer_folds[d1]]][[data_view_folds[d2]]] <- r1
    r1$data_type <- str_replace_all((paste(r1$Data.Category, r1$Data.Type, sep = "_")), " ", "")
    r1$cancer <- cancer_folds[d1]
    r1$dvfold <- data_view_folds[d2]
    if(!exists("r11")){
      r11 <- r1
    } else {
      r11 <- rbind(r11, r1)
    }
    data_types <- unique(c(data_types, str_replace_all(unique(paste(r1$Data.Category, r1$Data.Type, sep = "_")), " ", "")))
  }
}

view_df <- data.frame(samples = sample_unique, clust = cancer_types)
for(d3 in 1:length(data_types)){
  view_df[[data_types[d3]]] <- 0
}

#####################################

# Which modalities are available for which observations?
## Used for creating the full data matrix further down the pipeline.
for(d1 in 1:length(cancer_folds)){
  for(d2 in 1:length(data_view_folds)){
    full_dir <- paste(data_dir, cancer_folds[d1], data_view_folds[d2], sep = "/")
    samples <- list.files(full_dir, pattern = "sample_sheet", full.names = TRUE)
    r1 <- read.delim(samples[1])
    for(ind1 in 1:nrow(r1)){
      case_id <- unlist(strsplit(r1$Case.ID[ind1], ", "))[1]
      case_match <- which(view_df$samples == case_id)
      data_col <- str_replace_all(unique(paste(r1$Data.Category[ind1], r1$Data.Type[ind1], sep = "_")), " ", "")
      view_df[[data_col]][case_match] <- 1
    }
  }
}

write.csv(r11, paste0(data_dir, "/full_sample_sheet.csv")) 
write.csv(view_df, paste0(data_dir, "/view_miss.csv"))

#####################################
#####################################
#####################################

# Following loops create list of common features across all files of a modality.

#####################################

file_types <- list.files(data_dir, pattern = "level3betas.txt", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo], header=FALSE)
  if(length(all_betas) == 0){
    all_betas <- c(oop$V1)
  } else {
    all_betas <- intersect(all_betas, oop$V1)
  }
}

write.csv(all_betas, paste0(save_dir, "/all_dna_betas.csv"))


#####################################

file_types <- list.files(data_dir, pattern = "mirnas.quantification.txt", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  if(length(all_betas) == 0){
    all_betas <- c(oop$miRNA_ID)
  } else {
    all_betas <- intersect(all_betas, oop$miRNA_ID)
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_trans_mirna.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "isoforms.quantification.txt", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  if(length(all_betas) == 0){
    all_betas <- paste(oop$miRNA_ID, oop$isoform_coords, sep = "_")
  } else {
    all_betas <- intersect(all_betas, paste(oop$miRNA_ID, oop$isoform_coords, sep = "_"))
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_trans_isoforms.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "augmented_star_gene_counts.tsv", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo], comment.char="#")[-c(1:4), ]
  if(length(all_betas) == 0){
    all_betas <- c(oop$gene_id)
  } else {
    all_betas <- intersect(all_betas, oop$gene_id)
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_trans_genex.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "nocnv_grch38.seg.v2", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  if(length(all_betas) == 0){
    all_betas <- paste(oop$Chromosome, oop$Start, sep = "_")
  } else {
    all_betas <- intersect(all_betas, paste(oop$Chromosome, oop$Start, sep = "_"))
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_cnv_masked.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "\\.grch38.seg.v2", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  if(length(all_betas) == 0){
    all_betas <- paste(oop$Chromosome, oop$Start, sep = "_")
  } else {
    all_betas <- intersect(all_betas, paste(oop$Chromosome, oop$Start, sep = "_"))
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_cnv_seg.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "allelic_specific.seg", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  print(oo)
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo], comment.char="#")
  if(length(all_betas) == 0){
    all_betas <- paste(oop$Chromosome, oop$Start, sep = "_")
  } else {
    all_betas <- intersect(all_betas, paste(oop$Chromosome, oop$Start, sep = "_"))
    print(length(all_betas))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_cnv_allelic.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "gene_level_copy_number", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  rd <- which(is.na(oop$copy_number))
  single_file <- oop[-rd, ]
  if(length(all_betas) == 0){
    all_betas <- c(paste(oop$gene_id, oop$gene_name, sep = "_"))
  } else {
    all_betas <- base::intersect(all_betas, c(paste(oop$gene_id, oop$gene_name, sep = "_")))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_cnv_genelevel.csv"))

#####################################

file_types <- list.files(data_dir, pattern = "RPPA", recursive = TRUE, full.names = TRUE)

all_betas <- c()
for(oo in 1:length(file_types)){
  if(grepl("parcel", file_types[oo])){next}
  single_file <-  read.delim(file_types[oo])
  rd <- which(is.na(oop$copy_number))
  single_file <- oop[-rd, ]
  if(length(all_betas) == 0){
    all_betas <- c(paste(oop$gene_id, oop$gene_name, sep = "_"))
  } else {
    all_betas <- base::intersect(all_betas, c(paste(oop$gene_id, oop$gene_name, sep = "_")))
  }
}

write.csv(all_betas, paste0(save_dir, "/all_RPPA.csv"))
