# ClusterQuilting

This repository contains scripts for creating the results from "Cluster Quilting: Spectral Clustering for Patchwork Learning" by Zheng, Chang, and Allen, 2025. 

## Individual Scripts

### HomogeneousSNRPatchOrdering.R
Finds optimal homogeneous SNR patch ordering based on Algorithm 4.1
### HeterogeneousSNRPatchOrdering.R
Finds optimal heterogeneous SNR patch ordering based on Algorithm S4.2
### ClusterQuilting.R
Cluster quilting estimation algorithm.
### MosaicPatchSim.R
Create simulated data with mosaic patch observation setting.
### SeqPatchSim.R
Create simulated data with sequential patch observation setting.
### comprunner.m
Contains scripts for running other comparison incomplete spectral clustering methods. 
Requires Matlab; code confirmed to work in version 2025b, probably back-compatible.
Adapted from code provided by authors of those methods. 
### SpectralClustering.R
Code for spectral clustering after imputation with NN imputation.
### CompileARI.R
Wrapper code for Cluster Quilting estimation, along with compiling ARI evaluation across all methods.
### MakeResultsPlots.R
Code for making figures that appear in Zheng, Chang, and Allen 2025. 
### PredictionValidation.R
Code for prediction validation evaluation of hyperparameters choices.
Method for hyperparameter selection as discussed in Zheng, Chang, and Allen 2025

## Subfolders

### MICRONS
Contains instructions for downloading functional database from MICrONs repository.
Contains scripts for extracting fluoresence traces and metadata information from database.

### TCGA
Contains instructions for downloading individual files from TCGA website.
Contains scripts for building dataset from individual files.

### Workflows
Contains example workflows that can be used in analysis.
Examples for workflows used for the empirical studies in Zheng, Chang, and Allen 2025. 
Specific details for each file are discussed below.

## Workflow examples

