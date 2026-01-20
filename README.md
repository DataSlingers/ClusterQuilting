# ClusterQuilting

This repository contains scripts for creating the results from "Cluster Quilting: Spectral Clustering for Patchwork Learning" by Zheng, Chang, and Allen (2025). 

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
Code for making figures that appear in Zheng, Chang, and Allen (2025). 
### PredictionValidation.R
Code for prediction validation evaluation of hyperparameters choices.
Method for hyperparameter selection as discussed in Zheng, Chang, and Allen (2025).

## Subfolders

### MICRONS
Contains instructions for downloading functional database from MICrONs repository.
Contains scripts for extracting fluoresence traces and metadata information from database.

### TCGA
Contains instructions for downloading individual files from TCGA website.
Contains scripts for building dataset from individual files.

### Workflows
Contains example workflows that can be used in analysis.
Examples for workflows used for the empirical studies in Zheng, Chang, and Allen (2025). 
Specific details for each file are discussed below.

## Workflow examples

This folder contains code pipelines and functions for producing results in Zheng, Chang, and Allen (2025). The following example workflows are available in the `workflows` subdirectory:

### Basic Simulation Workflow (BasicSimulationWorkflow.R)

1. `seq_create.R`/`mos_create.R`: Simulate GMM data with single set of simulation parameters.
2. `comprunner.m`: Get results from comparison methods.
3. `patch_ordering.R`: Get patch ordering for Cluster Quilting.
4. `ClusterQuilting.R`: Get estimates for Cluster Quilting, comparison methods.

### Full Simulation Study Workflow (SimulationStudyWorkflow.R)

1. `seq_create.R`/`mos_create.R`: Simulate GMM data across multiple sets of simulation parameters.
2. `comprunner.m`: Get results from comparison methods.
3. `compile_ari_oracle.R` / `compile_ari_dd.R`: Get estimates for Cluster Quilting; compile quantitative eval metrics for all methods.
4. `plot_func.R`: Create figures for showing performances between methods.

### TCGA Workflow (TCGAWorkflow.R)

1. `fuse_data.R`: Create census of overall downloaded files.
2. `TCGA_datamake.R`: Make a single data frame of all data out of individual files; create synthetic missingness patterns.
3. `comprunner.m`: Get results from comparison methods.
4. `compile_ari_oracle.R` / `compile_ari_dd.R`: Get estimates for Cluster Quilting; compile quantitative eval metrics for all methods.
5. `plot_func.R`: Create figures for showing performances between methods.
6. `TCGA_heatmap.R`: Create figures showing data imputation results.

### MICrONS Workflow (MICRONSWorkflow.R)

1. `trace_extract.py`: Get functional activity data from downloaded data.
2. `make_full.R`: Create single dataframe from multiple patches; create objects needed for estimation downstream.
3. `comprunner.m`: Get results from comparison methods.
4. `CompileARIMicrons.R`: Get estimates from Cluster Quilting.
5. `MICRONSfigures.R`: Create figures for showing cluster estimates from different methods. 