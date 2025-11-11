# A Pipeline for Analysis of NanoString GeoMx Spatial Transcriptomics:

Steps Included in this analysis: 

1. Quality Control using NanoString GeoMx pipeline
2. Normalization using quantile normalization, TMM, or log cpm normalization provided by voom
3. Global visualization using UMAP, tSNE, and heatmap, with confidence interval ellipses between treatment groups
4. Differential Expression Analysis using TMM and voom normalization with generalized estimation equations (GEE)
5. Weighted Gene Correlation Network Analysis (WGCNA) using TMM and voom normalization with generalized estimation equation (GEE)
6. Visualizations for Differential Expression results and WGCNA, including gene-module membership

## To use and edit the code in this repository: 

```bash

git clone https://github.com/estagaman/NanoStringGeoMx.git

cd NanostringGeoMx

```

## What you will need: 

datadir: The path to a directory containing:

  - .dcc files from Nanostring GeoMx sequencer
  - .pkc file, also provided by NanoString. This can be downloaded from their website directly
  - .xlsx file of the sample metadata, containing column "Sample_ID" to match to DCC file identifiers

out_dir: for each script, you will specify the output directory you would like to save any files to

## Required R Packages: 

NanoStringNCTools
GeomxTools
GeoMxWorkflows
limma 
edgeR
preprocessCore
tidyverse
dplyr
glmtoolbox
geesmv
WGCNA
umap
rtsne
pheatmap
readxl
ggplotify
grid
data.table       
