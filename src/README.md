# Systematic comparison reveals that transcriptional differences between in vivo and ex vivo hematopoietic model systems modulate the outcomes of genetic perturbations.
This folder contains the R code to analyze Perturb-seq CRISPR screens. R scripts in this folder are organized by the type of analysis.
The scripts from steps SCRNA_01 to SCRNA_07 are modified from https://github.com/csbg/tfcf/tree/main/src
To run this code you need the Singularity container,  and then use renv to install the packages in the file renv.lock.

## Contents
- Scripts to obtain external data, start with EXT
- Scripts to analyze Perturb-seq data start with SCRNA
    - Basic analysis (QC, integration, duplet detection, cell type identification, analysis of marker genes) in files SCRNA_01 to SCRNA_07
    - SCRNA_08 to SCRNA_12 is used for differential expression analysis using limma and enrichment analysis
- SCRNA_13 to SCRNA_14 for comparison of gene expression within and between datasets
-  SCRNA_15 to SCRNA_16 (Differential expression analyses and correlation of actual(Perturb-seq dataset) versus biolord-predicted KO effects)
-  SCRNA_17 ( Differential expression analysis and enrichment analysis for splenic macrophages and T-cells)
-  Figures- numbered according to figures in the manuscript
- Additional functions are defined in 
    - FUNC_ProjecTILs_PLUS.R
    - Ag_top_filtered_genes.R
    - Ag_enrichR_mouse_genes.R
    - Ag_Optimized_theme_fig.R
