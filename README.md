# Introduction 

This github repository contains code to reproduce the analysis in our upcoming paper "Transcriptomic landscape analysis identifies regional biology and tumor heterogeneity in medulloblastoma and uncovers two major unrecognized ependymoma subtypes" from the Holland Lab at Fred Hutch Cancer Center.

## System Requirements 

All analaysis was done in R , using version 4.4.0

## Instructions 

The following scripts reproduce the analysis done in the paper. 

1. This [01_ProcessingBulkRNASeqData.R](01_ProcessingBulkRNASeqData.R) reads in data , batch corrects it and normalizes it using 2 approaches and creates 2d and 3d umaps
2. Once the umap is created, this [02_PlotLandscape.R](02_PlotLandscape.R) visualizes avaiable metadata and recreates the fig1 from the paper
3. In this file[03_Genewise_analysis_plots.R](03_Genewise_analysis_plots.R), we recreate gene wise plots for a given gene, we plot  Gene exp, Copy no, Gene fusion for a given gene in different figures for the paper.
4. GSVA analysis was done using various datasets such as Reactome, Biocarta and KEGG , this script (04_Pathway.Rmd)[04_Pathway.Rmd] shows how to visualize the scores over the umap 
5. This [05_clustering.R](05_clustering.R) shows how we sub-clustered the different subtypes for ependymoma and medulloblastoma
6. This script [06_FusionDotplot.R](06_FusionDotplot.R) looks at gene fusions across medulloblastoma and ependymoma and recreastes fig5 in the paper.
7. This script [07_DEGanalysis.R](07_DEGanalysis.R) contains code on how to do DEG analysis for the ependymoma datasets
8. This script [08_LandPatient.R](08_LandPatient.R) shows how you can use our estabilished landscape to land new patients onto the landscape , to better understand their biology.
9. Finally, this script [09_scRNASeq_medullo.Rmd](09_scRNASeq_medullo.Rmd) takes in publicly avaiable scRNASeq data to validate the novel pathways we find in medulloblastoma, using bulk RNASeq data.

## Data availability

The harmonized, gene expression data for the paper is available under [releases](https://github.com/sonali-bioc/MedulloEPNLandscape/releases)

## Oncoscape, an interactive website to visualize the dataset
The landscape is accessible via Oncoscape, an open-source platform. The figures to the paper can be visualized on Oncoscape in 3D here:

[Oncoscape, for Figures 1 through 5.](https://oncoscape.sttrcancer.org:&project=medulloblastomaumap33)

After the data finishes loading, you can click the blue "Figures from the Paper" dropdown menu to choose a sub-figure.

![Screenshot of figures menu](https://i.imgur.com/ZnVBJQd.png)

Figure 6 is an interactive volcano plot where you can look at individual genes of interest and where they land on the volcano plot
[Figure 6 (E1 vs. E2 Volcano Plot)](https://oncoscape-apps.vercel.app/presentations/medullo-volcanos)
