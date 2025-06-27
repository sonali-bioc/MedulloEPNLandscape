---
title: "Reproducing Pathway Figures from the paper"
author: "Sonali Arora"
date: "June 26, 2025"
output:
  html_document:
    toc: true
    theme: united
---

In this file, we reproduce the consensus clustering shown in the paper

# Set up
```{r}
rm(list=ls())
df = readRDS("raw_batch_corrected_with_only_fetal.Rds")
umap_data = read.delim("umap_data_vst_tpm_umap_9_16_2024.txt", header=T, stringsAsFactors=F)
df = df[ , match(umap_data$sampleName, colnames(df))]
identical( colnames(df), umap_data$sampleName)

# Load required libraries
library(DESeq2)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(ggplot2)
library(circlize )

# Set parameters
maxK <- 10  # Maximum number of clusters to test
reps <- 100  # Number of resampling iterations
pItem <- 0.8  # Proportion of items to sample in each iteration
pFeature <- 1  # Proportion of features to sample (use 1 for all genes)
clusterAlg <- "km"  # Hierarchical clustering (or use "km" for k-means)
distance <- "euclidean"  # Distance metric

#function for making Consensus Matrix Heatmap
plotConsensusMatrix <- function(results, k) {
  consensus_matrix <- results[[k]]$consensusMatrix
  Heatmap(consensus_matrix, 
          name = "Consensus Index", 
          col = colorRampPalette(c("white", "blue"))(50),
          show_row_names = FALSE, 
          show_column_names = FALSE,
          cluster_rows = TRUE, 
          cluster_columns = TRUE)
}
# Custom heights for annotated heatmap
annotation_height_custom <- unit(c(8, 8, 8), "mm")  # Make annotation bars taller
heatmap_height_custom <- unit(50, "mm")             # Shrink heatmap a bit
```

# Ependymoma
```{r}
epn_Idx = which(umap_data$disease=="EPN")
epn_raw = df[, epn_Idx]
epn_raw  = round(epn_raw)
epn_vst = DESeq2::vst(epn_raw)
chosen = umap_data[epn_Idx, ]
expr_matrix <- as.matrix(epn_vst)  # Ensure it's a matrix
mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

results3 <- ConsensusClusterPlus(
  expr_matrix,
  maxK = maxK,
  reps = reps,
  pItem = pItem,
  pFeature = pFeature,
  clusterAlg = clusterAlg,
  distance = distance,
  seed = 1234,  # Set seed for reproducibility
  plot = "png",
  title="epn_consensus_cluster_5kgenes",
)

clusters_result3 <- results3[[2]]$consensusClass
table(clusters_result3)

# Plot the consensus matrix for a chosen k (e.g., k=3)
plotConsensusMatrix(results3, k=2)

# Annotated figure for ependymoma 
k <- 2
consensus_matrix <- results3[[k]]$consensusMatrix

# Define color palettes
diagnosis_colors <- c(
  "Supratentorial EPN"  = "#E6AB02", 
  "Posterior Fossa EPN" = "#386CB0",
  "EPN Tumor"           = "coral",
  "Myxopapillary EPN"   = "black",
  "Anaplastic EPN"      = "red",
  "Spinal EPN"          = "purple"  
)
fusion_colors <- c(
  "ZFTA_RELA" = "#CC00FF", 
  "YAP1_MAMLD1" = "black", 
  "other ZFTA" = "chocolate", 
  "other YAP1" = "blue", 
  "Not Available" = "#D3D3D3"
)
chosen$fusion_clean <- ifelse(is.na(chosen$fusion), "Unknown", chosen$fusion)
cluster_colors <- c("E1" = "darkred", "E2" = "yellowgreen")
cluster_labels <- factor(results3[[2]]$consensusClass, levels = c(1,2), labels = c("E1", "E2"))

# Annotation block
ha <- HeatmapAnnotation(
  Diagnosis = chosen$diagnosis,
  Fusion = chosen$fusion_clean,
  Cluster = cluster_labels,
  col = list(
    Diagnosis = diagnosis_colors,
    Fusion = fusion_colors,
    Cluster = cluster_colors
  ),
  which = "column",
  annotation_height = annotation_height_custom
)
# Consensus heatmap
ht = Heatmap(
  consensus_matrix,
  name = "Consensus Index",
  col = colorRamp2(c(0, 1), c("white", "blue")),
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  height = heatmap_height_custom,
  heatmap_legend_param = list(
    title_position = "topcenter",  # Move title to the left
    legend_direction = "horizontal",
    legend_height = unit(40, "mm")
  )
)
pdf("EPN_consensus_heatmap_final.pdf", width =8, height = 6)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

```

# Medulloblastoma - shh
```{r}
shh_Idx = which(umap_data$med_subgroup=="SHH")
shh_raw = df[, shh_Idx]
shh_raw  = round(shh_raw)
shh_vst = DESeq2::vst(shh_raw)
chosen = umap_data[shh_Idx, ]

expr_matrix <- as.matrix(shh_vst)  # Ensure it's a matrix
mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

results2 <- ConsensusClusterPlus(
   expr_matrix,
   maxK = maxK,
   reps = reps,
   pItem = pItem,
   pFeature = pFeature,
   clusterAlg = clusterAlg,
   distance = distance,
   seed = 1234,  # Set seed for reproducibility
   plot = "png", 
   title="shh_consensus_cluster_5kgenes",
)
clusters_result2 <- results2[[3]]$consensusClass
table(clusters_result2)

k <- 3
consensus_matrix <- results2[[k]]$consensusMatrix

cluster_colors <- c("S1" = "chocolate", "S2" = "grey20",  "S3" = "forestgreen")

# Annotation block
ha <- HeatmapAnnotation(
  Diagnosis = chosen$diagnosis,
 # Taylor = chosen$Taylor_subgroups,
  Cluster = cluster_labels,
  col = list(
    Diagnosis = c("SHH" = "plum2"),
     Cluster = cluster_colors
  ),
  which = "column",
  annotation_height = annotation_height_custom
)
# Consensus heatmap
ht = Heatmap(
  consensus_matrix,
  name = "Consensus Index",
  col = colorRamp2(c(0, 1), c("white", "blue")),
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  height = heatmap_height_custom,
  heatmap_legend_param = list(
    title_position = "topcenter",  # Move title to the left
    legend_direction = "horizontal",
    legend_height = unit(40, "mm")
  )
)

pdf("shh_consensus_heatmap_final.pdf", width =8, height = 6)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
```

# Medulloblastoma - Group3 and Group4
```{r}
g3g4_Idx = which(umap_data$med_subgroup %in% c( "Group3", "Group4"))
g3g4_raw = df[, g3g4_Idx]
g3g4_raw  = round(g3g4_raw)
g3g4_vst = DESeq2::vst(g3g4_raw)
chosen = umap_data[g3g4_Idx, ]

expr_matrix <- as.matrix(g3g4_vst)  # Ensure it's a matrix
mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

results1 <- ConsensusClusterPlus(
   expr_matrix,
   maxK = maxK,
   reps = reps,
   pItem = pItem,
   pFeature = pFeature,
   clusterAlg = clusterAlg,
   distance = distance,
   seed = 1234,  # Set seed for reproducibility
   plot = "png", 
   title="g3g4_consensus_cluster_5kgenes",
)

clusters_result1 <- results1[[6]]$consensusClass
table(clusters_result1)


k <- 6
consensus_matrix <- results1[[k]]$consensusMatrix

cluster_colors = c( "#7FC97F", "#BEAED4" ,"#FDC086", "#FFFF99" ,"#386CB0" ,"#F0027F" )
names(cluster_colors) = c("C1", "C2", "C3", "C4", "C5", "C6")
cluster_labels <- factor(clusters_g3g4, levels = c("C1", "C2", "C3", "C4", "C5", "C6"))

# Annotation block
ha <- HeatmapAnnotation(
  Diagnosis = chosen$diagnosis,
  #Fusion = chosen$fusion_clean,
  Cluster = cluster_labels,
  col = list(
    Diagnosis = c("Group3" = "red", "Group4"= "black"),
    Cluster = cluster_colors
  ),
  which = "column",
  annotation_height = annotation_height_custom
)

# Consensus heatmap
ht = Heatmap(
  consensus_matrix,
  name = "Consensus Index",
  col = colorRamp2(c(0, 1), c("white", "blue")),
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  height = heatmap_height_custom,
  heatmap_legend_param = list(
    title_position = "topcenter",  # Move title to the left
    legend_direction = "horizontal",
    legend_height = unit(40, "mm")
  )
)

pdf("g3_g4_consensus_heatmap_final.pdf", width =8, height = 6)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
```
