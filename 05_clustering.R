---
title: "Clustering"
author: "Sonali Arora"
date: "Feb 20, 2025"
output:
  html_document:
    toc: true
    theme: united
---

In this vigentte, we show how we subclustered the ependymoma, shh medulloblastoma and group3/group4 medulloblastoma samples 

```{r setup}
library(ConsensusClusterPlus)
library(ggpubr)
library(mclust)
set.seed(123)  # For reproducibility

maindir = "C:/Users/sarora.FHCRC/Desktop/data_for_paper"
resdir = file.path(maindir, "figures_9_19_2024" )
setwd(maindir)
umap_data2 = read.delim("umap_data_vst_tpm_umap_9_16_2024.txt", 
                        header=T, stringsAsFactors = FALSE)
gene_exp_file = readRDS( file.path(maindir,"recalculated_log2tpm_for_med_epn.rds"))
# sanity check - ensure order is same for all matrices. 
gene_exp_file= gene_exp_file[ , match(umap_data2$sampleName, colnames(gene_exp_file))]
```

# consensus clustering 

```{r}
keep = which(umap_data2$disease  =="EPN")
chosen = umap_data2[keep, ]
your_data  = gene_exp_file[, keep]
expr_matrix <- as.matrix(your_data)

mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

# Set parameters
maxK <-15  # Maximum number of clusters to test
reps <- 100  # Number of resampling iterations
pItem <- 0.8  # Proportion of items to sample in each iteration
pFeature <- 1  # Proportion of features to sample (use 1 for all genes)
clusterAlg <- "hc"  # Hierarchical clustering (or use "km" for k-means)
distance <- "euclidean"  # Distance metric

# Run consensus clustering
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
    title="epn_consensus_cluster_5k_st_pf",
)
# k =2 chosen after looking at consensus matrix and delta area plot

# Gaussian mixture models
gmm_result <- Mclust(epn_umap_data2[, 1:2], G = 2)
finaldf_g <- cbind(epn_umap_data2, cluster = gmm_result$classification)
em1 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("Gaussian Mixture Models (GMM)") +
    scale_color_manual(values = c("red", "darkgreen", "grey90"))

# k-means 
kmeans_result <- kmeans(epn_umap_data2[, 1:2], centers = 2)
finaldf_g <- cbind(epn_umap_data2, cluster = kmeans_result$cluster)
em2 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("K-Means") +
    scale_color_manual(values = c("red", "darkgreen", "grey90"))

# hierarchial clustering
dist_matrix <- dist(epn_umap_data2[, 1:2])
hclust_result <- hclust(dist_matrix, method = "ward.D2")
finaldf_g <- cbind(epn_umap_data2, cluster = cutree(hclust_result, k = 2))
em3 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() +ggtitle("Hierarchial Clustering") +
    scale_color_manual(values = c("red", "darkgreen", "grey90"))

ggarrange(em1, em2, em3, ncol = 1, nrow =3)

````

## shh medulloblastoma

```{r}

keep = which(umap_data2$med_subgroup =="SHH")
chosen = umap_data2[keep, ]
your_data  = gene_exp_file[, keep]
expr_matrix <- as.matrix(your_data)  
mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

# Set parameters
maxK <- 20 # Maximum number of clusters to test
reps <- 100  # Number of resampling iterations
pItem <- 0.8  # Proportion of items to sample in each iteration
pFeature <- 1  # Proportion of features to sample (use 1 for all genes)
clusterAlg <- "hc"  # Hierarchical clustering (or use "km" for k-means)
distance <- "euclidean"  # Distance metric

# Run consensus clustering
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
    title="shh_consensus_cluster_5k_outlier",
)
# k =3 chosen after looking at consensus matrix and delta area plot

# Gaussian mixture models
gmm_result <- Mclust(chosen[, 1:2], G = 3)
finaldf_g <- cbind(try, cluster = gmm_result$classification)
t1 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("Gaussian Mixture Models (GMM)") +
    scale_color_manual(values = col_vector[1:length(unique(gmm_result$classification))])

# k-means 
kmeans_result <- kmeans(chosen[, 1:2], centers = 3)
finaldf_g <- cbind(try, cluster = kmeans_result$cluster)
t2 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("K-Means") +
    scale_color_manual(values =col_vector[1:length(unique(finaldf_g$cluster))])

# hierarchial clustering
dist_matrix <- dist(chosen[, 1:2])
hclust_result <- hclust(dist_matrix, method = "ward.D2")
finaldf_g <- cbind(try, cluster = cutree(hclust_result, k = 3))
t3 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() +ggtitle("Hierarchial Clustering") +
    scale_color_manual(values = col_vector[1:length(unique(finaldf_g$cluster))])

ggarrange(t1, t2, t3,  ncol = 1, nrow =3)

```


## Group3 / Group4 medulloblastoma

```{r}
keep = which(umap_data2$med_subgroup  %in% c("Group3", "Group4" ))
chosen = umap_data2[keep, ]
your_data  = gene_exp_file[, keep]
expr_matrix <- as.matrix(your_data)  # Ensure it's a matrix

mads=apply(expr_matrix, 1, mad)
expr_matrix=expr_matrix[rev(order(mads))[1:5000],]

# Set parameters
maxK <- 10  # Maximum number of clusters to test
reps <- 100  # Number of resampling iterations
pItem <- 0.8  # Proportion of items to sample in each iteration
pFeature <- 1  # Proportion of features to sample (use 1 for all genes)
clusterAlg <- "km"  # Hierarchical clustering (or use "km" for k-means)
distance <- "euclidean"  # Distance metric

# Run consensus clustering
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
    title="med_g3_g4_consensus_cluster_5k",
)

# k =6 chosen after looking at consensus matrix and delta area plot

# Gaussian mixture models
gmm_result <- Mclust(chosen[, 1:2], G = 6)
finaldf_g <- cbind(try, cluster = gmm_result$classification)
g1 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("Gaussian Mixture Models (GMM)") +
    scale_color_manual(values = col_vector[1:length(unique(gmm_result$classification))])

# k-means 
kmeans_result <- kmeans(chosen[, 1:2], centers = 6)
finaldf_g <- cbind(try, cluster = kmeans_result$cluster)
g2 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() + ggtitle("K-Means") +
    scale_color_manual(values =col_vector[1:length(unique(finaldf_g$cluster))])

# hierarchial clustering
dist_matrix <- dist(chosen[, 1:2])
hclust_result <- hclust(dist_matrix, method = "ward.D2")
finaldf_g <- cbind(try, cluster = cutree(hclust_result, k = 6))
g3 <- ggplot(finaldf_g, aes(UMAP1_2d, UMAP2_2d, color=factor(cluster))) +
    geom_point(size=2) + theme_bw() +ggtitle("Hierarchial Clustering") +
    scale_color_manual(values = col_vector[1:length(unique(finaldf_g$cluster))])


ggarrange(g1, g2, g3, ncol = 1, nrow =3)

```





