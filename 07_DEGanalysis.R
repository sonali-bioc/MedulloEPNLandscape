---
title: "Reproducing Differential expression analysis and volcano plots  from the paper"
author: "Sonali Arora"
date: "Sept 23, 2024"
output:
  html_document:
    toc: true
    theme: united
---
In this file, we first conduct a Differential expression analysis ( DEG) using E1 and E2 clusters using R package DESeq2. 
Next, we also make volcano plots which are included in the paper

# DEG analysis 

```{r}
count_matrix_raw = readRDS("raw_batch_corrected_with_only_fetal.Rds")
umap_data2 = read.delim("umap_data_vst_tpm_umap_9_16_2024.txt", 
                        header=T, stringsAsFactors = FALSE)

# sanity check to ensure samples are in correct order
identical(umap_data2$sampleName, colnames(count_matrix_raw))


my_deseq2_fun = function( rawdata, group1, group2, group1_idx, group2_idx, resdir){    
    coldata = data.frame(sampleGroup = c(rep("control", length(group1_idx)), 
                                         rep("test", length(group2_idx))), 
                         sampleName = colnames(rawdata) )  
    rownames(coldata) = colnames(rawdata)
    dds <- DESeqDataSetFromMatrix(countData =round (rawdata), 
                                  colData = coldata, design = ~ sampleGroup)    
    keep <- rowSums(counts(dds)) >= 4
    dds = dds[keep, ]    
    dds = DESeq(dds)    
    lfc = log2(1.25)
    res1 <- results(dds, alpha = 0.05,  contrast = c("sampleGroup", "test", "control"))    
    resdf = as.data.frame(res1)    
    up_genes = resdf[ which(resdf$log2FoldChange > lfc &  resdf$padj < 0.05), ]
    down_genes = resdf[which(resdf$log2FoldChange < -lfc &  resdf$padj < 0.05), ]      
    up_genes = up_genes[order(up_genes$log2FoldChange, decreasing=T), ]
    down_genes = down_genes[order(down_genes$log2FoldChange), ]  
  
    # write deseq2 results    
    lst = list(all_genes = resdf, up_genes = up_genes, down_genes = down_genes)
    fname =paste0(c("DESeq2_analysis_", group2, "_vs_", group1), collapse ="")
    write_xlsx(lst, file.path(resdir, paste0(fname  , ".xlsx")))
    message(fname, ":", nrow(up_genes), ":", nrow(down_genes))       
}

my_umap_data2 = umap_data2
adjusted = count_matrix_raw

group1 = "E2"
group2 = "E1"
group1_idx = which(my_umap_data2$chosen==group1)
group2_idx = which(my_umap_data2$chosen==group2)
rawdata = cbind(adjusted[,group1_idx], adjusted[, group2_idx])
folder_name = paste0(c("DESeq2_", group2, "_vs_", group1), collapse="")
resdir = file.path(bigresdir, folder_name)
e1 = my_deseq2_fun(rawdata, group1, group2, group1_idx, group2_idx, resdir )


group1 = "E2"
group2 = "E1"
group2a = which(my_umap_data2$chosen==group1 & umap_data2$fusion =="ZFTA_RELA")
group2b = which(my_umap_data2$chosen==group2 & umap_data2$fusion =="ZFTA_RELA")
rawdata = cbind(adjusted[,group2a], adjusted[, group2b])
folder_name = "DESeq2_E1_zftarela_vs_E2_zftarela"
resdir = file.path(bigresdir, folder_name)
e2 = my_deseq2_fun(rawdata, group1, group2, group2a, group2b, resdir )


group1 = "E1_ZFTARELA"
group2 = "E1_nonzftarela"
group3a = which(my_umap_data2$chosen=="E1" & umap_data2$fusion !="ZFTA_RELA")
group3b = which(my_umap_data2$chosen=="E1" & umap_data2$fusion =="ZFTA_RELA")
rawdata = cbind(adjusted[,group3a], adjusted[, group3b])
folder_name = "DESeq2_within_E1_zftarela_vs_non_zftarela"
resdir = file.path(bigresdir, folder_name)
e3 = my_deseq2_fun(rawdata, group1, group2, group3a, group3b, resdir )

```

# Volcano plots

```{r}

kinases = read_xlsx("Kinases_info.xlsx", sheet = 1)
kinases = as.data.frame(kinases)

lfc_cutoff = log2(1.25)
all_genes = read_xlsx("DESeq2_E1_vs_E2/DESeq2_analysis_E1_vs_E2.xlsx", sheet = 1)
up_genes = read_xlsx("DESeq2_E1_vs_E2/DESeq2_analysis_E1_vs_E2.xlsx", sheet = 2)
down_genes = read_xlsx("DESeq2_E1_vs_E2/DESeq2_analysis_E1_vs_E2.xlsx", sheet = 3)
up_genes = as.data.frame(up_genes)
down_genes = as.data.frame(down_genes)
up_kinases = up_genes[match( intersect(kinases[,1], up_genes[,1]), up_genes[,1]), ]
up_kinases = up_kinases[order(up_kinases$log2FoldChange, decreasing=T), ]
down_kinases = down_genes[match( intersect(kinases[,1], down_genes[,1]), down_genes[,1]), ]
down_kinases = down_kinases[order(down_kinases$log2FoldChange, decreasing=F), ]
res = as.data.frame(all_genes)
res$log10padj = -log10(res$padj)
res$delabel = NA_character_
res$delabel[ match(c(up_kinases[,1], down_kinases[,1]), res[,1])] = "KINASES"

chosen_genes = c("NTRK2", "NTRK3", "MERTK", "EPHB4",  "FGFR1", "FGFR3", "MET")

pdf("Volcano_plot_e1_e2.pdf", width = 5, height =5 )
with(res, plot(log2FoldChange, log10padj, pch=20, col ="grey80",main=paste0("Volcano plot: E1_vs_E2") ))
with(subset(res, padj<.05 & log2FoldChange>lfc_cutoff), points(log2FoldChange, log10padj, pch=20, col="pink"))
with(subset(res, padj<.05 & log2FoldChange < -lfc_cutoff), points(log2FoldChange, log10padj, pch=20, col="lightgreen"))
with(subset(res[which(!is.na(res$delabel)), ]), points(log2FoldChange, log10padj, pch=20, col="steelblue3"))
abline(v = lfc_cutoff , col = "black", lty = 2)
abline(v = -lfc_cutoff, col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)
temp = res[which(!is.na(res$delabel)), ]
temp$status =  "NOCHANGE"
temp$status[which(temp$padj < 0.05 & abs(temp$log2FoldChange) > lfc_cutoff )] ="sig"
with(temp[which(temp[,1] %in% chosen_genes), ], points(log2FoldChange, log10padj, pch=20,  col="black"))
with(temp[which(temp[,1] %in% chosen_genes), ], textxy(log2FoldChange,log10padj, labs=gene, cex=1) ) # label the genes
dev.off()
```
