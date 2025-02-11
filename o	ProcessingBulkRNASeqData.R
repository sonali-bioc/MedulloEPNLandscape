---
title: "Processing bulk RNASeq  data"
author: "Sonali Arora"
date: "Feb 10, 2025"
output:
  html_document:
    toc: true
    theme: united
---

In this vigentte, we show how we combine the datasets, batch correct , normalize and do dimension reduction to build the UMAP

```{r}
# read in the CBTN data
big_mat = readRDS("~/datasets/CBTN_medullo_raw_gene_Exp_93samples.rds")
sampleinfo_cbtn = read.delim("~/HollandLabShared/Sonali/Taylor_group/datasets/medullo_subtype_brain_umap_cbtn.txt", 
                             header=T, stringsAsFactors = FALSE)
sampleinfo_cbtn = sampleinfo_cbtn[match( colnames(big_mat), sampleinfo_cbtn$sampleName), 
                                  c("sampleName", "molecular_subtype","pathology", "os_days", "os_event",  "gender", "age_in_years", "race")]

# read in taylor group data
taylor_raw =readRDS( "~/datasets/all_Samples_combined_8_14_2023.rds")
colnames(taylor_raw) = gsub(".gene_id.exon.no.htseq.txt", "", colnames(taylor_raw))
sampleinfo = read.delim("~/datasets/Hendrikse_et_al_libraries_with_Scratch_path_filesize.txt")
rownames(sampleinfo) = sampleinfo[,2]
sampleinfo = sampleinfo[colnames(taylor_raw), ]
sampleinfo = sampleinfo[, c("Hendrikse_subgroup", "Histology", "OS", "death", "sex", "age") ] # race is missing 

## fetal
fetal = readRDS("~/datasets/fetal_reverse_protein_coding.Rds")

## get epn data 
felix_epn = readRDS("~/datasets/Felix_raw_protein_coding_genes.rds")
toronto_epn = readRDS("~/datasets/Toronto_raw_protein_coding_genes.rds")
heidelberg_epn =readRDS("~/datasets/Heidelberg_raw_protein_coding_genes.rds")
mack_epn = readRDS("~/datasets/Mack_stjude_raw_protein_coding_genes.rds")
epn_counts = cbind(felix_epn, toronto_epn, heidelberg_epn, mack_epn)

epn_cbtn = readRDS("~/datasets/CBTN_EPN_raw_gene_Exp_77samples.rds")
sampleinfo_cbtn2 = readRDS("~/datasets/epn_info.Rds")
sampleinfo_cbtn2 = sampleinfo_cbtn2[match( colnames(epn_cbtn), sampleinfo_cbtn2$sampleName), 
             c("tumor_type","disease", "os_days", "os_event",  "gender", "age_in_years", "race", "sampleName")]

# annotations
library(rtracklayer)
gtf = import("~/GATK_hg38_reference/gencode.v39.primary_assembly.annotation.gtf")
gtf = gtf[which(gtf$type=="gene") , ]
gtf = gtf[which(gtf$gene_type=="protein_coding"), ]
uni_gene = unique(gtf$gene_name)
gtf = gtf[ match(uni_gene, gtf$gene_name), ]

taylor_raw = taylor_raw[gtf$gene_id, ]
rownames(taylor_raw) =gtf$gene_name

# use common genes from each dataset
common_genes = intersect( rownames(big_mat), rownames(taylor_raw))

taylor_raw = taylor_raw[common_genes, ]
big_mat = big_mat[common_genes, ]
epn_counts = epn_counts[common_genes, ]
epn_cbtn = epn_cbtn[common_genes, ]
fetal = fetal[common_genes, ]

# batch correction 
library(sva)
count_matrix <- cbind( taylor_raw, big_mat, epn_counts, epn_cbtn, fetal)
count_matrix = data.matrix(count_matrix)
batch = c(rep("med-batch1", ncol(taylor_raw)), rep("med-batch2", ncol(big_mat)),
          paste0("epn-",epn_coldata$batch), rep("epn-batch5", ncol(epn_cbtn)) , 
           rep("fetal", ncol(fetal)))

adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)

dataset = c(rep("Med-Taylor", ncol(taylor_raw)), rep("Med-CBTN", ncol(big_mat)),
            paste0("EPN-", epn_coldata$dataset ), rep("EPN-CBTN", ncol(epn_cbtn)), 
           rep("fetal", ncol(fetal))) 

disease =  c(rep("Med", ncol(taylor_raw)), rep("Med", ncol(big_mat)),  
             rep("EPN", ncol(epn_counts)), rep("EPN", ncol(epn_cbtn)), 
             rep("fetal", ncol(fetal))) 

# VST normalization
dataMtrx = DESeq2::vst(adjusted)

#2d umap ( based on VST) 
umap_out <- umap(t(dataMtrx), random_state = 123, min_dist = 0.5) # Run UMAP
umap_2d = umap_out$layout
colnames(umap_2d) = c("UMAP1_2d", "UMAP2_2d")

# 3d umap ( based on VST) 
custom.settings = umap.defaults
custom.settings$n_components = 3
umap_out = umap(t(dataMtrx), config = custom.settings, random_state = 123, min_dist = 0.5)
umap_3d = umap_out$layout
colnames(umap_3d) = c("UMAP1_3d", "UMAP2_3d", "UMAP3_3d")

# TPM normalization 
getTPM = function(exon_mat, gtf){    
    gw = width(gtf)
    y <- DGEList(exon_mat)
    y$genes <- data.frame(Length=gw)
    rpkm_vals = rpkm(y )
    log2_rpkm_vals = log2(rpkm_vals +1)
    tpm_vals = apply(rpkm_vals, 2, function(x){
        (x/sum(x))*10^6
    })
    log2_tpm_vals = log2(tpm_vals +1)
    log2_tpm_vals
}
combat_edata2 = getTPM(adjusted, gtf)

#2d umap ( based on log2 TPM ) 
umap_out <- umap(t(combat_edata2), random_state = 123, min_dist = 0.5) # Run UMAP
umap_2d_tpm = umap_out$layout
colnames(umap_2d_tpm) = c("UMAP1_2d_TPM", "UMAP2_2d_TPM")

# 3d umap ( based on log2 TPM) 
custom.settings = umap.defaults
custom.settings$n_components = 3
umap_out = umap(t(combat_edata2), config = custom.settings, random_state = 123, min_dist = 0.5)
umap_3d_tpm = umap_out$layout
colnames(umap_3d_tpm) = c("UMAP1_3d_TPM", "UMAP2_3d_TPM", "UMAP3_3d_TPM")

# save result files
ans_df =data.frame(  umap_2d, umap_3d, umap_2d_tpm, umap_3d_tpm, sampleName = colnames(dataMtrx))
write.table(ans_df, "umap_data_vst_tpm_umap_9_16_2024.txt", sep ="\t", quote=F, row.names=F, col.names=T)

```
