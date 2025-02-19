---
title: "Reproducing Gene exp, Copy no, Gene fusion plots from the paper"
author: "Sonali Arora"
date: "Oct 9, 2024"
output:
  html_document:
    toc: true
    theme: united
---
In this file, we recreate gene wise plots for a given gene, we plot  Gene exp, Copy no, Gene fusion for a given gene in different figures for the paper.

```{r setup}

rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

s1 = 2# size for points in PCA plot
legend_pt_size =4
plot_title_size = 25
axis_text_size = 25
axis_title_size=25
legend_text_size=25
spacing=1
chosen_margin = c(0.5,1,0.5,1)# margins:top,right,bottom,left
theme_clean <- theme_bw() +
    theme(
        plot.title = element_text(hjust=0, vjust=0, 
                                  lineheight=.8, face="bold", size=plot_title_size ),
        plot.margin=unit(chosen_margin,"cm"), 
        legend.text=element_text(size=legend_text_size),
        legend.key.height = unit(spacing, "cm"),
        legend.position = "right",
        legend.justification = 'left',
        legend.spacing.y = unit(0.5, 'cm'),
        legend.title=element_blank() , 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

vals = c("chocolate",  "grey90")  # for gene fusions
names(vals) = c(1, 0)

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu"))) # for gene expression

vals2 = c("darkred",  "grey90",  "darkblue") # for copy no plots
names(vals2) = c(1, 0, -1)

maindir = "C:/Users/sarora.FHCRC/Desktop/data_for_paper"
resdir = file.path(maindir, "figures_9_19_2024" )
setwd(maindir)
finalChrMat = readRDS(file.path(maindir, "updated_finalChrMat_3_5_2024.Rds")) # copy no 
gene_fusions = readRDS(file.path(maindir, "gene_fusions_9_18_2024.rds"))
score = readRDS(file.path(maindir, "GSVA_score_v2_only_fetal.rds"))
umap_data2 = read.delim("umap_data_vst_tpm_umap_9_16_2024.txt", 
                        header=T, stringsAsFactors = FALSE)
gene_exp_file = readRDS( file.path(maindir,"recalculated_log2tpm_for_med_epn.rds"))

# sanity check - ensure order is same for all matrices. 
gene_exp_file= gene_exp_file[ , match(umap_data2$sampleName, colnames(gene_exp_file))]
score = score[, match(umap_data2$sampleName, colnames(score))]
finalChrMat = finalChrMat[match(umap_data2$sampleName, rownames(finalChrMat)), ]

identical(umap_data2$sampleName, colnames(gene_exp_file))
identical(rownames(finalChrMat) , umap_data2$sampleName)
identical(colnames(score), umap_data2$sampleName)

my_cn_plot  = function(arm_list){
    lapply(arm_list, function(arm){
        x  = match(arm, colnames(finalChrMat))
        cn = finalChrMat[, x]
        cn[is.na(cn)] =0
        ggdf = umap_data2
        ggdf$cn = 0
        ggdf$cn = cn        
        ggdf = ggdf[ order(abs(ggdf$cn), decreasing=FALSE), ]
        ggdf$cn = factor(ggdf$cn, levels = names(vals2))        
        p1= ggplot(ggdf, aes(UMAP1_2d, UMAP2_2d, color=cn )) +
            geom_point(size=3) +
            scale_color_manual(values = vals2, labels = names(vals2), na.value = "grey90")+
            ggtitle(paste0(colnames(finalChrMat)[x])) + 
            theme_clean + theme(legend.position ="none")
        p1
    })
}

my_gene_fun = function(goi){
    goi = intersect(goi, rownames(gene_exp_file))    
    l1 = lapply(goi, function(my_goi){
        message(my_goi)
        mat= umap_data2
        mat$gene_exp = unlist(gene_exp_file[ match(my_goi, rownames(gene_exp_file)), ])              
        mat = mat[ order(mat$gene_exp, decreasing=FALSE), ]
        ylim = max(mat$gene_exp)
        ylow = min(mat$gene_exp )
        sc <- scale_colour_gradientn(colours = myPalette(ylim), limits=c(ylow, ylim))        
        fig1 =  ggplot(mat, aes(UMAP1_2d, UMAP2_2d, color=gene_exp )) +
            geom_point(size=3) + ggtitle(paste0( my_goi)) + 
            sc + theme_clean +  
            theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2) )
        
        fig1
    })
}

my_fusion_plot  = function(chosen_partners){
    lapply(chosen_partners, function(x){
        gene1  = unlist(strsplit(x, "_"))[1]
        gene2  = unlist(strsplit(x, "_"))[2]
        rev_x = paste0(gene2, gene1, collapse="_")        
        gidx = c( which(gene_fusions$gene_partner==x), which(gene_fusions$gene_partner==rev_x) )
        my_fusions = gene_fusions[gidx, ]
        my_samples = unique(unlist(strsplit(my_fusions[, "sampleName"], ", ")) )        
        ggdf = umap_data2
        ggdf$fusion = 0
        midx = match(my_samples, umap_data2$sampleName)
        ggdf$fusion[midx] = 1        
        my_title = paste0(gene1, ":", gene2)
        ggdf = ggdf[ order(abs(ggdf$fusion), decreasing=FALSE), ]
        ggdf$fusion = factor(ggdf$fusion, levels = names(vals))
        p1= ggplot(ggdf, aes(UMAP1_2d, UMAP2_2d, color=fusion)) +
            geom_point(size=3) + xlab("") + ylab("") +
            scale_color_manual(values = vals, labels = names(vals), na.value = "grey90")+
            ggtitle(my_title) + 
            theme_clean + theme(legend.position ="none") 
        p1
    })
}    

# main figure2 

a = my_cn_plot("9q")
b = my_cn_plot("17p")
c = my_cn_plot("17q")
d = my_gene_fun("ATOH1")
e = my_gene_fun("MYC")
plot_f =  my_gene_fun("KCNA1")
g  = my_fusion_plot("ZFTA_RELA")
g = lapply(g, remove_titles)
h = my_fusion_plot("YAP1_MAMLD1")
h = lapply(h, remove_titles)
i =  my_cn_plot("9p")
j = my_cn_plot("1q")
k = my_cn_plot("6q")
l = my_cn_plot("22q")
m = my_gene_fun("WNT5A")
n = my_gene_fun("TGFB1")
o = my_gene_fun("IGF2")
pdf(file.path(resdir, "Multiplot_fig2_validate_9_16_2024.pdf"), width = 35, height = 21)
ggarrange( a[[1]], b[[1]], c[[1]], d[[1]], e[[1]], plot_f[[1]], 
           g[[1]], h[[1]], i[[1]], j[[1]],k[[1]], l[[1]], m[[1]],
           n[[1]],o[[1]], nrow =3, ncol =5)
dev.off()

# main fig 6 
a1 = my_gene_fun(c("NTRK2", "NTRK3", "MERTK", "EPHB4"))
pdf(file.path(resdir, "Multiplot_fig6_9_16_2024_recalc_tpm.pdf"), width = 7, height  =21)
ggarrange( plotlist = a1, nrow =4, ncol =1)
dev.off()

# main fig 7 
o = my_gene_fun("EIF4EBP1")[[1]]
pdf(file.path(resdir,"main_fig7_9_16_2024.pdf"), width = 7, height = 5)
print(o)
dev.off()

## supp fig2 : 
row1 = my_gene_fun( c("SFRP1", "HHIP", "GABRA5","IMPG2", "EOMES" ))
row2 = my_gene_fun(  c( "RBM24", "HOXB2", "L1CAM","CCND1","CREBBP"))
row3 =my_fusion_plot(c("CCDC196_LINC02290", "GJE1_VTA1" ,"PVT1_PCAT1", "TUBB2B_LMAN2L", "ELP4_IMMP1L"))
row3 = lapply(row3, remove_titles)
pdf(file.path(resdir,"Multiplot_Supp_fig2_validate_9_16_2024.pdf"), width = 35, height = 21)
ggarrange( plotlist = c(row1, row2, row3), nrow =3, ncol =5)
dev.off()

## Supp Fig 5 
supp_fig5 = c( "SALL2_METTL3", "SERPINI2_WDR49" , "SLC1A1_SPATA6L", toupper("frmpd2_ptpn20cp"),  
               "TIAM2_SCAF8", "NEDD1_CFAP54",
               "EP400_ZNF471", "CRACR2B_GTF2H1", "DCP2_COMT", "EPB41L4A_NREP", #"BARHL2_SNORDG2", 
               "LINC01419_SCAPER", "DDX31_GFI1B", "EOMES_ADGRV1", #"FSCN1P1_DOLDA8A", 
               "PRDM6_C11orf71", "BTG4_POU2AF1", "EPHA8_ABCF1")
b1 =  my_fusion_plot(supp_fig5)
pdf(file.path(resdir,"Multiplot_Supp_fig5_other_fusions.pdf"), width = 28, height = 21)
ggarrange( plotlist = b1, nrow =4, ncol =4)
dev.off()

## Supp Fig 6 - synaptic gene expression plots
syn_genes  = c("GRIN1" , "CHRNB1", "CHRNA4", "GRM2", "P2RX5", 
               "GABRA5", "DRD1", "GRIN2C", "P2RX7")
syn_plot = my_gene_fun( syn_genes)
pdf(file.path(maindir,"Multiplot_Supp_fig6_syn_genes_9_16_2024.pdf"), width = 21, height = 18)
ggarrange( plotlist = syn_plot, nrow =3, ncol =3)
dev.off()
```
