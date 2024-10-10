---
title: "Reproducing Gene Fusion Dotplots from the paper"
author: "Sonali Arora"
date: "Oct 9, 2024"
output:
  html_document:
    toc: true
    theme: united
---
In this file, we reproduce the dotplots made using gene fusions results from the paper.

```{r}
rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggpubr)

maindir = "C:/Users/sarora.FHCRC/Desktop/data_for_paper"
resdir = file.path(maindir, "figures_9_19_2024" )
gene_fusions = readRDS(file.path(maindir, "gene_fusions_9_18_2024.Rds"))
umap_data2 = read.delim("umap_data_vst_tpm_umap_9_16_2024.txt", 
                        header=T, stringsAsFactors = FALSE)

ggdf = umap_data2
tada1 = sapply(ggdf$sampleName, function(x){
    idx = grep(x, gene_fusions$sampleName)
    length(idx)
})
cut_want = cut(as.numeric(as.character(tada1)), breaks = c(0, 5,15, 25 , 50,  150, 410))
cut_want[which(is.na(cut_want)) ] = "(0,5]"
ggdf$gene_fusions_freq2 =cut_want
ggdf= ggdf[order(ggdf$gene_fusions_freq2, decreasing=FALSE), ]
gf_freq =  ggplot(ggdf, aes(UMAP1_2d, UMAP2_2d, color=gene_fusions_freq2 )) +
    geom_point(size=s1) + 
    scale_colour_manual(values = c("grey90", "green", "forestgreen", 
                                   "pink","mediumvioletred",  "black" ), na.value ='grey90') +
    theme_bw() + theme(legend.position = "inside", legend.position.inside = c(0.7, 0.2) )


pdf("Fig3_gene_fusions_freq_6_14_2024.pdf", width = 7, height = 5)
print(gf_freq)
dev.off()


chosen_gene_fusions = read.delim("chosen_gene_fusions.txt", header=T, stringsAsFactors=F)
l3 =lapply(chosen_gene_fusions[-1,1], function(x){
    gene1  = unlist(strsplit(x, "_"))[1]
    gene2  = unlist(strsplit(x, "_"))[2]
    rev_x = paste0(gene2, gene1, collapse="_")
    
    gidx = c( which(gene_fusions$gene_partner==x), which(gene_fusions$gene_partner==rev_x) )
    my_fusions = gene_fusions[gidx, ]
    my_samples = unique(unlist(strsplit(my_fusions[, "sampleName"], ", ")) )
        
    ## all st-epn in e1
    want = ggdf[which(ggdf$EPN_subtype=="Supratentorial EPN" & ggdf$chosen =="E1"), ]$sampleName
    a1 = length(intersect(want, my_samples))
    
    # all st-epn in e2
    want2 = ggdf[which(ggdf$EPN_subtype=="Supratentorial EPN" & ggdf$chosen =="E2"), ]$sampleName
    a2 = length(intersect(want2, my_samples))
    
    # st-eon in e2 with rela-zfta fusion
    want3 = ggdf[e2_zfta_Rela_st_epn, ]$sampleName
    a3 = length(intersect(want3, my_samples))
    
    ## All posterior fossa FUSIONS?
    want4 = ggdf[which(ggdf$EPN_subtype=="Posterior Fossa EPN" & ggdf$chosen =="E2"), ]$sampleName
    a4 = length(intersect(want4, my_samples) )
    
    ## All posterior fossa FUSIONS?
    want5 = ggdf[which(ggdf$EPN_subtype=="Posterior Fossa EPN" & ggdf$chosen =="E1"), ]$sampleName
    a5 = length(intersect(want5, my_samples))
    
    ## pf-b FUSIONS?
    want6= ggdf[which(ggdf$EPN_PF_subgroup=="PF-B"), ]$sampleName
    a6 = length(intersect(want6, my_samples))
    
    ## pf-a FUSIONS?
    want7 = ggdf[which(ggdf$EPN_PF_subgroup=="PF-A"), ]$sampleName
    a7 = length(intersect(want7, my_samples))
    
    ## anaplastic 
    want8 = ggdf[which(ggdf$EPN_subtype=="Myxopapillary EPN" ), ]$sampleName
    a8 = length(intersect(want8, my_samples))
    
    want9 = ggdf[which(ggdf$EPN_subtype=="Spinal EPN"), ]$sampleName
    a9 = length(intersect(want9, my_samples))
    
    want10 = ggdf[which(ggdf$EPN_subtype=="Anaplastic EPN" & ggdf$chosen =="E2"), ]$sampleName
    a10 =length(intersect(want10, my_samples))
    
    want11 = ggdf[which(ggdf$EPN_subtype=="Anaplastic EPN" & ggdf$chosen =="E1"), ]$sampleName
    a11 = length(intersect(want11, my_samples))
    
    want12 = ggdf[which(ggdf$med_subgroup=="Group3"), ]$sampleName
    a12 = length(intersect(want12, my_samples))
    
    want13 = ggdf[which(ggdf$med_subgroup=="Group4"), ]$sampleName
    a13 = length(intersect(want13, my_samples))
    
    want14 = ggdf[which(ggdf$med_subgroup=="SHH"), ]$sampleName
    a14 = length(intersect(want14, my_samples))
    
    want15 = ggdf[which(ggdf$dataset=="fetal"), ]$sampleName
    a15 = length(intersect(want15, my_samples))
    
    c(a1/length(want)*100, a2/length(want2)*100, a3/length(want3)*100, a4/length(want4)*100, a5/length(want5)*100, 
      a6/length(want6)*100, a7/length(want7)*100, a8/length(want8)*100, a9/length(want9)*100, 
      a10/length(want10)*100, a11/length(want11)*100, a12/length(want12)*100, 
      a13/length(want13)*100, a14/length(want14)*100, a15/length(want15)*100)    
})

l3 = do.call(rbind, l3)
l3[is.na(l3)] =0
colnames(l3) = c( "ST.EPN_E1"  ,    "ST.EPN_E2"  ,  "ST.EPN.ZFTA.RELA.E2", "PF.EPN_E2",   "PF.EPN_E1"   ,       
                  "PF.B"   ,"PF.A"  , "Myxopapillary_EPN", "Spinal_EPN" ,"Anaplastic_EPN_E2"  ,
                  "Anaplastic_EPN_E1"  ,  "Group3" ,"group4"  ,   "shh" , "fetal")
rownames(l3) = chosen_gene_fusions[-1,1]
adf = melt(l3)
adf$group = "EPN E1"
colnames(adf)[1]="fusion"
colnames(adf)[2]="variable"

adf$group[grep("E2", adf$variable)] ="EPN E2"
adf$group[grep("PF.A|PF.B", adf$variable)] ="EPN E2"
adf$group[grep("shh|group4|Group3", adf$variable)] ="Medullo"
adf$group[grep("Myxopapillary_EPN", adf$variable)] ="EPN E2"
adf$group[grep("fetal", adf$variable)] ="Fetal"

adf$disease = "Medullo"
adf$disease[grep("shh", adf$variable)] = "SHH"
adf$disease[grep("group4", adf$variable)] = "Group4"
adf$disease[grep("Group3", adf$variable)] = "Group3"
adf$disease[grep("ST.EPN", adf$variable)] = "ST-EPN"
adf$disease[grep("PF.EPN", adf$variable)] = "PF-EPN"
adf$disease[grep("PF.A", adf$variable)] = "PF-A"
adf$disease[grep("PF.B", adf$variable)] = "PF-B"
adf$disease[grep("Spinal", adf$variable)] = "Spinal-EPN"
adf$disease[grep("Anaplastic", adf$variable)] = "Anaplastic-EPN"
adf$disease[grep("Myxopapillary", adf$variable)] = "Myxopapillary-EPN"
adf$disease[grep("fetal", adf$variable)] = "Fetal"

adf$disease = factor(adf$disease, levels = c("ST-EPN", "PF-EPN", "PF-A", "PF-B", "Anaplastic-EPN", "Myxopapillary-EPN", "Spinal-EPN", 
                                             "SHH", "Group3", "Group4", "Fetal"))

adf$group = factor(adf$group, levels = c("EPN E1", "EPN E2", "Medullo", "Fetal"))

f0 = ggplot(adf , aes(x = disease, y = fusion, size = value)) + geom_point( color = "red") + 
    xlab("") + ylab("") +
    facet_grid(~group, scales = "free", space="free")+ 
    guides(size=guide_legend(title="% of patients \nwith fusion")) +
    theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1, size = 16),
          axis.text.y= element_text(size = 16),
          legend.title = element_text(size = 16),
          strip.text = element_text(size = 14,  face="bold"))+ scale_size(range = c(0,8))

pdf(file.path(resdir, "Fig5_dotplot_gene_fusion_9_16_2024.pdf"), width = 10, height = 12)
print(f0)
dev.off()

```
