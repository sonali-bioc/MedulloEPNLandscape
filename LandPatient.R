---
title: "Reproducing Figure 7  from the paper"
author: "Sonali Arora"
date: "Sept 23, 2024"
output:
  html_document:
    toc: true
    theme: united
---
In this file, we reproduce Fig7  from the paper. 

First, we calculated vst counts for the 12 unclassified MB tumors. 
Next, we ran the nearest neighbor pipeline from Thirimane et al. 
The code for which can be obatined from  https://github.com/FredHutch/MeningiomaLandscape-HollandLab
Finally, the python code produced a file called "medullo_centroid.csv" which we overlayed and plotted with our reference landscape. 

```{r}
umap = read.delim("umap_data_vst_tpm_umap_6_14_2024.txt", header=T, stringsAsFactors = FALSE)
centroid = read.csv("medullo_centroid.csv")
rownames(centroid) = centroid[,1]
colnames(centroid) = c("sampleName", "UMAP1_2d", "UMAP2_2d")

centroid$new_name = toupper(letters[1:12])
umap$new_name = ""

umap$color = "grey70"
centroid$color = "black"
centroid$med_subgroup = "black"
centroid$Taylor_subgroups  = "black"
newdf = rbind(centroid, umap[, c("sampleName", "UMAP1_2d", "UMAP2_2d", "new_name", "color", "med_subgroup", "Taylor_subgroups")])
newdf = newdf[order(newdf$color, decreasing=T), ]

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

newdf$med_subgroup = factor(newdf$med_subgroup, 
                            levels = c("Group3", "Group4", "SHH","WNT", "black" , 
                                       "EPN",  "fetal"))
mydf = newdf
mydf$new_name[which(mydf$new_name =="A")] = "A,G"
mydf$new_name[which(mydf$new_name =="G")] = ""

plot2 = ggplot(mydf,aes(UMAP1_2d, UMAP2_2d, color=med_subgroup ,label = new_name, alpha = color, size  = color )) +
    geom_point() +
    scale_color_manual(values = c(col_vector[c(2,3, 5,6)] ,"black","lightblue" , "forestgreen" ) ) +
    scale_size_manual(values = c( 3, 2)) +
    scale_alpha_manual(values = c( 1, 0.5)) +
    theme_bw() +  
    ggtitle(paste0("Medullo subgroups"))  +
    guides(color=guide_legend(ncol=1)) + 
    geom_text(hjust=0, vjust=-0.5, size = 6)

pdf(file.path(resdir,"main_fig7_WINNER_UMAP.pdf"), width =9, height = 7)
print(plot2)
dev.off()


```
In the supplemental figure 7, we show 1 umap for each of the unclassified sample - showing the ground truth as well as the predicted coordinate.

```{r}
my_samples = rownames(centroid)
y1 = lapply(my_samples, function(x){
    ggdf = newdf
    ggdf$mycolor = "none"
    ggdf$mycolor[which( ggdf$sampleName==x & ggdf$new_name =="")] = "original"
    ggdf$mycolor[which( ggdf$sampleName==x & ggdf$new_name !="")] = "predicted"
    ggdf$mycolor = factor(ggdf$mycolor, levels = c("original", "predicted", "none"))
    
    title_name = paste0(centroid[which(rownames(centroid)==x), ]$new_name, ":", x)
    
    p1 = ggplot(ggdf %>% arrange(desc(is.na(mycolor))), aes(UMAP1_2d, UMAP2_2d, color=mycolor , alpha = mycolor, size  = mycolor )) +
        geom_point() +
        scale_color_manual(values = c(  "red", "black", "grey80")) +
        scale_size_manual(values = c( 4, 4,2)) +
        scale_alpha_manual(values = c( 1, 1, 0.5)) +
        theme_bw() +  
        ggtitle(title_name)  + xlab("UMAP1") + ylab("UMAP2") +
        guides(color=guide_legend(ncol=1)) +  theme_clean
    p1
})

pdf(file.path(resdir, "Supp_Fig7_prediction_9_25_2024.pdf"), width = 18, height = 14)
ggarrange(plotlist = y1 ,
          nrow = 3, ncol = 3 , 
          common.legend = TRUE, legend="right")
dev.off()

```
Finally , we show the survival analysis for MYC and EIF4EBP1

```{r}
library(survminer)
library(survival)
log2_tpm =gene_exp_file

idx = which(umap_data2$med_subgroup=="Group3")
temp_gene = log2_tpm[, idx]
ridx = grep("^MYC$", rownames(log2_tpm)) 
higher = which(temp_gene[ridx,] > 7)
lower = which(temp_gene[ridx,] < 7)
higher_samples = colnames(temp_gene)[higher]
lower_samples = colnames(temp_gene)[lower]

group1_idx = match(higher_samples, colnames(log2_tpm))
group2_idx = match(lower_samples, colnames(log2_tpm))

my_chosen_tpm_data = cbind(log2_tpm[,group1_idx], log2_tpm[, group2_idx])
chosen_umap_data = rbind(umap_data2[group1_idx, ], umap_data2[group2_idx, ])
chosen_umap_data$myc_status = "Myc_low_group3"
chosen_umap_data$myc_status [match(higher_samples, chosen_umap_data$sampleName)] = "Myc_amp_group3"

mysurvdata = chosen_umap_data
mysurvdata$os.months = as.numeric( mysurvdata$os.months)
mysurvdata$status =mysurvdata$os_event
mysurvdata$status = as.numeric(mysurvdata$status)
fit1 <- surv_fit(Surv(os.months, status) ~ myc_status, data = mysurvdata)

p1 = ggsurvplot(fit1, size = 2,         
                palette = c("red", "orange", "grey70"),                
                risk.table.height = 0.3, 
                ylab = "", xlab = "Time (months)",
                ggtheme = theme_bw()      
) + ggtitle(paste0("Myc amp vs Myc low: group3"))

```
Similar analysis for EIF4EB1
```{r}
idx = which(umap_data2$med_subgroup=="Group3")
temp_gene = log2_tpm[, idx]
ridx = grep("^EIF4EBP1$", rownames(log2_tpm)) 
higher = which(temp_gene[ridx,] > 7)
lower = which(temp_gene[ridx,] < 7)
higher_samples = colnames(temp_gene)[higher]
lower_samples = colnames(temp_gene)[lower]

group1_idx = match(higher_samples, colnames(log2_tpm))
group2_idx = match(lower_samples, colnames(log2_tpm))

my_chosen_tpm_data = cbind(log2_tpm[,group1_idx], log2_tpm[, group2_idx])
chosen_umap_data = rbind(umap_data2[group1_idx, ], umap_data2[group2_idx, ])
chosen_umap_data$myc_status = "EIF4EBP1_low_group3"
chosen_umap_data$myc_status [match(higher_samples, chosen_umap_data$sampleName)] = "EIF4EBP1_amp_group3"

mysurvdata = chosen_umap_data
mysurvdata$os.months = as.numeric( mysurvdata$os.months)
mysurvdata$status =mysurvdata$os_event
mysurvdata$status = as.numeric(mysurvdata$status)
fit1 <- surv_fit(Surv(os.months, status) ~ myc_status, data = mysurvdata)

p2 = ggsurvplot(fit1, size = 2,         
                palette = c("red", "orange", "grey70"),                
                risk.table.height = 0.3, 
                ylab = "", xlab = "Time (months)",
                ggtheme = theme_bw()      
) + ggtitle(paste0("EIF4EBP1 amp vs EIF4EBP1 low: group3"))

pdf(file.path(resdir, "Fig7_survival_analysis.pdf"), width = 5 , height =4)
print(p1)
print(e1)
dev.off()

```
    

