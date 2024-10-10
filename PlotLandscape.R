---
title: "Reproducing Fig1 from Paper"
author: "Sonali Arora"
date: "Oct 9, 2024"
output:
  html_document:
    toc: true
    theme: united
---
In this file, we reproduce Fig1 from the paper

```{r setup}

rm(list=ls())

library(edgeR)
library(readxl)
library(umap)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(data.table)
library(magrittr)
library(ggpubr)

1 = 2# size for points in PCA plot
legend_pt_size =4
plot_title_size = 25
axis_text_size = 25
axis_title_size=25
legend_text_size=15
plot_text = 10
spacing=1
chosen_margin = c(0.5,1,0.5,1)# margins:top,right,bottom,left
theme_clean <- theme_bw() +
    theme(
        plot.title = element_text(hjust=0, vjust=0, 
                                  lineheight=.8, face="bold", size=plot_title_size ),
        plot.margin=unit(chosen_margin,"cm"), 
        legend.text=element_text(size=legend_text_size),
        legend.key.height = unit(spacing, "cm"),
        legend.position = "inside", legend.position.inside = c(0.8, 0.2) ,
        legend.title=element_blank() , 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

maindir = "C:/Users/sarora.FHCRC/Desktop/data_for_paper"
resdir = file.path(maindir, "figures_9_19_2024" )
setwd(resdir)

umap_data2 = read.delim(file.path(maindir, "umap_data_vst_tpm_umap_9_16_2024.txt"), 
                        header=T, stringsAsFactors = FALSE)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


umap_data2$dataset[which(umap_data2$dataset=="EPN-Heidelberg(Mack)")] = "EPN-Heidelberg"
umap_data2$dataset[which(umap_data2$dataset=="EPN-StJude(Mack)")] = "EPN-StJude"
umap_data2$dataset[which(umap_data2$dataset=="EPN-Toronto(Mack)")] = "EPN-Toronto"
dataset_plot= ggplot(umap_data2 %>% arrange(desc(dataset)),  aes(UMAP1_2d, UMAP2_2d, color=dataset )) +
    geom_point(size=s1) + 
    scale_color_manual(values = c(col_vector[1:7], "lightpink")) +
    theme_clean + 
    guides(color=guide_legend(ncol=1)) 

pub_plot= ggplot(umap_data2 %>% arrange(desc(Source)),  
                 aes(UMAP1_2d, UMAP2_2d, color=Source )) +
    geom_point(size=s1) + 
    scale_color_manual(values = c(col_vector[1:7], "lightpink")) +
    theme_clean +  
    guides(color=guide_legend(ncol=1)) 


umap_data2$Characteristics.organism = factor(umap_data2$Characteristics.organism, 
                                             levels = c("forebrain", "hindbrain", "EPN", "Med"))
fetal_type_plot= ggplot(umap_data2,aes(UMAP1_2d, UMAP2_2d, color=Characteristics.organism )) +
    geom_point(size=s1) + 
    scale_color_manual(values = c("purple", "yellowgreen", "lightblue", "lightsalmon")) +
    theme_clean +  
    guides(color=guide_legend(ncol=1))


umap_data2$diagnosis = umap_data2$med_subgroup
umap_data2$diagnosis = factor(umap_data2$diagnosis, 
                               levels = c("fetal", "Group3", "Group4", "SHH","WNT", "To be classified" , 
                                          "EPN"))
med_subtype_plot_v2= ggplot(umap_data2, aes(UMAP1_2d, UMAP2_2d, color=diagnosis)) +
    geom_point(size=2) +
    scale_color_manual(values = c("forestgreen", "#BEAED4", "#FDC086",  "plum2", "#F0027F" , "grey90",
                                  "#E6AB02" ,"#386CB0",  "coral",
                                  "black",  "red" ,"purple"  )) + theme_clean 
 
umap_data2$EPN_PF_subgroup[which(umap_data2$EPN_PF_subgroup=="Not available")] =NA_character_
umap_data2$EPN_PF_subgroup[which(umap_data2$EPN_PF_subgroup=="fetal")] = NA_character_
umap_data2$EPN_PF_subgroup = factor(umap_data2$EPN_PF_subgroup, levels = c( "PF-A", "PF-B", "PF", "PF_SE",
                                                                            "Medulloblastoma"))
pf_plot= ggplot(umap_data2  %>% arrange(desc(is.na(EPN_PF_subgroup))),
                aes(UMAP1_2d, UMAP2_2d, color=EPN_PF_subgroup )) +
    geom_point(size=s1) + 
    scale_colour_manual(values = c("red", "blue", "black", "magenta" ,
                                   "lightsalmon"), na.value = "grey90") + 
    theme_clean + 
    guides(color=guide_legend(ncol=1))

umap_data2$ST_EPN_subtype = as.character(umap_data2$EPN_subtype)
umap_data2$ST_EPN_subtype[which(umap_data2$ST_EPN_subtype!="Supratentorial EPN")] = "Others" 
umap_data2$ST_EPN_subtype[which(umap_data2$ST_EPN_subtype=="Supratentorial EPN")] = "ST-EPN" 
only_st_epn_plot = ggplot(umap_data2  %>% arrange(ST_EPN_subtype),
                          aes(UMAP1_2d, UMAP2_2d, color=ST_EPN_subtype )) +
    geom_point(size=s1) + scale_color_manual(values = c("grey90", col_vector[14]), na.value = "grey90") +
    theme_clean  #ggtitle(paste0("Supratentorial EPN"))  +

umap_data2$gender[which(umap_data2$gender=="Not available")]=NA_character_
umap_data2$gender = factor(umap_data2$gender, levels = c("F", "M", "fetal" ,NA_character_))
gender_plot= ggplot(umap_data2  %>% arrange(desc(is.na(gender))), 
                    aes(UMAP1_2d, UMAP2_2d, color=gender )) +
    geom_point(size=s1) + scale_color_manual(values =c("magenta","skyblue1", "forestgreen", "grey70"), na.value ="grey90" ) +
    theme_clean 

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
sc2 <- scale_colour_gradientn(colours = myPalette(12), limits=c(0, 40), na.value = "grey90")
umap_data2 = umap_data2[ order(umap_data2$age_in_years), ]
tumor_age_plot= ggplot(umap_data2,aes(UMAP1_2d, UMAP2_2d, color=age_in_years )) +
    geom_point(size=s1) + 
    geom_point(data =umap_data2[which(!is.na(umap_data2$age_in_years)), ], aes(UMAP1_2d, UMAP2_2d ))+
    sc2 + 
    theme_clean

fetal_myPalette <- colorRampPalette(rev(brewer.pal(9, "Purples")))
fetal_sc <- scale_colour_gradientn(colours = fetal_myPalette(25), limits=c(4, 25), na.value = "lightgoldenrod1")
ggdf2 = umap_data2
ggdf2$alpha_order = "red"
ggdf2$alpha_order [which(is.na(umap_data2$fetal_age_forebrain))] = "black"

ggdf2 = ggdf2[order(ggdf2$fetal_age_forebrain, decreasing=T), ]
ggdf2$alpha_order  = factor(ggdf2$alpha_order , levels = c("red", "black"))
fetal_age_plot1= ggplot(ggdf2 %>% arrange(desc(is.na(fetal_age_forebrain))), 
                        aes(UMAP1_2d, UMAP2_2d, color=fetal_age_forebrain)) +
    geom_point(size=s1) + 
    fetal_sc + theme_clean + guides(alpha = "none") 

fetal_myPalette <- colorRampPalette(rev(brewer.pal(9, "Greens")))
fetal_sc <- scale_colour_gradientn(colours = fetal_myPalette(25), limits=c(4, 25), na.value = "lightgoldenrod1")
fetal_age_plot2= ggplot(ggdf2 %>% arrange(desc(is.na(fetal_age_hindbrain))), 
                        aes(UMAP1_2d, UMAP2_2d, color=fetal_age_hindbrain ))+ 
    geom_point(size=s1) + 
    fetal_sc + 
    theme_clean

pdf(file.path(resdir, "Multiplot_fig1_new_umap_colors_9_16_2024_9panel.pdf"),width =22, height = 16)
ggarrange(pub_plot , fetal_type_plot ,  med_subtype_plot, 
          pf_plot , only_st_epn_plot , gender_plot ,
          tumor_age_plot , fetal_age_plot1 ,fetal_age_plot2,
          ncol = 3 ,  nrow =3)
dev.off()
```
