#REQUIRED PACKAGES 
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
devtools::install_github(repo = 'mojaveazure/seurat-disk', ref = 'develop')
remotes::install_github("FASTGenomics/fgread-r")
install.packages("Seurat")
devtools::install_github("immunogenomics/harmony") #if you have problem installing it, see: https://github.com/immunogenomics/symphony/issues/2
BiocManager::install("celldex")

#LIBRARIES
library(fgread)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(celldex)
library(dplyr)
library(scales)
source("/Users/otavio.cmarques/Documents/Otavio/Orientacoes /Lena/Covid Datasets/10x PBMC/Usage/Scripts/SEURAT VISUALISATION.R")

#Import Data Set
# load data from rds file
seurat_c1_full_FG <- readRDS("/Users/otavio.cmarques/Documents/Otavio/Orientacoes /Lena/Covid Datasets/10x PBMC/Data/seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds")

#Group vs Sex
ggplot(seurat_c1_full_FG@meta.data, aes(x = group_per_sample , fill= sex))+
  geom_bar()+
  geom_text(stat='count', aes(label=..count..), position = position_stack(vjust = 0.5))+
  theme_classic()+
  theme(
    panel.grid=element_blank(),
    legend.text=element_text(size=10),
    text = element_text(size=12),
    legend.title = element_blank(),
    axis.title.x = element_blank()
  )+  
  ylab("# of cells")+ 
  RotatedAxis()

#QC (FALTA INSERIR O PLOT PARA MOSTRAR A PORCENTAGEM MITOCONDRIAL)
#seurat_c1_full_FG[["percent.mt"]] <- PercentageFeatureSet(seurat_c1_full_FG, pattern = "^MT-")
#seurat_c1_full_FG <- subset(seurat_c1_full_FG, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## NORMALIZATION
#seurat_c1_full_FG <- NormalizeData(object = seurat_c1_full_FG, normalization.method = "LogNormalize", scale.factor = 1e4)
## DEFINE VARIABLE GENES
#seurat_c1_full_FG <- FindVariableFeatures(object = seurat_c1_full_FG, assay="RNA", selection.method = 'vst', nfeatures = 10000)
## SCALING
#seurat_c1_full_FG <- ScaleData(object = seurat_c1_full_FG, vars.to.regress = c("nCount_RNA"))
#seurat_c1_full_FG <- ScaleData(object = seurat_c1_full_FG, features = features)

# RUN SCTRANSFORM (INSTEAD OF NORMALIZATION AND SCALING)
#seurat_c1_full_FG <- SCTransform(seurat_c1_full_FG, vars.to.regress = "nCount_RNA", verbose = FALSE, return.only.var.genes = T, ncells = 500)

## PCA
#seurat_c1_full_FG <- RunPCA(object = seurat_c1_full_FG,  verbose = FALSE)

## UMAP
#seurat_c1_full_FG <- RunUMAP(seurat_c1_full_FG, reduction = "pca", dims = 1:20, seed.use = 42)

## HARMONY
#seurat_c1_full_FG <- RunHarmony(seurat_c1_full_FG, group.by.vars = "orig.ident", reduction = "pca", dims.use = 1:20)

## HARMONY UMAP
#seurat_c1_full_FG <- RunUMAP(seurat_c1_full_FG, reduction = "harmony", dims = 1:20, seed.use = 42)

## FIND NEIGHBORS
#seurat_c1_full_FG <- FindNeighbors(object = seurat_c1_full_FG, dims = 1:20, reduction="harmony", force.recalc = TRUE)

## CALCULATE CLUSTERS
#seurat_c1_full_FG <- FindClusters(object = seurat_c1_full_FG, resolution = 0.4, algorithm = 1)
 
 DimPlot(object = seurat_c1_full_FG, 
        reduction = 'umap',
        label = F,
        group.by = "group_per_sample",
        cols=c("#eaac43","#38c3e1","#1657e1"))+ggtitle("UMAP")+
   scale_color_manual("Group", values = c(mild = "#38c3e1", severe= "#1657e1",  alpha(c(control="#eaac43"),0.5)))+
   theme(axis.text.x = element_text(size = 17),
         axis.text.y = element_text(size = 17),  
         axis.title.x = element_text(size = 17),
         axis.title.y = element_text(size = 17),
         title =element_text(size=17, face='bold'))+NoLegend()

 DimPlot(object = seurat_c1_full_FG,
        reduction = 'umap',
        group.by="RNA_snn_res.0.4",
        label = TRUE, repel = TRUE) + scale_color_manual(values=color_clusters) + NoLegend()
   

#Figure 2A
 Idents(seurat_c1_full_FG)<-"RNA_snn_res.0.4"
 seurat_c1_full_FG <- RenameIdents(seurat_c1_full_FG,
                                   `0`= "0",
                                   `7`= "1",
                                   `12`="2",
                                   `6`= "3",
                                   `8`= "4",
                                   `11`="5",
                                   `13`="6",
                                   `17`="7",
                                   `18`="8" ,              
                                   `1`= "9",
                                   `3`= "10",
                                   `10`="11",
                                   `2`= "12",
                                   `15`="13",
                                   `9`= "14",
                                   `4`= "15" ,     
                                   `5`= "16",
                                   `19`="17",
                                   `20`="18",
                                   `16`="19" ,      
                                   `14`="20",
                                   `21`="21",
                                   `22`="22"
 )
 seurat_c1_full_FG$new.order <- Idents(seurat_c1_full_FG)
 
 seurat_c1_full_FG <- RenameIdents(seurat_c1_full_FG,
                                   
                                   `0`="Classical Monocytes",
                                   `7`="HLA-DR CD83 Monocytes",
                                   `12`="CD163 Monocytes",
                                   `6`="HLA-DR S100A Monocytes",
                                   `8`="Non-classical Monocytes",
                                   `11`="Neutrophils",  
                                   `13`="Immature Neutrophils",
                                   `17`="mDCs",
                                   `18`="pDCs",
                                   `1`="CD4 T",
                                   `3`="CD4 T",
                                   `10`="CD4 T",
                                   `2`="CD8 T",
                                   `15`="CD8 T",
                                   `9`="CD8 T",
                                   `4`="NK",
                                   `5`="B",
                                   `19`="B",
                                   `20`="B",
                                   `16`="Plasmablasts",
                                   `14`="Megakaryocyte",
                                   `21`="mixed",
                                   `22`="undefined"
 )
 seurat_c1_full_FG$cluster_labels_res.0.4 <- Idents(seurat_c1_full_FG)
 
 
 umap.colors<- c(
   "#7D6C86", #classic
   "#59BF30",  #mono zfp
   "#9D43BB", #IFI6+ mono
   "#5E2870", #s100+ mono
   "#DB8A0F",  #cd16+
   "#638B83", #LD1
   "#A6E9DB", #LD2
   "#FBD64A", #mDC
   "#B4DC49", #21 pDC
   "#DE342F", #cd4
   "#CF5046", #cd4
   "#FB6C46", #cd4
   "#0F95B9",#8 CD8 TnTcm
   "#0F95DA",#7 CD8
   "#6FD6E8", # prol.T
   "#1F405C", #NK
   "#A34F23", #B
   "#A35A33", #B
   "#A33A43", #B
   "#F2895E", #plas
   "#E7C595", #18 PPBP+
   "#D0D0D0", #mix
   "#EC5ECE" #20 undefined
 )
 
DimPlot(seurat_c1_full_FG, group.by="new.order", label=TRUE)+ scale_color_manual(values=umap.colors) + NoLegend() +ggtitle("UMAP")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))
 DimPlot(seurat_c1_full_FG, group.by="cluster_labels_res.0.4", label=T, repel = T)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("UMAP")


#Dimplot divided by group
colorControl <- c("#E2E2E2","#E2E2E2","#eaac43")
colorMild <- c("#E2E2E2","#E2E2E2","#38c3e1")
colorCritical <- c("#E2E2E2","#E2E2E2","#1657e1")

DimPlot(seurat_c1_full_FG,group.by="group_per_sample",cols=colorControl,order=c("control"))+NoLegend()+ggtitle("Control")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))
DimPlot(seurat_c1_full_FG,group.by="group_per_sample",cols=colorMild,order=c("mild"))+NoLegend()+ggtitle("Mild")+
theme(axis.text.x = element_text(size = 17),
      axis.text.y = element_text(size = 17),  
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      title =element_text(size=17, face='bold'))
DimPlot(seurat_c1_full_FG,group.by="group_per_sample",cols=colorCritical,order=c("severe"))+NoLegend()+ggtitle("Severe")+
theme(axis.text.x = element_text(size = 17),
      axis.text.y = element_text(size = 17),  
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      title =element_text(size=17, face='bold'))

#Cluster Marker Genes
Idents(seurat_c1_full_FG) <- "new.order"
markers <- FindAllMarkers(object = seurat_c1_full_FG)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

features=c("CCL2","CXCL8","CD83","CCL4","TNF",
           "PTGS2","IL6","GYPB","AREG","AZU1","MCEMP1",
           "LCN2","CTSG","CEACAM8","ELANE","DEFA4","MMP8",
           "GYPA","OLFM4","CD177","ARG1")

library(ggplot2)
#Run it and order
table(seurat_c1_full_FG$cluster_labels_res.0.4)
seurat_c1_full_FG$cluster_labels_res.0.4 <- factor(seurat_c1_full_FG$cluster_labels_res.0.4, 
                               levels =  c("CD4+ T","Classical Monocytes", "CD8+ T", "NK", 
                                           "B", "HLA-DR- S100A+ Monocytes", "HLA-DR+ CD83+ Monocytes",
                                           "Non-classical Monocytes", "Neutrophils", "CD163+ Monocytes",  
                                           "Immature Neutrophils", "Megakaryocyte", "Plasmablasts",
                                           "mDCs", "pDCs", "mixed", "undefined"))

#Dotplot
DotPlot(seurat_c1_full_FG,features = features,
  split.by = "group_per_sample",cols = c(mild = "#38c3e1", severe= "#1657e1", control="#eaac43"))+
  coord_flip()+ theme(axis.title = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  

#Heatmap
cols.use <- list(group_per_sample=c("#eaac43","#38c3e1","#1657e1"))

DoMultiBarHeatmap(seurat_c1_full_FG, features=features, group.by = "cluster_labels_res.0.4", additional.group.by = "group_per_sample", additional.group.sort.by = "group_per_sample", size = 0, cols.use = cols.use, draw.lines = TRUE, lines.width = NULL)+
  theme(text = element_text(size = 17, color = "black"))+
  scale_fill_gradientn(colours=pals::parula(25), na.value = "white")

   
#Feature plot
FeaturePlot(seurat_c1_full_FG, features=c("CCL2","IFNG","CXCL8","IRG1","CD83","CCL4","TNF",
                                          "SDC2","PTGS2","IL6","GYPB","AREG","AZU1","MCEMP1",
                                          "LCN2","IL1R1","CTSG","CEACAM8","ELANE","DEFA4","MMP8",
                                          "GYPA","OLFM4","CD177","ARG1"), order = T, ncol = 5)  

#ViolinPlot
VlnPlot(seurat_c1_full_FG, group.by="cluster_labels_res.0.4", features=c("CTSG", "ELANE", "ARG1", "IL1R1"), ncol = 2)

##### DIVIDED BY GROUPS #####

##CONTROL
#subset
seurat_c1_casec_FG <- subset(seurat_c1_full_FG, subset=group_per_sample %in% c("control"))
seurat_c1_ctrl_FG <- subset(seurat_c1_casec_FG)
remove(seurat_c1_casec_FG)

#MILD
#subset
seurat_c1_casem_FG <- subset(seurat_c1_full_FG, subset=group_per_sample %in% c("mild"))
seurat_c1_mild_FG <- subset(seurat_c1_casem_FG)
remove(seurat_c1_casem_FG)

#SEVERE
#subset
seurat_c1_cases_FG <- subset(seurat_c1_full_FG, subset=group_per_sample %in% c("severe"))
seurat_c1_severe_FG <- subset(seurat_c1_cases_FG)
remove(seurat_c1_cases_FG)

#FeaturePlot
FeaturePlot(seurat_c1_ctrl_FG, features=features, 
            order = T, ncol = 1, cols = c('lightgrey','#009A49'))

FeaturePlot(seurat_c1_mild_FG, features=features, 
            order = T, ncol = 1, cols = c('lightgrey','#F7B64F'))  

FeaturePlot(seurat_c1_severe_FG, feature=features, 
            order = T, ncol = 1, cols = c('lightgrey','#E33630')) 

#Dimplot
DimPlot(seurat_c1_ctrl_FG, group.by="new.order", label=TRUE, label.size = 5.5)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("Control")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))
DimPlot(seurat_c1_mild_FG, group.by="new.order", label=TRUE, label.size = 5.5)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("Mild")+ 
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))
DimPlot(seurat_c1_severe_FG, group.by="new.order", label=TRUE, label.size = 5.5)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("Severe")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))


#UMAP for genes
#AREG
FeaturePlot(candida_Stim, features=c("AREG"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)

FeaturePlot(candida_Unstim, features=c("AREG"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#ARG1
FeaturePlot(seurat_c1_ctrl_FG, features=c("ARG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("ARG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("ARG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#AZU1
FeaturePlot(seurat_c1_ctrl_FG, features=c("AZU1"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("AZU1"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("AZU1"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#CCL2
FeaturePlot(seurat_c1_ctrl_FG, features=c("CCL2"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CCL2"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CCL2"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#CCL4
FeaturePlot(seurat_c1_ctrl_FG, features=c("CCL4"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,5),limits=c(0,5), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CCL4"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,5),limits=c(0,5), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CCL4"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,5),limits=c(0,5), oob=squish)


#CD177
FeaturePlot(seurat_c1_ctrl_FG, features=c("CD177"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CD177"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CD177"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#CD83
FeaturePlot(seurat_c1_ctrl_FG, features=c("CD83"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CD83"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CD83"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#CEACAM8
FeaturePlot(seurat_c1_ctrl_FG, features=c("CEACAM8"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CEACAM8"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CEACAM8"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#CTSG
FeaturePlot(seurat_c1_ctrl_FG, features=c("CTSG"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CTSG"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CTSG"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#CXCL8
FeaturePlot(seurat_c1_ctrl_FG, features=c("CXCL8"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ NoLegend()+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))+
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,4),limits=c(0,4), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("CXCL8"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ NoLegend()+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))+
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,4),limits=c(0,4), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("CXCL8"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ NoLegend()+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))+
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,4),limits=c(0,4), oob=squish)

#DEFA4
FeaturePlot(seurat_c1_ctrl_FG, features=c("DEFA4"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,4),limits=c(0,4), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("DEFA4"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,4),limits=c(0,4), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("DEFA4"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,4),limits=c(0,4), oob=squish)

#ELANE
FeaturePlot(seurat_c1_ctrl_FG, features=c("ELANE"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("ELANE"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("ELANE"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#GYPA
FeaturePlot(seurat_c1_ctrl_FG, features=c("GYPA"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("GYPA"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("GYPA"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#GYPB
FeaturePlot(seurat_c1_ctrl_FG, features=c("GYPB"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("GYPB"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("GYPB"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#IFNG
FeaturePlot(seurat_c1_ctrl_FG, features=c("IFNG"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("IFNG"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("IFNG"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#IL1R1
FeaturePlot(seurat_c1_ctrl_FG, features=c("IL1R1"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("IL1R1"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("IL1R1"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#IL6
FeaturePlot(seurat_c1_ctrl_FG, features=c("IL6"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("IL6"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("IL6"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#IRG1
FeaturePlot(seurat_c1_ctrl_FG, features=c("IRG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("IRG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("IRG1"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,1),limits=c(0,1), oob=squish)

#LCN2
FeaturePlot(seurat_c1_ctrl_FG, features=c("LCN2"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("LCN2"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("LCN2"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#MCEMP1
FeaturePlot(seurat_c1_ctrl_FG, features=c("MCEMP1"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("MCEMP1"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("MCEMP1"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#MMP8
FeaturePlot(seurat_c1_ctrl_FG, features=c("MMP8"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("MMP8"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,2),limits=c(0,2), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("MMP8"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,2),limits=c(0,2), oob=squish)

#OLFM4
FeaturePlot(seurat_c1_ctrl_FG, features=c("OLFM4"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("OLFM4"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("OLFM4"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,1),limits=c(0,1), oob=squish)

#PTGS2
FeaturePlot(seurat_c1_ctrl_FG, features=c("PTGS2"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("PTGS2"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("PTGS2"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)

#SDC2
FeaturePlot(seurat_c1_ctrl_FG, features=c("SDC2"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("SDC2"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,1),limits=c(0,1), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("SDC2"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,1),limits=c(0,1), oob=squish)

#TNF
FeaturePlot(seurat_c1_ctrl_FG, features=c("TNF"), 
            order = T, ncol = 1, cols = c('lightgrey','#eaac43'))+ 
  scale_color_gradient(low='lightgrey',high = '#eaac43',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_mild_FG, features=c("TNF"), 
            order = T, ncol = 1, cols = c('lightgrey','#38c3e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#38c3e1',breaks=c(0,3),limits=c(0,3), oob=squish)
FeaturePlot(seurat_c1_severe_FG, features=c("TNF"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1',breaks=c(0,3),limits=c(0,3), oob=squish)


#############################
#This can be used after having a seurat object and transformed data

#BOXPLOT OF CELL TYPE PERCENTAGE BY DONOR AND CONDITION

library(plyr)
library(dplyr)

a <- seurat_c1_full_FG@meta.data
a <- ddply(a, .(donor), mutate, n_donor = length(donor))
a <- ddply(a, .(donor, cluster_labels_res.0.4), mutate, n_cell = length(donor))
df <- mutate(a, percent_cell= n_cell/n_donor) #df means dataframe


my_comparisons <- list( c("control", "mild"), c("mild", "severe"), c("control", "severe") )

ggplot(df, aes(x = group_per_sample, y = percent_cell, fill= group_per_sample))+
  facet_wrap(.~cluster_labels_res.0.4, scales = "free_y")+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.01, position = position_jitter(.1))+
  scale_fill_manual("Group", values = c(control="#eaac43", mild = "#38c3e1", severe= "#1657e1"))+
  labs(x = "Group", y= "CXCL8 expression")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_line(),
        strip.background = element_rect(color="white", fill="white"),
        strip.placement = "outside",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(2,"line"))


#BOXPLOT: GENE EXPRESSION BY CELL TYPE

tabela <- seurat_c1_full_FG@assays$RNA@scale.data
genes <- c("CCL16",	"CCL2",	"CCL25",	"CCL4",	"CCL7",	"CD177",	"CX3CL1",	"CXCL10",	"CXCL11",	"CXCL13",	"CXCL3",	"CXCL5",	"CXCL6",	"CXCL9",	"DYSF",	"EDN1",	"IL1F10",	"ITGA1",	"ITGA9",	"JAM3",	"MDK",	"PPBP",	"PRTN3",	"S100A12",	"S100A8",	"S100A9",	"SAA1",	"SLAMF8",	"SLIT2",	"XG")
tabela <- tabela[genes,]
tabela <- t(tabela)
tabela<-as.data.frame(tabela)

#Build an object that is composed by $ (coluna com título CCL16) da tabela
CCL16<- tabela$CCL16
CCL2<- tabela$CCL2
CCL25<- tabela$CCL25
CCL4<- tabela$CCL4
CCL7<- tabela$CCL7
CD177<- tabela$CD177
CX3CL1<- tabela$CX3CL1
CXCL10<- tabela$CXCL10
CXCL11<- tabela$CXCL11
CXCL13<- tabela$CXCL13
CXCL3<- tabela$CXCL3
CXCL5<- tabela$CXCL5
CXCL6<- tabela$CXCL6
CXCL9<- tabela$CXCL9
DYSF<- tabela$DYSF
EDN1<- tabela$EDN1
IL1F10<- tabela$IL1F10
ITGA1<- tabela$ITGA1
ITGA9<- tabela$ITGA9
JAM3<- tabela$JAM3
MDK<- tabela$MDK
PPBP<- tabela$PPBP
PRTN3<- tabela$PRTN3
S100A12<- tabela$S100A12
S100A8<- tabela$S100A8
S100A9<- tabela$S100A9
SAA1<- tabela$SAA1
SLAMF8<- tabela$SLAMF8
SLIT2<- tabela$SLIT2
XG<- tabela$XG

a <- seurat_c1_full_FG@meta.data
a <- cbind(a, AZU1,MCEMP1,LCN2,CTSG,CEACAM8,ELANE,DEFA4,MMP8,OLFM4,CD177, ARG1, CXCL8)
df <- aggregate(a[,genes], list(a$donor, a$cluster_labels_res.0.4), mean)
columns <- c("group_per_sample", "donor", "cluster_labels_res.0.4", "percent_cell")
df1 <- unique(a[, columns])
df <- df[order(df$Group.1),]
df1 <- df1[order(df1$donor),]
group_per_sample <- df1$group_per_sample
df <- cbind(df, group_per_sample)

library(reshape2)
dfmelted <- melt(df)
dfmelted_IN<-subset(dfmelted, Group.2=="Immature Neutrophils")

my_comparisons <- list( c("control", "mild"), c("mild", "severe"), c("control", "severe") )

library(lemon)
library(ggpubr)

ggplot(dfmelted_IN, aes(x = group_per_sample, y = value, fill= group_per_sample))+
  facet_rep_wrap(.~variable, ncol=6)+
  geom_boxplot(outlier.shape = NA, lwd = 0.2)+
  ylim(-0.5, 1.1)+
  stat_compare_means(method= "t.test", label = "p.signif", hide.ns = T, label.y = c(0.6,0.7,0.8),
                     aes(group = group_per_sample), comparisons = my_comparisons, vjust = 0.5)+
  geom_jitter(size=0.6, position = position_jitter(.1), alpha=0.5)+
  scale_fill_manual("Group", values = c(control="#eaac43", mild = "#38c3e1", alpha(c(severe= "#1657e1"),0.5)))+
  labs(x = "Genes", y= "Average Expression")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line.x=element_line(),
        strip.background = element_rect(color="white", fill="white"),
        strip.placement = "outside",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5,"line"),
        legend.background = element_blank(),
        axis.line.y = element_line(),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text  = element_text(size = 11, color = "black"))+
  scale_y_continuous(breaks=seq(-0.4, 1, 0.2), limits=c(-0.5, 1.1))+NoLegend()




###############################
#Divided by gender

#MALE SEVERE
seurat_c1_casema_FG <- subset(seurat_c1_severe_FG, subset=sex %in% c("male"))
seurat_c1_male_FG <- subset(seurat_c1_casema_FG)
remove(seurat_c1_casema_FG)


#FEMALE SEVERE
seurat_c1_casefe_FG <- subset(seurat_c1_severe_FG, subset=sex %in% c("female"))
seurat_c1_female_FG <- subset(seurat_c1_casefe_FG)
remove(seurat_c1_casefe_FG)


DimPlot(seurat_c1_severe_FG,group.by="sex", split.by = "sex")
DimPlot(seurat_c1_male_FG,group.by="cluster_labels_res.0.4", label=T)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("Male")
DimPlot(seurat_c1_female_FG,group.by="cluster_labels_res.0.4", label=T)+ scale_color_manual(values=umap.colors) + NoLegend()+ggtitle("Female")


FeaturePlot(seurat_c1_severe_FG, features=c("CTSG", "ELANE", "ARG1", "IL1R1"), 
            order = F, ncol = 2, cols = c('lightgrey','#E33630'), split.by = "sex") 

FeaturePlot(seurat_c1_severe_FG, features=c("CCL2","IFNG","CXCL8","CD83","CCL4","TNF"), 
            order = F, ncol = 2, cols = c('lightgrey','#E33630'), split.by = "sex")

FeaturePlot(seurat_c1_severe_FG, features=c("SDC2","PTGS2","IL6","GYPB","AREG","AZU1","MCEMP1"), 
            order = F, ncol = 2, cols = c('lightgrey','#E33630'), split.by = "sex")

FeaturePlot(seurat_c1_severe_FG, features=c("LCN2","CEACAM8","DEFA4","MMP8", "GYPA","OLFM4","CD177"), 
            order = F, ncol = 2, cols = c('lightgrey','#E33630'), split.by = "sex")


VlnPlot(seurat_c1_severe_FG, group.by="sex", features=c("CTSG", "ELANE", "ARG1", "IL1R1"), ncol = 4)
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("CTSG", "ELANE", "ARG1", "IL1R1"), ncol = 4)
VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("CTSG", "ELANE", "ARG1", "IL1R1"), ncol = 4)

#Doplot splited
DotPlot(seurat_c1_male_FG,features = features, group.by="cluster_labels_res.0.4")+
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(n =11, name = "GnBu"))+
  theme(axis.title = element_blank())+ coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(seurat_c1_female_FG,features = features, group.by="cluster_labels_res.0.4")+
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(n =11, name = "GnBu"))+
  theme(axis.title = element_blank())+ coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Dotplot male and female together
DotPlot(seurat_c1_full_FG,features = features, group.by="cluster_labels_res.0.4", split.by = "sex", 
        cols = c("royalblue", "greeen"))+
  theme(axis.title = element_blank())+ coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Violinplots
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("CCL2","IFNG","CXCL8","CD83"), ncol=2)
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("SDC2","PTGS2","IL6","GYPB"), ncol=2)
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("AZU1","MCEMP1", "TNF","OLFM4"), ncol=2)
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("LCN2","CEACAM8","DEFA4","MMP8"), ncol=2)
VlnPlot(seurat_c1_male_FG, group.by="cluster_labels_res.0.4", features=c("GYPA","CD177","AREG","CCL4"), ncol=2)

VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("CCL2","IFNG","CXCL8","CD83"), ncol=2)
VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("SDC2","PTGS2","IL6","GYPB"), ncol=2)
VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("AZU1","MCEMP1", "TNF","OLFM4"), ncol=2)
VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("LCN2","CEACAM8","DEFA4","MMP8"), ncol=2)
VlnPlot(seurat_c1_female_FG, group.by="cluster_labels_res.0.4", features=c("GYPA","CD177","AREG","CCL4"), ncol=2)
