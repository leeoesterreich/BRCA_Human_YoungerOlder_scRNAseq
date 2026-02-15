---
title: "AgeStudy_Wuetal_scRNAseq"
format: html
editor: visual
R version: 4.5.2 (2025-10-31) -- "Not Part in a Rumble"
---

# Human scRNA-seq data analysis information, dataset from Wu et al. 2021. Nat. Genetics

## scRNA-seq data available:

Human single cell RNA sequencing (scRNA-seq) of ER+ breast cancer, publically available at Wu et al. 2021. Nat. Genetics

## Preparation1. Install and load the necessary R packages.

```{r}
library(Seurat)
library(ggplot2)
library(tidyr)
library(data.table)
library(cowplot)  # to use plot_grid()
library(ggrastr) # CRAN
library(dittoSeq)  # bioconductor - it requires "nloptr" CRAN package     
library(glmGamPoi)
#library(DoubletFinder) 
library(patchwork)  # CRAN. to use plot_annotation()
library(harmony)  #  install.packages("harmony") to use RunHarmony() https://github.com/immunogenomics/harmony
library(reticulate)
library(gridExtra)  # for 'grid.arrange' fun

```

## Preparation2. Set working directory and direct your input data files

```{r}
dir<-dirname(rstudioapi::getSourceEditorContext()$path); print(dir);  
setwd(dir);  
SeuratObj_GSE176 <- "SeuratObj_GSE176078_ERpos_NewMeta_AfterQCSCT.rds"  
# This data file was made while publishing Lee et al. Oncoimmunology 2024. 
# This file size is 1.7Gb. If you need this file, contact the corresponding author.
ClinicalDataFile <- "ClinicalData_Wu_scRNAseq_26p.txt"  # I downloaded clinical data file from Wu et al. and processed it. I posted this file in this Github

```

# Section I. Read scRNA-seq data and metadata, and process them.

## Step1. Read clinical data and seurat object RDS data

```{r}
## Clinical data. 
ClinicalData <- fread(ClinicalDataFile, header=TRUE, stringsAsFactors = FALSE);
colnames(ClinicalData) <- gsub("\\ ", "", colnames(ClinicalData));  dim(ClinicalData); head(ClinicalData) # 26 7
ClinicalData$AgeGroup <- ifelse(ClinicalData$Age > 80, "Elderly",  ifelse(ClinicalData$Age <= 50, "Young", "MidAge"))

## Seurat object data and run PCA and UMAP
SeuratObject_GSE176 <- readRDS(SeuratObj_GSE176)
SeuratObject_GSE176_PCA <- RunPCA(SeuratObject_GSE176, npcs = 30, verbose = F)
SeuratObject_GSE176_PCA <- RunUMAP(SeuratObject_GSE176_PCA, reduction = "pca", dims = 1:30, verbose = F)

# Visualization by UMAP by CaseID - Before Harmony
table(SeuratObject_GSE176_PCA$orig.ident); length(table(SeuratObject_GSE176_PCA$orig.ident)); #  26 patients total but only ER+ 14 patients. 

## plot_annotation() requires 'patchwork' R package
MyDimplot_SeuratMerge_UMAP <- DimPlot(SeuratObject_GSE176_PCA,reduction="umap") + plot_annotation(title="Before Harmony integration")
#ggsave(MyDimplot_SeuratMerge_UMAP, height=8,width=8, dpi=300, filename=paste0("OutUMAPplot_GSE176078_MergeSCT_PreTreatBfHarmony.pdf"), useDingbats=FALSE)

```

## Step2. Subset only Young (n=3), MidAge (n=4) and Elderly patients (n=3)

```{r}
### Subset by orig.ident - I need only Young (n=3), MidAge (n=4), and Elderly (n=3) samples 
Seurat_GSE176078_YoungElderly <- subset(x=SeuratObject_GSE176_PCA, subset=orig.ident %in%  c("CID3941","CID4530N","CID4535", "CID3948","CID4067","CID4290A", "CID4463","CID4040","CID4471","CID4461")) # 18,832 cells
table(Seurat_GSE176078_YoungElderly$orig.ident) # 10 patients
## Add "Young" and "Elderly" labels 
Seurat_GSE176078_YoungElderly@meta.data$AgeGroup <- ifelse(Seurat_GSE176078_YoungElderly@meta.data$orig.ident %in% c("CID3941","CID4530N","CID4535"), "Young",
                           ifelse(Seurat_GSE176078_YoungElderly@meta.data$orig.ident %in% c("CID3948","CID4067","CID4290A"), "Elderly",
                           ifelse(Seurat_GSE176078_YoungElderly@meta.data$orig.ident %in% c("CID4463","CID4040","CID4471","CID4461"), "MidAge", "")))

table(Seurat_GSE176078_YoungElderly@meta.data$orig.ident, Seurat_GSE176078_YoungElderly@meta.data$AgeGroup)
Seurat_GSE176078_YoungElderly <- SetIdent(Seurat_GSE176078_YoungElderly, value=Seurat_GSE176078_YoungElderly$orig.ident)   # change active.ident
Seurat_GSE176078_YoungElderly@meta.data[1:2,]; table(Seurat_GSE176078_YoungElderly$AgeGroup) # Elderly: 10713 # MidAge: 12127  # Young 8119 #

saveRDS(Seurat_GSE176078_YoungElderly, file="Seurat_GSE176078_PreTreat_Young3sMidAge4sElderly3s.rds")  
```

## Step3. UMAP Reduction by Harmony and clustering

```{r}
table(Seurat_GSE176078_YoungElderly$CellTypeSubset)
table(Seurat_GSE176078_YoungElderly$CellTypeMajor, Seurat_GSE176078_YoungElderly$CellTypeMinor)

Seurat_GSE176078_YoungElderly$CellTypeMajor <- gsub("-","",Seurat_GSE176078_YoungElderly$CellTypeMajor)
table(Seurat_GSE176078_YoungElderly$CellTypeMajor)

## Metadata 
Metadata_GSE176078 <- Seurat_GSE176078_YoungElderly@meta.data %>% data.frame
Metadata_GSE176078$CellTypeAnnotSH <- ifelse(Metadata_GSE176078$CellTypeMajor %in% c("Myeloid","Tcells"), Metadata_GSE176078$CellTypeMinor, Metadata_GSE176078$CellTypeMajor)

### Remove special characters in cell type
Metadata_GSE176078$CellTypeAnnotSH <- gsub("-|_|\\+","",Metadata_GSE176078$CellTypeAnnotSH)
Seurat_GSE176078_YoungElderly$CellTypeAnnotSH <- Metadata_GSE176078$CellTypeAnnotSH
table(Seurat_GSE176078_YoungElderly@meta.data$CellTypeAnnotSH )

saveRDS(Seurat_GSE176078_YoungElderly, file="Seurat_GSE176078_PreTreat_Young3sMidAge4sElderly3s_CellTypeAnnotSH.rds")  

MedtaData <- Seurat_GSE176078_YoungElderly@meta.data[, c("CaseID","CellTypeAnnotSH")]; dim(MedtaData); MedtaData[1:2,] # 30959     3
MedtaData_CellID <- MedtaData %>% tibble::rownames_to_column("CellID")
fwrite(MedtaData_CellID, file="Metadata_GSE176078_YoungElderly_30959c.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

####### =========== inner_join between Metadata and Clinical data.  ========== #########
Metadata_GSE176078_Proc <- Metadata_GSE176078 %>% tibble::rownames_to_column("CellID")
ClinicalData_Metadata <- dplyr::inner_join(ClinicalData,Metadata_GSE176078_Proc ); dim(ClinicalData_Metadata); ClinicalData_Metadata[1:2,]

### UMAP
MyDimplot_SeuratObject_GSE176_UMAP <- DimPlot(Seurat_GSE176078_YoungElderly,reduction="umap", label=TRUE, label.size=6, group.by="CellTypeAnnotSH") + plot_annotation(title="Two datasets UMAP_Harmony and Clustering") +
           theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratObject_GSE176_UMAP, height=8,width=11, dpi=300, filename=paste0("OutUMAPplot_GSE176078_CellTypeAnnotSH.pdf"), useDingbats=FALSE)

### UMAP by patient ID
MyDimplot_SeuratObject_OrigIdent_UMAP <- DimPlot(Seurat_GSE176078_YoungElderly,reduction="umap", label=TRUE, label.size=6, group.by="orig.ident") + plot_annotation(title="Two datasets UMAP_Harmony and Clustering") +
          theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratObject_OrigIdent_UMAP, height=8,width=11, dpi=300, filename=paste0("OutUMAPplot_GSE176078_ByOrigIdent.pdf"), useDingbats=FALSE)

### =========== Subset CancerEpithelial or Macrophage and label cells by AgeGroup ============== ##
Seurat_GSE176078_YoungElderly_CancerEpi <- subset(Seurat_GSE176078_YoungElderly, subset=CellTypeAnnotSH=="CancerEpithelial"); table(Seurat_GSE176078_YoungElderly_CancerEpi$CellTypeAnnotSH)
MyDimplot_UMAP_CanEpi_ByAgeGroup <- DimPlot(Seurat_GSE176078_YoungElderly_CancerEpi,reduction="umap", label=TRUE, label.size=6, group.by="AgeGroup") + plot_annotation(title="datasets UMAP_Harmony and Clustering") +
                  theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))

## Macrophage + Monocyte
Seurat_GSE176078_YoungElderly_MonoMacro <- subset(Seurat_GSE176078_YoungElderly, subset=CellTypeAnnotSH %in% c("Monocyte", "Macrophage")); table(Seurat_GSE176078_YoungElderly_Macrophage$CellTypeAnnotSH)
MyDimplot_UMAP_MonoMacro_ByAgeGroup <- DimPlot(Seurat_GSE176078_YoungElderly_MonoMacro,reduction="umap", label=TRUE, label.size=6, group.by="AgeGroup") + plot_annotation(title="datasets UMAP_Harmony and Clustering") +
            theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))

# Seurat_GSE176078_YoungElderly_Macrophage_UMAP <- Seurat_GSE176078_YoungElderly_MonoMacro %>% RunUMAP(reduction="umap", verbose=F, dims=1:2) %>% FindNeighbors(reduction="umap", k.param=15, dims=1:2)
Seurat_GSE176078_YoungElderly_Macrophage_PCA <- RunPCA(Seurat_GSE176078_YoungElderly_MonoMacro, features=VariableFeatures(object = Seurat_GSE176078_YoungElderly_MonoMacro))
MyDimplot_UMAP_Macrophage_ByAgeGroup <- DimPlot(Seurat_GSE176078_YoungElderly_Macrophage_PCA,reduction="pca", label=TRUE, label.size=6, group.by="AgeGroup") + plot_annotation(title="datasets UMAP_Harmony and Clustering") +
            theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), axis.text.y=element_text(vjust=0.6, size=25,angle=0))

Seurat_GSE176078_YoungElderlyTSNE <- RunTSNE(Seurat_GSE176078_YoungElderly_Macrophage_UMAP, dims = 1:10)
MyDimplot_UMAP_Macrophage_ByAgeGroup_TSNE <- DimPlot(Seurat_GSE176078_YoungElderlyTSNE, label=TRUE, label.size=6, group.by="AgeGroup", reduction="tsne")


```

## Step4. Add Age data to metadata

```{r}
MyMetaData <- Seurat_GSE176078_YoungElderly@meta.data; dim(MyMetaData); table(MyMetaData$orig.ident) # 30959    12
# CID3586  CID3838  CID3921  CID3941  CID3946  CID3948  CID3963  CID4040  CID4066  CID4067 CID4290A  CID4398 CID44041  CID4461  CID4463  CID4465  CID4471  CID4495 CID44971 CID44991 
# 0        0        0      585        0     2289        0     2417        0     3370     5054        0        0      476     1086        0     8148        0        0        0 
# CID4513  CID4515 CID45171  CID4523 CID4530N  CID4535 
# 0        0        0        0     4063     3471 

MyMetadataAge <- dplyr::inner_join(MyMetaData, ClinicalData[, c("CaseID","Age")]); MyMetadataAge[1:2,]; table(MyMetadataAge$CaseID, MyMetadataAge$Age)
Seurat_GSE176078_YoungElderly_Age <- AddMetaData(object=Seurat_GSE176078_YoungElderly, metadata=c(MyMetadataAge$Age), col.name=c("Age")); table(Seurat_GSE176078_YoungElderly_Age$Age)
table(Seurat_GSE176078_YoungElderly_Age$CaseID, Seurat_GSE176078_YoungElderly_Age$Age)

## Violin plot by CaseID, but CaseID should be ordered by age. 
## Remove cells of patients who don't have age data,
origIdent_Age <- Seurat_GSE176078_YoungElderly_Age@meta.data[,c("CaseID","Age")]
origIdent_Age_sort <- origIdent_Age[with(origIdent_Age, order(origIdent_Age$Age, decreasing=FALSE)), ]; head(origIdent_Age_sort); tail(origIdent_Age_sort)
Seurat_GSE176078_YoungElderly_Age$CaseID <- factor(Seurat_GSE176078_YoungElderly_Age$CaseID, levels=unique(origIdent_Age_sort$CaseID))

### Violin plot of intersting gene expression in each cell type. 
MyViolinplot_ByCelltype <- VlnPlot(object=Seurat_GSE176078_YoungElderly_Age, features=c("SHH","KLF4","GLI1","KIF7","SMO"), pt.size=0.00001, group.by="CellTypeAnnotSH", ncol=3)
ggsave(MyViolinplot_ByCelltype, height=9,width=12, dpi=300, filename=paste0("OutVlnplot_InterestingGene_By16CellType_20240215.pdf"), useDingbats=FALSE)

```

## Step5. Making RNA assay "scale.data"

```{r}
DefaultAssay(Seurat_GSE176078_YoungElderly_Age) <- "RNA"
Seurat_GSE176078_YoungElderly_Age <- NormalizeData(Seurat_GSE176078_YoungElderly_Age)
# saveRDS(Seurat_object_MainRename_NKT_Rmv_Normal, file="Seurat_object_MainRename_NKT_Rmv_Normal.rds")

all.genes <- rownames(Seurat_GSE176078_YoungElderly_Age)
Seurat_GSE176078_YoungElderly_Age_Scale <- ScaleData(Seurat_GSE176078_YoungElderly_Age, features = all.genes) ### This requires a lot of memory. Go to CRC
dim(Seurat_GSE176078_YoungElderly_Age_Scale[["RNA"]]@data) # 29377g 30959c

saveRDS(Seurat_GSE176078_YoungElderly_Age_Scale, file="Seurat_GSE176078_YoungElderly_Age_Scale.rds")

MyViolinplot_MainSub <- VlnPlot(object=Seurat_GSE176078_YoungElderly_Age_Scale, features=c("CTLA4","PDCD1","PDCD1LG2","CD274","LAG3","HAVCR2"), group.by="CaseID", stack=FALSE,flip=TRUE )
ggsave(MyViolinplot_MainSub, height=9,width=12, dpi=300, filename=paste0("OutViolin_ImmuneCheckpoint_ByPatient_SortYoungMidElderly_20231218.pdf"), useDingbats=FALSE)

##### Subset Macrophage only and make violin plot to show immunecheckpoint gene expression per patient. 
Seurat_GSE176078_YoungElderly_Age_Scale_MacSubset <- subset(Seurat_GSE176078_YoungElderly_Age_Scale, subset=CellTypeAnnotSH=="Macrophage"); table(Seurat_GSE176078_YoungElderly_Age_Scale_MacSubset$CellTypeAnnotSH)
MyViolinplot_MacroSubset <- VlnPlot(object=Seurat_GSE176078_YoungElderly_Age_Scale_MacSubset, features=c("CTLA4","PDCD1","PDCD1LG2","CD274","LAG3","HAVCR2"), group.by="CaseID", stack=FALSE,flip=TRUE )
ggsave(MyViolinplot_MainSub, height=9,width=12, dpi=300, filename=paste0("OutViolin_ImmuneCheckpoint_ByPatient_SortYoungMidElderly_Macrophage_20231218.pdf"), useDingbats=FALSE)
```

# Section 2. Visualize scRNA-seq data analysis with cell type fraction bar plot or gene expression violin plots in different age group.

## Step6. Make feature plots and cell type fraction bar plot

```{r}
## FeaturePlot 
MyFeaturePlot <- FeaturePlot(Seurat_GSE176078_YoungElderly,raster=FALSE, pt.size=1,features = c("EPCAM","KRT19","KRT18","FCGR3A", "CD3D","CD68","MS4A1", "PECAM1","IL7R",'CD8A'))   #  raster.dpi = c(512, 512)   
ggsave(MyFeaturePlot, height=13,width=18, dpi=300, filename=paste0("OutFeatureUMAP_ByMainMarker_HarmoneyGSE176EGA6608_Res1.5PC30KP15_10g.pdf"), useDingbats=FALSE)

## Cell type fraction per samples. Make dotplot for macrophage fraction. 
MyDittoPlot <- dittoBarPlot(object = Seurat_GSE176078_YoungElderly, var="CellTypeAnnotSH", group.by="orig.ident", x.reorder=c(1,9,10,  7,3,8,6,  2,4,5))  #, x.reorder=c(12,4,1,8,11,5,7,2,6,9,10,3))  #  x.reorder=c(1,4,5,8,11,12,2,3,6,7,9,10)
## https://github.com/satijalab/seurat/issues/962
ggsave(MyDittoPlot, height=5,width=8, dpi=300, filename=paste0("OutFractionBarplot_YoungMidElderly.pdf"), useDingbats=FALSE)

CellTypeFraction <- table(Seurat_GSE176078_YoungElderly$CellTypeAnnotSH, Seurat_GSE176078_YoungElderly@meta.data$CaseID)  ## This code does the same thing. 

## To make input data for cell type fraction correlation heatmap
CellTypeFraction_DF <- CellTypeFraction %>% matrix(ncol=ncol(CellTypeFraction)) %>% data.frame %>% t
colnames(CellTypeFraction_DF) <- rownames(CellTypeFraction)
rownames(CellTypeFraction_DF) <- colnames(CellTypeFraction)
head(CellTypeFraction_DF)
#CellTypeFraction_MacroDC <- CellTypeFraction_DF %>% dplyr::select(-c("MonocyteMacrophage")) # %>% dplyr::mutate(Macrophage = Mono_Macro + Macrophage) %>% dplyr::select(-c(Mono_Macro))
#fwrite(CellTypeFraction_MacroDC, file="CellTypeFraction_CellCountTcellSubtyp_GSE176EGA6608_TcellRes1.0.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

## To make input data for cell type fraction correlation heatmap  

ProportionCal <- CellTypeFraction_DF/rowSums(CellTypeFraction_DF);  dim(ProportionCal) #    10 16 

## Boxplot no filling in the box. Put dots for samples 
ProportionCal_Tibble <- ProportionCal %>% data.frame %>% tibble::rownames_to_column("CaseID") ; dim(ProportionCal_Tibble) # 10 17

## inner_join between ProportionCal and Clinicaldata
Clinicaldata_ProportionCal <- dplyr::inner_join( ClinicalData[,c("CaseID","Age","AgeGroup")],  ProportionCal_Tibble); 
dim(Clinicaldata_ProportionCal); Clinicaldata_ProportionCal[1:2,] # 10 19
fwrite(Clinicaldata_ProportionCal, file="ClinicalData_CelltypeProportion.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## Make boxplot input data 
CellTypeFraction_Tidyr <- Clinicaldata_ProportionCal %>% tidyr::gather( key="CellType", value="CellFraction", colnames(Clinicaldata_ProportionCal)[4]:colnames(Clinicaldata_ProportionCal)[ncol(Clinicaldata_ProportionCal)]); 
CellTypeFraction_Tidyr_Proc <- CellTypeFraction_Tidyr %>% dplyr::mutate(CellType_AgeGroup = paste0(CellType, "_", AgeGroup)) %>% dplyr::filter(AgeGroup != "MidAge")
head(CellTypeFraction_Tidyr_Proc)  # use column names for tidyr range.
table(CellTypeFraction_Tidyr_Proc$CaseID)

BoxPlot <-ggplot(CellTypeFraction_Tidyr_Proc, aes(x=CellType_AgeGroup, y=CellFraction))+ ##, fill=expValue, colour=expValue)) + 
  ##https://stackoverflow.com/questions/22808219/ggplot2-how-do-you-control-point-color-with-geom-jitter-and-geom-point
  stat_boxplot(geom='errorbar', width=0.5, lwd=0.5) + geom_boxplot(color="black", varwidth=FALSE, lwd=0.5)+  # geom_boxplot(fill="white")  #lwd=1.2  To make box lines thicker
  #If I put geom_boxplot + stat_boxplot, it makes vertical line inside the box. 
  scale_y_continuous(expand=c(0,0), limits=c(-0.05,1.0))+     # TP63 limits=c(-20,55))
  theme_classic()
BoxPlot

MyPlot_Dot <- BoxPlot + geom_point(aes(fill=AgeGroup), pch=21,size=2, position=position_jitterdodge()) + 
  scale_fill_manual(values=rep(c("#009900", "#E69F00"),9))+ 
  theme(legend.position="none") +   # or pch=21  
  theme(axis.text.x=element_text(colour="Black",size=10,angle=90,hjust=.5,vjust=.5,face="bold"),  #face="plain" 
        axis.text.y=element_text(colour="Black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
        axis.title.y=element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"), 
        axis.line.y=element_line(color="black", linewidth=1.0),
        axis.line.x=element_line(color="black", linewidth =1.0),
        axis.ticks=element_line(colour="black",linewidth=1.0),
        axis.ticks.length=unit(.22, "cm")) 
MyPlot_Dot    
ggsave(MyPlot_Dot, height=6,width=10, dpi=300, filename=paste0("OutFractionBoxDotplot_YoungEldery.pdf"), useDingbats=FALSE)

```

## Step7. Violin plot for interesting genes across 10 smp based on cell type fraction in each cell type

```{r}
Clinicaldata_ProportionCal_Sort <- Clinicaldata_ProportionCal[with(Clinicaldata_ProportionCal, order(Age, decreasing=FALSE)), ]; Clinicaldata_ProportionCal_Sort[,1:3]

CaseID_SortByMacroFraction <- Clinicaldata_ProportionCal_Sort$CaseID; print(CaseID_SortByMacroFraction)  # [1] "CID4530N" "CID4535"  "CID3941"  "CID4461"  "CID4471"  "CID4040"  "CID4463"  "CID3948"  "CID4067"  "CID4290A"

SeuratObj_ERpos_CellTypeMacrophage <- Seurat_GSE176078_YoungElderly
SeuratObj_ERpos_CellTypeMacrophage <- SetIdent(SeuratObj_ERpos_CellTypeMacrophage, value=SeuratObj_ERpos_CellTypeMacrophage$CellTypeAnnotSH)   
table(SeuratObj_ERpos_CellTypeMacrophage@active.ident)


saveRDS(SeuratObj_ERpos_CellTypeMacrophage, file="SeuratObj_ERpos_CellTypeMacrophage_AgeGroup.rds")

## subset by cell type
MyCellType <- names(table(SeuratObj_ERpos_CellTypeMacrophage$CellTypeAnnotSH))
MyCellType_NoNA <- MyCellType; print(MyCellType_NoNA)  # 16 cell types   # MyCellType[MyCellType!="NA"]; 
# [1] "Bcells"           "CAFs"             "CancerEpithelial" "CyclingMyeloid"   "CyclingTcells"    "DCs"              "Endothelial"      "Macrophage"       "Monocyte"         "NKcells"          "NKTcells"         "NormalEpithelial"
# [13] "Plasmablasts"     "PVL"              "TcellsCD4"        "TcellsCD8"

MyCellType_NoNASelc <- MyCellType_NoNA[c(1,2,3,4,5,8,15,16)]; print(MyCellType_NoNASelc)
#[1] "Bcells"           "CancerEpithelial" "CyclingMyeloid"   "CyclingTcells"    "Macrophage"       "TcellsCD4"        "TcellsCD8"      

Feature_CanEpi <- c("SAA1","HMGA1P7","RPL12P37","PAH","TRPM8","JCHAIN","DMRTC1","OR8A1","KRT37","PMS2P6",
                    "SLC6A3","HAAO","RAB19","VIL1","RP1","ACP7","SLX1B","KRT20","XKR7","CDHR2"); length(Feature_CanEpi) # 20
Feature_CanEpi_v2 <- c("USP14","CCT8","CPR180","KIF20A","EIF2AK1","EXTL2","CLASP2","BRWD3","ANLN","STIL","ZNF770","SETD7","BBS9","GTPBP4","CS","COPA"); length(Feature_CanEpi) # 16
Feature_Mac <- c("IL1A","IL1B","IL1RN","IL2","IL3","IL4","IL5","IL6","IL7","CXCL8","IL9","IL10","IL12B","IL13","IL15","IL17A","IL33","CD40LG","EGF",
                 "CCL11","FGF2","CSF3","CSF2","IFNA17","IFNG","CXCL1","CXCL2","CXCL3"); length(Feature_Mac) # 28
Feature_Mac_v2 <- c ("CXCL5","CXCL9","CXCL10","CXCL13","CCL2","CCL3","CCL4","CCL7","CCL8","CCL11","CCL21","CCL22","CCL24","CCL26","CCL27",
                     "TNF","TGFB1","TNFSF10","VEGFA","CD274","PDCD1LG2","CD163","MRC1","NOS2","HLA-A"); length(Feature_Mac_v2); #  25


DefaultAssay(SeuratObj_ERpos_CellTypeMacrophage) <- "RNA"
SeuratObj_ERpos_CellTypeMacrophage_Normal <- NormalizeData(SeuratObj_ERpos_CellTypeMacrophage)
# saveRDS(SeuratObj_ERpos_CellTypeMacrophage_Normal, file="SeuratObj_ERpos_CellTypeMacrophage_Normal.rds")

all.genes <- rownames(SeuratObj_ERpos_CellTypeMacrophage_Normal)
SeuratObj_ERpos_CellTypeMacrophage_Scale <- ScaleData(SeuratObj_ERpos_CellTypeMacrophage_Normal, features = all.genes) ### This requires a lot of memory. Go to CRC
saveRDS(SeuratObj_ERpos_CellTypeMacrophage_Scale, file="SeuratObj_ERpos_CellTypeMacrophage_Scale.rds")
# SeuratObj_ERpos_CellTypeMacrophage_Scale <- readRDS("SeuratObj_ERpos_CellTypeMacrophage_Scale.rds")


SeuratObj_ERpos_CellTypeMacrophage_Scale$AgeGroup <- factor(SeuratObj_ERpos_CellTypeMacrophage_Scale$AgeGroup, levels=c("Young","MidAge","Elderly"))
## Make violin plot per patient. 
for(EachCellType in MyCellType_NoNA) { # MyCellType_NoNASelc
    # EachCellType <- MyCellType_NoNASelc[6]; print(EachCellType)
    SeuratObj_ERpos_SubsetMyeloid <- subset(x=SeuratObj_ERpos_CellTypeMacrophage_Scale, idents=c(EachCellType)); table(SeuratObj_ERpos_SubsetMyeloid@active.ident)
    SeuratObj_ERpos_SubsetMyeloid$orig.ident <- factor(x=SeuratObj_ERpos_SubsetMyeloid$orig.ident, levels=CaseID_SortByMacroFraction)
    
    # MyFeature <- c("CCL2","CCR2","CSF1","CXCL9","SHH","HMGB1","AGER","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","HSD17B7")
    MyFeature <- c("CCR1","CCR3","CCR4","CD74","APP","CXCL10","DPP4","EGFR","AGER","GPR75","CCL5","IL2RG","IL15","S100A11")
     
    MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_SubsetMyeloid, features = c(MyFeature), group.by="AgeGroup", slot = 'data')  # raster=TRUE requires ggrastr r packge.  # group.by = "orig.ident"
    OutFile <- paste0("OutVlnplot_CellPhoneDBGene_", EachCellType,  "_ByAgeGroup_20231208.pdf"); print(OutFile)
    ggsave(MyViolinplot, height=11,width=12, dpi=300, filename=OutFile, useDingbats=FALSE)  #  height=26,width=19,
}
```

## Step8. Violin plot DefaultAssay RNA for "scale.data"

```{r}
DefaultAssay(SeuratObj_ERpos_CellTypeMacrophage) <- "RNA"
SeuratObj_ERpos_CellTypeMacrophage_Normal <- NormalizeData(SeuratObj_ERpos_CellTypeMacrophage)
saveRDS(SeuratObj_ERpos_CellTypeMacrophage_Normal, file="SeuratObj_ERpos_CellTypeMacrophage_Normal.rds")

all.genes <- rownames(Seurat_HarmonyUMAP_NewMainCellType_Normal)
# Seurat_HarmonyUMAP_NewMainCellType_Scale <- ScaleData(Seurat_HarmonyUMAP_NewMainCellType_Normal, features = all.genes) ### This requires a lot of memory. Go to CRC
SeuratObj_ERpos_CellTypeMacrophage_Scale <- readRDS("SeuratObj_ERpos_CellTypeMacrophage_Scale.rds")

## Make violin plot per AgeGroup 
for(EachCellType in MyCellType_NoNA) { # MyCellType_NoNASelc
  # EachCellType <- MyCellType_NoNASelc[3]; print(EachCellType)
  SeuratObj_ERpos_SubsetMyeloid <- subset(x=SeuratObj_ERpos_CellTypeMacrophage_Scale, idents=c(EachCellType)); table(SeuratObj_ERpos_SubsetMyeloid@active.ident)
  SeuratObj_ERpos_SubsetMyeloid$orig.ident <- factor(x=SeuratObj_ERpos_SubsetMyeloid$orig.ident, levels=CaseID_SortByMacroFraction)
  
  MyFeature <- c("CCL2","CCL3","CCL4","TNF","TGFB1","CD163","MRC1")
  
  MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_SubsetMyeloid, features = c(MyFeature), group.by="AgeGroup", slot = 'scale.data')  # raster=TRUE requires ggrastr r packge.  # group.by = "orig.ident"
  OutFile <- paste0("OutVlnplot_CCL2_TGFB1_CD163_", EachCellType,  "_ByAgeGroup_20231030.pdf"); print(OutFile)
  ggsave(MyViolinplot, height=11,width=12, dpi=300, filename=OutFile, useDingbats=FALSE)  #  height=26,width=19,
}


## Sort CaseID by age
SeuratObj_ERpos_CellTypeMacrophage_Scale$CaseID <- factor(SeuratObj_ERpos_CellTypeMacrophage_Scale$CaseID, levels=c(Clinicaldata_ProportionCal_Sort$CaseID))

## Make violin plot per each patients sorted by age. 
for(EachCellType in MyCellType_NoNA) { # MyCellType_NoNASelc
  # EachCellType <- MyCellType_NoNASelc[6]; print(EachCellType)
  SeuratObj_ERpos_SubsetMyeloid <- subset(x=SeuratObj_ERpos_CellTypeMacrophage_Scale, idents=c(EachCellType)); table(SeuratObj_ERpos_SubsetMyeloid@active.ident)
  SeuratObj_ERpos_SubsetMyeloid$orig.ident <- factor(x=SeuratObj_ERpos_SubsetMyeloid$orig.ident, levels=CaseID_SortByMacroFraction)
  
  MyFeature <- c("CCL2","CCL3","CCL4","TNF","TGFB1","CD163","MRC1")
  
  MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_SubsetMyeloid, features = c(MyFeature), group.by="CaseID", slot = 'scale.data')  # raster=TRUE requires ggrastr r packge.  # group.by = "orig.ident"
  OutFile <- paste0("OutVlnplot_CCL2_TGFB1_CD163_", EachCellType,  "_BySortedCaseID_20231030.pdf"); print(OutFile)
  ggsave(MyViolinplot, height=11,width=12, dpi=300, filename=OutFile, useDingbats=FALSE)  #  height=26,width=19,
}

```

## Step9. Calculate the ration CD163+ or MRC1+ macrophage over total macrophage in individual patients.

```{r}
ScaledExpDataCount <- SeuratObj_ERpos_SubsetMyeloid@assays$SCT$counts; dim(ScaledExpDataCount); ScaledExpDataCount[1:5,1:5] # 2439g1   974c
ScaledExpDataCount_CD163MRC1 <- ScaledExpDataCount %>% data.frame %>% dplyr::filter(grepl("CD163$|MRC1", rownames(ScaledExpDataCount))) %>% t %>% data.frame;
dim(ScaledExpDataCount_CD163MRC1); ScaledExpDataCount_CD163MRC1[1:2,] # 974 2

LoopNumb<-0
for(EachGene in colnames(ScaledExpDataCount_CD163MRC1)) {
    # EachGene <- colnames(ScaledExpDataCount_CD163MRC1)[2]
    LoopNumb <- LoopNumb+1; 
    ScaledExpDataCount_SubsetByGene <- ScaledExpDataCount_CD163MRC1 %>% dplyr::select_if(grepl(EachGene, colnames(ScaledExpDataCount_CD163MRC1)))
    dim(ScaledExpDataCount_SubsetByGene);head(ScaledExpDataCount_SubsetByGene) # 974 1

    ## subset by patient ID
    AllPt_PosCellRatio<-c()
    for (EachPt in Clinicaldata_ProportionCal_Sort$CaseID){
        # EachPt <- Clinicaldata_ProportionCal_Sort$CaseID[1]; print(EachPt)
        ScaledExpDataCount_SubsetByPt <- ScaledExpDataCount_SubsetByGene %>% dplyr::filter(grepl(EachPt, rownames(ScaledExpDataCount_SubsetByGene)))
        dim(ScaledExpDataCount_SubsetByPt) # 45 1
        
        colnames(ScaledExpDataCount_SubsetByPt)[1] <- "CurrentGene"
        # ExpData_NonZeroExp <- ScaledExpDataCount_SubsetByPt[ScaledExpDataCount_SubsetByPt$CurrentGene!=0, ] %>% data.frame
        
        CountExpPosCell <- table(ScaledExpDataCount_SubsetByPt$CurrentGene==0) %>% data.frame
        NumbCell_ExpPos <- CountExpPosCell$Freq[CountExpPosCell$Var1=="FALSE"]; print(paste0("Numb cell of exp Positive: ", NumbCell_ExpPos))
        
        TotalCellNumb <- nrow(ScaledExpDataCount_SubsetByPt); print(paste0("Total cell Numb: ", TotalCellNumb))
        
        PosCellRatio <- NumbCell_ExpPos/TotalCellNumb; print(paste0(EachPt, " has ratio: ", PosCellRatio)); 
        names(PosCellRatio) <- EachPt; print(PosCellRatio)
        
        AllPt_PosCellRatio <- c(AllPt_PosCellRatio, PosCellRatio)
    }
    
    AllPt_PosCellRatio_Df <- data.frame(AllPt_PosCellRatio); colnames(AllPt_PosCellRatio_Df)[1]<-paste0(EachGene,"_Ratio")
    
    if(LoopNumb==1) {
          AllPtAllGene_PosCellRatio <- AllPt_PosCellRatio_Df
    } else {
          AllPtAllGene_PosCellRatio <- cbind(AllPtAllGene_PosCellRatio, AllPt_PosCellRatio_Df)
    }
}

dim(AllPtAllGene_PosCellRatio); print(AllPtAllGene_PosCellRatio) # 10 2. 
fwrite(AllPtAllGene_PosCellRatio, file="PosCellRatio_MRC1CD163.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

## Step10. Make violin plot per cell types

```{r}
for(EachCellType in MyCellType_NoNASelc) {
  # EachCellType <- MyCellType_NoNASelc[1]
  if (EachCellType == "Macrophage") { 
    MyFeature <- Feature_Mac  
  } else if (EachCellType == "CD4Tcell") { 
    MyFeature <- Feature_CD4Tcell  
  } else if (EachCellType == "CD8Tcell") { 
    MyFeature <- Feature_CD8Tcell
  } else if (EachCellType == "CancerEpithelial") { 
    MyFeature <- Feature_CanEpi 
  } else if (EachCellType == "NKcell") { 
    MyFeature <- Feature_NKcell 

  }
  
  MyViolinplot <- VlnPlot(object = SeuratObj_ERpos_CellTypeMacrophage, features = c(MyFeature), group.by = "CellTypeMacroTcell_GSE176EGA6608")  # raster=TRUE requires ggrastr r packge.
  OutFile <- paste0("OutVlnplot_CytoskeletonExp_PerCellType.pdf"); print(OutFile)
  ggsave(MyViolinplot, height=15,width=22, dpi=300, filename=OutFile, useDingbats=FALSE)
}
```

## Step6. Make CellPhoneDB input files - Baseline count data and cell metadata for YoungElderly Rich vs Poor

```{r}
# ## 1) Young patients
# Remove zero count genes and save the count data for CellPhoneDB input
CountData_Young_CellType <- SeuratObj_ERpos_Young[["RNA"]]@counts; dim(CountData_Young_CellType)  # 29733  8119
CountData_Young_CellType_NoLowCountGene <- CountData_Young_CellType[rowSums(CountData_Young_CellType)>20, ]; dim(CountData_Young_CellType_NoLowCountGene); #  17249g  8119c
colnames(CountData_Young_CellType_NoLowCountGene) <- gsub("-","_",colnames(CountData_Young_CellType_NoLowCountGene))
write.table(CountData_Young_CellType_NoLowCountGene, "CountDataCellPhoneDB_GSE176078_Young_MainCellType_17249g8119c.txt", col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)  #

# Make CellPhoneDB Metadata for MacroPoor tumor cells
CellPhoneDB_MetaData_Young <- data.frame(SeuratObj_ERpos_Young@active.ident) %>% tibble::rownames_to_column("Cell");
colnames(CellPhoneDB_MetaData_Young)[2] <- "cell_type"; dim(CellPhoneDB_MetaData_Young) # 8119   2  
CellPhoneDB_MetaData_Young$Cell <- gsub("-","_",CellPhoneDB_MetaData_Young$Cell)   ## if there is "-" in the cell ID, cellphonedb gets error.
table(CellPhoneDB_MetaData_Young$Cell == colnames(CountData_Young_CellType_NoLowCountGene))  # Should be all TRUE
CellPhoneDB_MetaData_Young$cell_type <- gsub("-","_",CellPhoneDB_MetaData_Young$cell_type)
write.table(CellPhoneDB_MetaData_Young, "Metadata_CellPhoneDB_GSE176078_Young3s_MainCellType_8119c.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)  #

# ## Make metadata file for SCENIC
# SCENIC_MetaData_MacPoor <- SeuratObj_ERpos_CellTypeYoungElderly_Poor@meta.data[, c("nCount_SCT","nFeature_SCT","CellTypeMacroTcell_GSE176EGA6608")]
# SCENIC_MetaData_MacPoor_Reorder <- SCENIC_MetaData_MacPoor[, c(3,2,1)]
# colnames(SCENIC_MetaData_MacPoor_Reorder) <- c("CellType","nGene","nUMI")  # 26352  3
# fwrite(SCENIC_MetaData_MacPoor_Reorder, "Metadata_SCENIC_GSE176078EGA6608_MacPoor_MainCellType_26352c.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
# 
## 2) ElderlyRich  -- Can be used for SCENIC too.
# Remove zero count genes and save the count data for CellPhoneDB input.  Can be used for SCENIC too.
CountData_Elderly_CellType <- SeuratObj_ERpos_Elderly[["RNA"]]@counts; dim(CountData_Elderly_CellType)  # 29733g 10713 cells for 7 rich smp #  35929 13419c for 5 rich smp
CountData_Elderly_CellType_NoLowCountGene <- CountData_Elderly_CellType[rowSums(CountData_Elderly_CellType)>20, ]; dim(CountData_Elderly_CellType_NoLowCountGene); #  16423g 10713cells.  # 17106 13419 for 5 rich smp
colnames(CountData_Elderly_CellType_NoLowCountGene) <- gsub("-","_",colnames(CountData_Elderly_CellType_NoLowCountGene))
write.table(CountData_Elderly_CellType_NoLowCountGene, "CountDataCellPhoneDB_GSE176078_Elderly_MainCellType_16423g10713c.txt", col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)  #
#saveRDS(CountData_Elderly_CellType_NoLowCountGene, "CountData_GSE176078_Elderly_MainCellType_16423g10713c_dgCMatrix.rds" )

# Make CellPhoneDB Metadata for MacroRich tumor cells
CellPhoneDB_MetaData_Elderly <- data.frame(SeuratObj_ERpos_Elderly@active.ident) %>% tibble::rownames_to_column("Cell");
colnames(CellPhoneDB_MetaData_Elderly)[2] <- "cell_type"; dim(CellPhoneDB_MetaData_Elderly) # 10713   2  
CellPhoneDB_MetaData_Elderly$Cell <- gsub("-","_",CellPhoneDB_MetaData_Elderly$Cell)
table(CellPhoneDB_MetaData_Elderly$Cell == colnames(CountData_Elderly_CellType_NoLowCountGene)) ## Should be all TRUE
CellPhoneDB_MetaData_Elderly$cell_type <- gsub("-","_",CellPhoneDB_MetaData_Elderly$cell_type)
write.table(CellPhoneDB_MetaData_Elderly, "Metadata_CellPhoneDB_GSE176078_Elderly3s_MainCellType_10713c.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)  #

# ## Make metadata file for SCENIC
# SCENIC_MetaData_MacRich <- SeuratObj_ERpos_CellTypeYoungElderly_Rich@meta.data[, c("nCount_SCT","nFeature_SCT","CellTypeMacroTcell_GSE176EGA6608")]
# SCENIC_MetaData_MacRich_Reorder <- SCENIC_MetaData_MacRich[, c(3,2,1)]
# colnames(SCENIC_MetaData_MacRich_Reorder) <- c("CellType","nGene","nUMI")
# fwrite(SCENIC_MetaData_MacRich_Reorder, "Metadata_SCENIC_GSE176078EGA6608_MacRich_MainCellType_18897c.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

## 3) Whole subjects   -- Can be used for GSVA and SCENIC too. 
# Remove zero count genes and save the count data for CellPhoneDB input
CountData_YoungElderlyWhole_CellType <- as.data.frame(SeuratObj_ERpos_CellTypeYoungElderly_Whole[["RNA"]]@counts); dim(CountData_YoungElderlyWhole_CellType)  # 29733 30959 cells
CountData_YoungElderlyWhole_CellType_NoLowCountGene <- CountData_YoungElderlyWhole_CellType[rowSums(CountData_YoungElderlyWhole_CellType)>10, ]; dim(CountData_YoungElderlyWhole_CellType_NoLowCountGene); #   150: 15979 30959  # 10: 22057 30959
## remove genes of non-proteincoding
CountData_YoungElderlyWhole_CellType_NoNonProtCoding <- CountData_YoungElderlyWhole_CellType_NoLowCountGene %>% dplyr::filter(!grepl("\\.", rownames(CountData_YoungElderlyWhole_CellType_NoLowCountGene))); dim(CountData_YoungElderlyWhole_CellType_NoNonProtCoding) # 18063 30959
CountData_YoungElderlyWhole_CellType_NoLowCountGene <- CountData_YoungElderlyWhole_CellType_NoNonProtCoding[, colSums(CountData_YoungElderlyWhole_CellType_NoNonProtCoding)>1000 ]; dim(CountData_YoungElderlyWhole_CellType_NoLowCountGene); 
                                                                                               # 10/1000: 18063g 28732c
colnames(CountData_YoungElderlyWhole_CellType_NoLowCountGene) <- gsub("-","_",colnames(CountData_YoungElderlyWhole_CellType_NoLowCountGene))
fwrite(data.frame(CountData_YoungElderlyWhole_CellType_NoLowCountGene), "CountDataForGSVA_GSE176078_YoungElderlyWhole_MainCellType_18063g28732c.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)  # 200/2000 19221g86138c

CountData_YoungElderlyWhole_CellType_NoLowCountGene_Tp <- t(CountData_YoungElderlyWhole_CellType_NoLowCountGene)
saveRDS(CountData_YoungElderlyWhole_CellType_NoLowCountGene_Tp, file="CountData_YoungElderlyWhole_MainCellType_28732c18063g_Tp.rds")   ## I can't write a big file. So, I moved to Bridges2
fwrite(data.frame(CountData_YoungElderlyWhole_CellType_NoLowCountGene_Tp), "CountDataCellPhoneDB_GSE176078_YoungElderlyWhole_MainCellType_28732c18063g_Tp.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)  # 200/2000 18063g28732c

# Make CellPhoneDB Metadata for MacroWhole tumor cells
CellPhoneDB_MetaData_YoungElderlyWhole <- data.frame(SeuratObj_ERpos_CellTypeYoungElderly_Whole@active.ident) %>% tibble::rownames_to_column("Cell"); 
colnames(CellPhoneDB_MetaData_YoungElderlyWhole)[2] <- "cell_type"; dim(CellPhoneDB_MetaData_YoungElderlyWhole) # 30959   2
CellPhoneDB_MetaData_YoungElderlyWhole$Cell <- gsub("-","_",CellPhoneDB_MetaData_YoungElderlyWhole$Cell)
### Subset CellID based on countdata matrix. 
CellPhoneDB_MetaData_YoungElderlyWhole <- CellPhoneDB_MetaData_YoungElderlyWhole %>% dplyr::filter(CellPhoneDB_MetaData_YoungElderlyWhole$Cell %in% colnames(CountData_YoungElderlyWhole_CellType_NoLowCountGene)); dim(CellPhoneDB_MetaData_YoungElderlyWhole)  # 150/1000  105184
# 28732 2

## Make metadata file for MacroWhole SCENIC
SCENIC_MetaData_MacWhole <- SeuratObj_ERpos_CellTypeYoungElderly_Whole@meta.data[, c("nCount_SCT","nFeature_SCT","CellTypeAnnotSH")]
SCENIC_MetaData_MacWhole_Reorder <- SCENIC_MetaData_MacWhole[, c(3,2,1)]
colnames(SCENIC_MetaData_MacWhole_Reorder) <- c("CellType","nGene","nUMI"); dim(SCENIC_MetaData_MacWhole_Reorder) # 30959   3
fwrite(SCENIC_MetaData_MacWhole_Reorder, "Metadata_SCENIC_GSE176078_MacWhole_MainCellType_30959c.txt", col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

table(CellPhoneDB_MetaData_YoungElderlyWhole$Cell == colnames(CountData_YoungElderlyroWhole_CellType_NoLowCountGene))
CellPhoneDB_MetaData_YoungElderlyWhole$cell_type <- gsub("-","_",CellPhoneDB_MetaData_YoungElderlyWhole$cell_type)
write.table(CellPhoneDB_MetaData_YoungElderlyWhole, "Metadata_CellPhoneDB_GSE176078_YoungElderlyWhole_MainCellType_28732c.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)  # 
```

## To run Gene Set Variation Analysis (GSVA) and Pathway RespOnsive GENes for activity inference (PROGENy)

Please refer to our code at

-   https://github.com/leeoesterreich/BRCA_OlderWomen_BulkRNAseq

-   https://github.com/leeoesterreich/BRCA_Rat_YoungerOlder_snRNAseq

## To run WCSEA, please refer to our code at

-   https://github.com/wangxlab/indepthPathway
