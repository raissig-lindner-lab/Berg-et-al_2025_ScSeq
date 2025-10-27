### R script B. distachyon scSeq atlas


###### 1. Use SoupX to filter out ambient RNA ------------------------------------------------------------------
library(SoupX)

### Load data
wt_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea1_2/wt_2022_genome_v1_2/outs/")
sid_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea1_2/sid_2022_genome_v1_2/outs/")

wt170523_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea2023/wt_2023-17-05_genome_v1_2/outs/")
wt160523_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea2023/wt_2023-16-05_genome_v1_2/outs/")

wtines_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_ines/wt_2021_genome_v1_2/outs/")
sidines_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_ines/sid_2021_genome_v1_2/outs/")

wt2024_1_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea2024/wt_2024-1_genome_v1_2/outs/")
wt2024_2_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea2024/wt_2024-2_genome_v1_2/outs/")
wt2024_3_sc = load10X("E:/Scseq_brachy_dev/Scseq-data_Lea2024/wt_2024-3_genome_v1_2/outs/")




### Estimate rho
wt_sc = autoEstCont(wt_sc)
#833 genes passed tf-idf cut-off and 48 soup quantile filter.  Taking the top 48.
#Using 79 independent estimates of rho.
#Estimated global rho of 0.10

sid_sc = autoEstCont(sid_sc)
#947 genes passed tf-idf cut-off and 81 soup quantile filter.  Taking the top 81.
#Using 172 independent estimates of rho.
#Estimated global rho of 0.06

wt170523_sc = autoEstCont(wt170523_sc)
#533 genes passed tf-idf cut-off and 169 soup quantile filter.  Taking the top 100.
#Using 402 independent estimates of rho.
#Estimated global rho of 0.06

wt160523_sc = autoEstCont(wt160523_sc)
#275 genes passed tf-idf cut-off and 80 soup quantile filter.  Taking the top 80.
#Using 132 independent estimates of rho.
#Estimated global rho of 0.06

wtines_sc = autoEstCont(wtines_sc)
#576 genes passed tf-idf cut-off and 305 soup quantile filter.  Taking the top 100.
#Using 495 independent estimates of rho.
#Estimated global rho of 0.10

sidines_sc = autoEstCont(sidines_sc)
#263 genes passed tf-idf cut-off and 102 soup quantile filter.  Taking the top 100.
#Using 269 independent estimates of rho.
#Estimated global rho of 0.08

wt2024_1_sc = autoEstCont(wt2024_1_sc)
#300 genes passed tf-idf cut-off and 97 soup quantile filter.  Taking the top 97.
#Using 204 independent estimates of rho.
#Estimated global rho of 0.17

wt2024_2_sc = autoEstCont(wt2024_2_sc)
#347 genes passed tf-idf cut-off and 78 soup quantile filter.  Taking the top 78.
#Using 219 independent estimates of rho.
#Estimated global rho of 0.10

wt2024_3_sc = autoEstCont(wt2024_3_sc)
#353 genes passed tf-idf cut-off and 75 soup quantile filter.  Taking the top 75.
#Using 213 independent estimates of rho.
#Estimated global rho of 0.07




### Clean the data
wt_out = adjustCounts(wt_sc)
sid_out = adjustCounts(sid_sc)

wt170523_out = adjustCounts(wt170523_sc)
wt160523_out = adjustCounts(wt160523_sc)

wtines_out = adjustCounts(wtines_sc)
sidines_out = adjustCounts(sidines_sc)

wt2024_1_out = adjustCounts(wt2024_1_sc)
wt2024_2_out = adjustCounts(wt2024_2_sc)
wt2024_3_out = adjustCounts(wt2024_3_sc)








###### 2. Initiate Seurat objects -----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(metap)


### Initialise the Seurat object with the SoupX SoupChannel objects
wt <- CreateSeuratObject(wt_out)
sid <- CreateSeuratObject(sid_out)

wt170523 <- CreateSeuratObject(wt170523_out)
wt160523 <- CreateSeuratObject(wt160523_out)

wtines <- CreateSeuratObject(wtines_out)
sidines <- CreateSeuratObject(sidines_out)

wt2024_1 <- CreateSeuratObject(wt2024_1_out)
wt2024_2 <- CreateSeuratObject(wt2024_2_out)
wt2024_3 <- CreateSeuratObject(wt2024_3_out)





###### 3. Quality control to select cells for further analysis -----------------------------------------------------
### load files with mitochondrial and chloroplast genes
mitochondrial_genes <- read.table(file = "E:/Scseq_brachy_dev/scseq_general/brachy_mitochondrial_gene_annot.csv", header = FALSE, sep = ";")
mitochondrial_genes <- mitochondrial_genes$V1
chloroplast_genes <- read.table(file= "E:/Scseq_brachy_dev/scseq_general/brachy_chloroplast_gene_annot.csv", header = FALSE, sep = ";")
chloroplast_genes <- chloroplast_genes$V1
### in genome v.1.2, it cannot find all the genes of the chloroplast_genes list in the data and gives an error, so we need to take those out
chloroplast_genes_v2 <- chloroplast_genes[-c(25, 153, 194, 196, 200)]




### add metadata that gives percentages of mitochondrial and chloroplast DNA
wt[["percent.mt"]] <- PercentageFeatureSet(object = wt, features = mitochondrial_genes)
sid[["percent.mt"]] <- PercentageFeatureSet(object = sid, features = mitochondrial_genes)
wt[["percent.chl"]] <- PercentageFeatureSet(object = wt, features = chloroplast_genes_v2)
sid[["percent.chl"]] <- PercentageFeatureSet(object = sid, features = chloroplast_genes_v2)

wt170523[["percent.mt"]] <- PercentageFeatureSet(object = wt170523, features = mitochondrial_genes)
wt160523[["percent.mt"]] <- PercentageFeatureSet(object = wt160523, features = mitochondrial_genes)
wt170523[["percent.chl"]] <- PercentageFeatureSet(object = wt170523, features = chloroplast_genes_v2)
wt160523[["percent.chl"]] <- PercentageFeatureSet(object = wt160523, features = chloroplast_genes_v2)

wtines[["percent.mt"]] <- PercentageFeatureSet(object = wtines, features = mitochondrial_genes)
sidines[["percent.mt"]] <- PercentageFeatureSet(object = sidines, features = mitochondrial_genes)
wtines[["percent.chl"]] <- PercentageFeatureSet(object = wtines, features = chloroplast_genes_v2)
sidines[["percent.chl"]] <- PercentageFeatureSet(object = sidines, features = chloroplast_genes_v2)

wt2024_1[["percent.mt"]] <- PercentageFeatureSet(object = wt2024_1, features = mitochondrial_genes)
wt2024_2[["percent.mt"]] <- PercentageFeatureSet(object = wt2024_2, features = mitochondrial_genes)
wt2024_3[["percent.mt"]] <- PercentageFeatureSet(object = wt2024_3, features = mitochondrial_genes)
wt2024_1[["percent.chl"]] <- PercentageFeatureSet(object = wt2024_1, features = chloroplast_genes_v2)
wt2024_2[["percent.chl"]] <- PercentageFeatureSet(object = wt2024_2, features = chloroplast_genes_v2)
wt2024_3[["percent.chl"]] <- PercentageFeatureSet(object = wt2024_3, features = chloroplast_genes_v2)




### Filtering out cells
## We filter out cells with less than 500 and more than 10,000 genes (nFeature_RNA), 
## less than 1250 and more than 50,000 transcripts (nCount_RNA), 
## a % of mitochondrial counts (percent.mt) > 5% and of chloroplast counts (percent.chl) > 10%:
wt_str <- subset(wt, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)
sid_str <- subset(sid, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)

wt170523_str <- subset(wt170523, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)
wt160523_str <- subset(wt160523, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)

wt_strines <- subset(wtines, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)
sid_strines <- subset(sidines, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)

wt2024_1_str <- subset(wt2024_1, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)
wt2024_2_str <- subset(wt2024_2, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)
wt2024_3_str <- subset(wt2024_3, nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1250 & nCount_RNA < 50000 & percent.mt < 5 & percent.chl < 10)




###### 4. Normalizing the data (expression normalized within each cell) ------------------------------------------------
wt_str <- NormalizeData(wt_str, normalization.method = "LogNormalize", scale.factor = 10000)
sid_str <- NormalizeData(sid_str, normalization.method = "LogNormalize", scale.factor = 10000)

wt170523_str <- NormalizeData(wt170523_str, normalization.method = "LogNormalize", scale.factor = 10000)
wt160523_str <- NormalizeData(wt160523_str, normalization.method = "LogNormalize", scale.factor = 10000)

wt_strines <- NormalizeData(wt_strines, normalization.method = "LogNormalize", scale.factor = 10000)
sid_strines <- NormalizeData(sid_strines, normalization.method = "LogNormalize", scale.factor = 10000)

wt2024_1_str <- NormalizeData(wt2024_1_str, normalization.method = "LogNormalize", scale.factor = 10000)
wt2024_2_str <- NormalizeData(wt2024_2_str, normalization.method = "LogNormalize", scale.factor = 10000)
wt2024_3_str <- NormalizeData(wt2024_3_str, normalization.method = "LogNormalize", scale.factor = 10000)




###### 5. Identification of highly variable features (feature selection) -----------------------------------------------
wt_str <- FindVariableFeatures(wt_str, selection.method = "vst", nfeatures = 2000)
sid_str <- FindVariableFeatures(sid_str, selection.method = "vst", nfeatures = 2000)

wt170523_str <- FindVariableFeatures(wt170523_str, selection.method = "vst", nfeatures = 2000)
wt160523_str <- FindVariableFeatures(wt160523_str, selection.method = "vst", nfeatures = 2000)

wt_strines <- FindVariableFeatures(wt_strines, selection.method = "vst", nfeatures = 2000)
sid_strines <- FindVariableFeatures(sid_strines, selection.method = "vst", nfeatures = 2000)

wt2024_1_str <- FindVariableFeatures(wt2024_1_str, selection.method = "vst", nfeatures = 2000)
wt2024_2_str <- FindVariableFeatures(wt2024_2_str, selection.method = "vst", nfeatures = 2000)
wt2024_3_str <- FindVariableFeatures(wt2024_3_str, selection.method = "vst", nfeatures = 2000)




###### 6. Scale data (gene expression normalized across cells) --------------------------------------------------------
wt_str <- ScaleData(wt_str, features = rownames(wt_str), order = T)
sid_str <- ScaleData(sid_str, features = rownames(sid_str), order = T)

wt170523_str <- ScaleData(wt170523_str, features = rownames(wt170523_str), order = T)
wt160523_str <- ScaleData(wt160523_str, features = rownames(wt160523_str), order = T)

wt_strines <- ScaleData(wt_strines, features = rownames(wt_strines), order = T)
sid_strines <- ScaleData(sid_strines, features = rownames(sid_strines), order = T)

wt2024_1_str <- ScaleData(wt2024_1_str, features = rownames(wt2024_1_str), order = T)
wt2024_2_str <- ScaleData(wt2024_2_str, features = rownames(wt2024_2_str), order = T)
wt2024_3_str <- ScaleData(wt2024_3_str, features = rownames(wt2024_3_str), order = T)




###### 7. run PCA -----------------------------------------------------------------------------------------------
wt_str <- RunPCA(wt_str, features = VariableFeatures(object = wt_str), npcs = 50)
sid_str <- RunPCA(sid_str, features = VariableFeatures(object = sid_str), npcs = 50)

wt170523_str <- RunPCA(wt170523_str, features = VariableFeatures(object = wt170523_str), npcs = 50)
wt160523_str <- RunPCA(wt160523_str, features = VariableFeatures(object = wt160523_str), npcs = 50)

wt_strines <- RunPCA(wt_strines, features = VariableFeatures(object = wt_strines), npcs = 50)
sid_strines <- RunPCA(sid_strines, features = VariableFeatures(object = sid_strines), npcs = 50)

wt2024_1_str <- RunPCA(wt2024_1_str, features = VariableFeatures(object = wt2024_1_str), npcs = 50)
wt2024_2_str <- RunPCA(wt2024_2_str, features = VariableFeatures(object = wt2024_2_str), npcs = 50)
wt2024_3_str <- RunPCA(wt2024_3_str, features = VariableFeatures(object = wt2024_3_str), npcs = 50)




###### 8. run UMAP ------------------------------------------------------------------------------------------------------------------
### check elbow plots first
ElbowPlot(wt_str, ndims = 50)  # --> 30 PCs
ElbowPlot(sid_str, ndims = 50) # --> 30 PCs

ElbowPlot(wt170523_str, ndims = 50)  # --> 30 PCs
ElbowPlot(wt160523_str, ndims = 50) # --> 30 PCs

ElbowPlot(wt_strines, ndims = 50) # --> 30 PCs
ElbowPlot(sid_strines, ndims = 50) # --> 30 PCs

ElbowPlot(wt2024_1_str, ndims = 50)  # --> 30 PCs
ElbowPlot(wt2024_2_str, ndims = 50)  # --> 30 PCs
ElbowPlot(wt2024_3_str, ndims = 50)  # --> 30 PCs



### run UMAP with the chosen number of dimensions
wt_str <- RunUMAP(wt_str, reduction = "pca", dims = 1:30)
sid_str <- RunUMAP(sid_str, reduction = "pca", dims = 1:30)

wt170523_str <- RunUMAP(wt170523_str, reduction = "pca", dims = 1:30)
wt160523_str <- RunUMAP(wt160523_str, reduction = "pca", dims = 1:30)

wt_strines <- RunUMAP(wt_strines, reduction = "pca", dims = 1:30)
sid_strines <- RunUMAP(sid_strines, reduction = "pca", dims = 1:30)

wt2024_1_str <- RunUMAP(wt2024_1_str, reduction = "pca", dims = 1:30)
wt2024_2_str <- RunUMAP(wt2024_2_str, reduction = "pca", dims = 1:30)
wt2024_3_str <- RunUMAP(wt2024_3_str, reduction = "pca", dims = 1:30)




###### 9. Check for doublets and exclude them from the data ---------------------------------------------------------
### 9.1. Find the right pK values for my datasets (no ground-truth) 

## 2022 datasets
sweep.res.list_wt_str <- paramSweep(wt_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt_str <- summarizeSweep(sweep.res.list_wt_str, GT = FALSE)
bcmvn_wt_str <- find.pK(sweep.stats_wt_str)
pK_value <- bcmvn_wt_str %>% filter(bcmvn_wt_str$BCmetric == max(bcmvn_wt_str$BCmetric))
pK_value # --> pK 0.26

sweep.res.list_sid_str <- paramSweep(sid_str, PCs = 1:30, sct = FALSE)
sweep.stats_sid_str <- summarizeSweep(sweep.res.list_sid_str, GT = FALSE)
bcmvn_sid_str <- find.pK(sweep.stats_sid_str)
pK_value <- bcmvn_sid_str %>% filter(bcmvn_sid_str$BCmetric == max(bcmvn_sid_str$BCmetric))
pK_value # --> pK 0.07


## 2023 datasets
sweep.res.list_wt170523_str <- paramSweep(wt170523_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt170523_str <- summarizeSweep(sweep.res.list_wt170523_str, GT = FALSE)
bcmvn_wt170523_str <- find.pK(sweep.stats_wt170523_str)
pK_value <- bcmvn_wt170523_str %>% filter(bcmvn_wt170523_str$BCmetric == max(bcmvn_wt170523_str$BCmetric))
pK_value # --> pK 0.19

sweep.res.list_wt160523_str <- paramSweep(wt160523_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt160523_str <- summarizeSweep(sweep.res.list_wt160523_str, GT = FALSE)
bcmvn_wt160523_str <- find.pK(sweep.stats_wt160523_str)
pK_value <- bcmvn_wt160523_str %>% filter(bcmvn_wt160523_str$BCmetric == max(bcmvn_wt160523_str$BCmetric))
pK_value # --> pK 0.14


## 2021 datasets
sweep.res.list_wt_strines <- paramSweep(wt_strines, PCs = 1:30, sct = FALSE)
sweep.stats_wt_strines <- summarizeSweep(sweep.res.list_wt_strines, GT = FALSE)
bcmvn_wt_strines <- find.pK(sweep.stats_wt_strines)
pK_value <- bcmvn_wt_strines %>% filter(bcmvn_wt_strines$BCmetric == max(bcmvn_wt_strines$BCmetric))
pK_value # --> pK 0.01

sweep.res.list_sid_strines <- paramSweep(sid_strines, PCs = 1:30, sct = FALSE)
sweep.stats_sid_strines <- summarizeSweep(sweep.res.list_sid_strines, GT = FALSE)
bcmvn_sid_strines <- find.pK(sweep.stats_sid_strines)
pK_value <- bcmvn_sid_strines %>% filter(bcmvn_sid_strines$BCmetric == max(bcmvn_sid_strines$BCmetric))
pK_value # --> pK 0.24


## 2024 datasets
sweep.res.list_wt2024_1_str <- paramSweep(wt2024_1_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt2024_1_str <- summarizeSweep(sweep.res.list_wt2024_1_str, GT = FALSE)
bcmvn_wt2024_1_str <- find.pK(sweep.stats_wt2024_1_str)
pK_value <- bcmvn_wt2024_1_str %>% filter(bcmvn_wt2024_1_str$BCmetric == max(bcmvn_wt2024_1_str$BCmetric))
pK_value # --> pK 0.3

sweep.res.list_wt2024_2_str <- paramSweep(wt2024_2_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt2024_2_str <- summarizeSweep(sweep.res.list_wt2024_2_str, GT = FALSE)
bcmvn_wt2024_2_str <- find.pK(sweep.stats_wt2024_2_str)
pK_value <- bcmvn_wt2024_2_str %>% filter(bcmvn_wt2024_2_str$BCmetric == max(bcmvn_wt2024_2_str$BCmetric))
pK_value # --> pK 0.21

sweep.res.list_wt2024_3_str <- paramSweep(wt2024_3_str, PCs = 1:30, sct = FALSE)
sweep.stats_wt2024_3_str <- summarizeSweep(sweep.res.list_wt2024_3_str, GT = FALSE)
bcmvn_wt2024_3_str <- find.pK(sweep.stats_wt2024_3_str)
pK_value <- bcmvn_wt2024_3_str %>% filter(bcmvn_wt2024_3_str$BCmetric == max(bcmvn_wt2024_3_str$BCmetric))
pK_value # --> pK 0.23




### 9.2. Homotypic Doublet Proportion Estimate 
## we need the expected doublet formation rate for the data set -- I chose those based approximately on the 10X V3 user guide:
## https://assets.ctfassets.net/an68im79xiti/7HMWNGXDSdJXbXNLtbEvCZ/5c4989e1990fc3abcf4dc811161b9f21/CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevE.pdf
## wt_str: 8% multiplet rate for about 16500 cells loaded
## sid_str: 8% multiplet rate for about 16500 cells loaded
## wt_strines: 6% (between 5.6% and 6.4% for about 12000 cells loaded so I decided on 6%)
## sid_strines: 6% multiplet rate for about 12000 cells loaded
## wt170523_str: 10% multiplet rate
## wt160523_str: 10% multiplet rate
## 2024 datasets: 10% multiplet rate

nExp_poi_wt_str <- round(0.08*nrow(wt_str@meta.data))
nExp_poi_sid_str <- round(0.08*nrow(sid_str@meta.data))

nExp_poi_wt170523_str <- round(0.1*nrow(wt170523_str@meta.data))
nExp_poi_wt160523_str <- round(0.1*nrow(wt160523_str@meta.data))

nExp_poi_wt_strines <- round(0.06*nrow(wt_strines@meta.data))
nExp_poi_sid_strines <- round(0.06*nrow(sid_strines@meta.data))

nExp_poi_wt2024_1_str <- round(0.1*nrow(wt2024_1_str@meta.data))
nExp_poi_wt2024_2_str <- round(0.1*nrow(wt2024_2_str@meta.data))
nExp_poi_wt2024_3_str <- round(0.1*nrow(wt2024_3_str@meta.data))




### 9.3. Run DoubletFinder
wt_str <- doubletFinder(wt_str, PCs = 1:30, pK = 0.26, nExp = nExp_poi_wt_str) #Soupx filtered
sid_str <- doubletFinder(sid_str, PCs = 1:30, pK = 0.07, nExp = nExp_poi_sid_str) #Soupx filtered

wt170523_str <- doubletFinder(wt170523_str, PCs = 1:30, pK = 0.19, nExp = nExp_poi_wt170523_str) #Soupx filtered
wt160523_str <- doubletFinder(wt160523_str, PCs = 1:30, pK = 0.14, nExp = nExp_poi_wt160523_str) #Soupx filtered

wt_strines <- doubletFinder(wt_strines, PCs = 1:30, pK = 0.01, nExp = nExp_poi_wt_strines) #Soupx filtered
sid_strines <- doubletFinder(sid_strines, PCs = 1:30, pK = 0.24, nExp = nExp_poi_sid_strines) #Soupx filtered

wt2024_1_str <- doubletFinder(wt2024_1_str, PCs = 1:30, pK = 0.3, nExp = nExp_poi_wt2024_1_str)
wt2024_2_str <- doubletFinder(wt2024_2_str, PCs = 1:30, pK = 0.21, nExp = nExp_poi_wt2024_2_str)
wt2024_3_str <- doubletFinder(wt2024_3_str, PCs = 1:30, pK = 0.23, nExp = nExp_poi_wt2024_3_str)




### 9.4 Visualize the doublets
## you can find the right column to plot by typing wt_str$ and one of the options starts with "DF.classification", that should be your doublet or singlet classifier)
DimPlot(wt_str, reduction = "umap", group.by = "DF.classifications_0.25_0.26_600")
DimPlot(sid_str, reduction = "umap", group.by = "DF.classifications_0.25_0.07_647")

DimPlot(wt170523_str, reduction = "umap", group.by = "DF.classifications_0.25_0.19_297")
DimPlot(wt160523_str, reduction = "umap", group.by = "DF.classifications_0.25_0.14_555")

DimPlot(wt_strines, reduction = "umap", group.by = "DF.classifications_0.25_0.01_170")
DimPlot(sid_strines, reduction = "umap", group.by = "DF.classifications_0.25_0.24_219")

DimPlot(wt2024_1_str, reduction = "umap", group.by = "DF.classifications_0.25_0.3_1493")
DimPlot(wt2024_2_str, reduction = "umap", group.by = "DF.classifications_0.25_0.21_1283")
DimPlot(wt2024_3_str, reduction = "umap", group.by = "DF.classifications_0.25_0.23_1845")




###### 10. Exclude doublets from the datasets and save the new objects --------------------------------------------
wt_str_nod <- subset(wt_str, subset = DF.classifications_0.25_0.26_600 == "Singlet")
sid_str_nod <- subset(sid_str, subset = DF.classifications_0.25_0.07_647 == "Singlet")

wt170523_str_nod <- subset(wt170523_str, subset = DF.classifications_0.25_0.19_297 == "Singlet")
wt160523_str_nod <- subset(wt160523_str, subset = DF.classifications_0.25_0.14_555 == "Singlet")

wt_strines_nod <- subset(wt_strines, subset = DF.classifications_0.25_0.01_170 == "Singlet")
sid_strines_nod <- subset(sid_strines, subset = DF.classifications_0.25_0.24_219 == "Singlet")

wt2024_1_str_nod <- subset(wt2024_1_str, subset = DF.classifications_0.25_0.3_1493 == "Singlet")
wt2024_2_str_nod <- subset(wt2024_2_str, subset = DF.classifications_0.25_0.21_1283 == "Singlet")
wt2024_3_str_nod <- subset(wt2024_3_str, subset = DF.classifications_0.25_0.23_1845 == "Singlet")

dim(wt_str_nod)
# 39068 genes across 6904 cells
dim(sid_str_nod)
# 39068 genes across 7437 cells

dim(wt170523_str_nod)
# 39068 genes across 2674 cells
dim(wt160523_str_nod)
# 39068 genes across 4995 cells

dim(wt_strines_nod)
# 39068 genes across 2656 cells
dim(sid_strines_nod)
# 39068 genes across 3427 cells

dim(wt2024_1_str_nod)
# 39068 genes across 13441 cells
dim(wt2024_2_str_nod)
# 39068 genes across 11545 cells
dim(wt2024_3_str_nod)
# 39068 genes across 16607 cells


### add a column "origin" specifying what dataset it is
wt_str_nod$origin <- "WT_2022"
sid_str_nod$origin <- "sid_2022"
wt170523_str_nod$origin <- "WT_2023-05-17"
wt160523_str_nod$origin <- "WT_2023-05-16"
wt_strines_nod$origin <- "WT_2021"
sid_strines_nod$origin <- "sid_2021"
wt2024_1_str_nod$origin <- "WT_2024_1"
wt2024_2_str_nod$origin <- "WT_2024_2"
wt2024_3_str_nod$origin <- "WT_2024_3"


### Save as .rds file
saveRDS(wt_str_nod, "wt_str_nod_genome_v1_2.rds")
saveRDS(sid_str_nod, "sid_str_nod_genome_v1_2.rds")
saveRDS(wt170523_str_nod, "wt170523_str_nod_genome_v1_2.rds")
saveRDS(wt160523_str_nod, "wt160523_str_nod_genome_v1_2.rds")
saveRDS(wt_strines_nod, "wt_strines_nod_genome_v1_2.rds")
saveRDS(sid_strines_nod, "sid_strines_nod_genome_v1_2.rds")
saveRDS(wt2024_1_str_nod, "wt2024_1_str_nod_genome_v1_2.rds")
saveRDS(wt2024_2_str_nod, "wt2024_2_str_nod_genome_v1_2.rds")
saveRDS(wt2024_3_str_nod, "wt2024_3_str_nod_genome_v1_2.rds")


### Reload data
wt_str_nod <- readRDS("wt_str_nod_genome_v1_2.rds")
sid_str_nod <- readRDS("sid_str_nod_genome_v1_2.rds")
wt170523_str_nod <- readRDS("wt170523_str_nod_genome_v1_2.rds")
wt160523_str_nod <- readRDS("wt160523_str_nod_genome_v1_2.rds")
wt_strines_nod <- readRDS("wt_strines_nod_genome_v1_2.rds")
sid_strines_nod <- readRDS("sid_strines_nod_genome_v1_2.rds")
wt2024_1_str_nod <- readRDS("wt2024_1_str_nod_genome_v1_2.rds")
wt2024_2_str_nod <- readRDS("wt2024_2_str_nod_genome_v1_2.rds")
wt2024_3_str_nod <- readRDS("wt2024_3_str_nod_genome_v1_2.rds")




###### 11. Merge datasets and run find variable features, scaling, PCA and UMAP on the merged data --------------------
wt_sid_nocca <- merge(wt_str_nod, y=list(sid_str_nod, wt_strines_nod, sid_strines_nod, wt170523_str_nod, wt160523_str_nod,
                                         wt2024_1_str_nod, wt2024_2_str_nod, wt2024_3_str_nod))

wt_sid_nocca <- FindVariableFeatures(wt_sid_nocca, selection.method = "vst", nfeatures = 2000)

wt_sid_nocca <- ScaleData(wt_sid_nocca)

wt_sid_nocca <- RunPCA(wt_sid_nocca, features = VariableFeatures(object = wt_sid_nocca), npcs = 100)

ElbowPlot(wt_sid_nocca, ndims = 50)

wt_sid_nocca <- FindNeighbors(wt_sid_nocca, reduction = "pca", dims = 1:40)

### If you change the resolution, more or less clusters will be found but the underlying plot will look the same
wt_sid_nocca <- FindClusters(wt_sid_nocca, resolution = 0.8)

wt_sid_nocca <- RunUMAP(wt_sid_nocca, reduction = "pca", dims = 1:40)


###### 12. Visualize the dataset --------------------------------------------------------------------------------------
DimPlot(wt_sid_nocca, reduction = "umap", group.by = "origin", shuffle = T)
DimPlot(wt_sid_nocca, reduction = "umap", group.by = "seurat_clusters", label = T)

### add meta data that distinguishes between WT and sid (bdmute)
## for soupx-filtered data
wt_sid_nocca$genotype <- ifelse(test = wt_sid_nocca$origin %in% "sid_2021", yes = "bdmute", 
                                no = ifelse(test = wt_sid_nocca$origin %in% "sid_2022", yes = "bdmute", no = "WT"))

DimPlot(wt_sid_nocca, reduction = "umap", group.by = "genotype", shuffle = T)


### add meta data that distinguishes between different tissue datasets (2021+2022+2024 (leaf dev. zone) vs 2023 (vSAM + leaf primordia))
## for soupx-filtered data
wt_sid_nocca$tissue <- ifelse(test = wt_sid_nocca$origin %in% "WT_2023-05-16", yes = "SAM + leaf primordia", 
                              no = ifelse(test = wt_sid_nocca$origin %in% "WT_2023-05-17", yes = "SAM + leaf primordia", 
                                          no = "leaf developmental zone"))

DimPlot(wt_sid_nocca, reduction = "umap", group.by = "tissue", shuffle = T)


### Add metadata on cell cycle phases
## save g2/m cell cycle genes
g2m.genes <- c("BdiBd21-3.4G0363900", "BdiBd21-3.3G0530200","BdiBd21-3.2G0144200",
               "BdiBd21-3.2G0281100","BdiBd21-3.2G0675300","BdiBd21-3.1G0398300", 
               "BdiBd21-3.5G0241300")
s.genes <- c("BdiBd21-3.1G1030800", "BdiBd21-3.3G0027100", "BdiBd21-3.3G0090800", 
             "BdiBd21-3.3G0550700", "BdiBd21-3.4G0249100")

## Join layers to enable cell cycle scoring
wt_sid_nocca <- JoinLayers(wt_sid_nocca)

### give each cell a score based on cell cycle genes
wt_sid_nocca <- CellCycleScoring(wt_sid_nocca, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(wt_sid_nocca, reduction = "umap", group.by = "Phase")


### save dataset
saveRDS(wt_sid_nocca, "wt_sid_nocca_all_genome_v1_2.rds")


### reload the data 
wt_sid_nocca <- readRDS("wt_sid_nocca.rds")




###### 13. Find out some numerical parameters -----------------------------------------------------------------------
### Number of cells
dim(wt_sid_nocca)
# 69687 cells


### Mean number of UMIs
mean(wt_sid_nocca$nCount_RNA) # 8205

### Median number of UMIs
median(wt_sid_nocca$nCount_RNA) # 5673


### Mean number of features/genes
mean(wt_sid_nocca$nFeature_RNA) # 2944

### Median number of features/genes
median(wt_sid_nocca$nFeature_RNA) # 2660


### values for each of the included datasets
stats_scseq %>% group_by(ori) %>% summarise(mean_umi=mean(umi), median_umi=median(umi), 
                                     mean_fea=mean(fea), median_fea=median(fea), 
                                     cells=length(fea))



### How many genes are expressed in total?
## calculate number of cells expressing each individual gene
all_genes <- Matrix::rowSums(wt_sid_nocca@assays$RNA$counts > 0)
all_genes <- as.data.frame(all_genes)
## filter out genes that are not expressed by any cell
expressed_genes <- all_genes %>% filter(all_genes>0) # 35584 out of 39068 genes are expressed in at least one cell
expressed_genes_strict <- all_genes %>% filter(all_genes>=5) # 31988 out of 39068 genes are expressed in at least 5 cells




###### 14. Visualize quality control plots -------------------------------------------------------------------------
FeaturePlot(wt_sid_nocca, features = "nFeature_RNA") 
FeaturePlot(wt_sid_nocca, features = "nCount_RNA") 
FeaturePlot(wt_sid_nocca, features = "percent.mt") 
FeaturePlot(wt_sid_nocca, features = "percent.chl")



###### 15. Tabulate cell statistics -----------------------------------------------------------------------------------
### How many cells are in each cluster?
table(wt_sid_nocca$seurat_clusters)

### How many cells are in each replicate?
table(wt_sid_nocca$origin)

### How does cluster membership vary by genotype?
table(Idents(wt_sid_nocca), wt_sid_nocca$genotype)




###### 16. Check marker gene expression on the dataset -------------------------------------------------------------
### Using a dot plot 
## if it says "-like" it means there were either several orthologues or it was not a 1:1 orthology
DotPlot(wt_sid_nocca, features = c("BdiBd21-3.5G0316500" # WOX4
                                   , "BdiBd21-3.5G0168500" # FCP1
                                   , "BdiBd21-3.1G0402200" # FON1/CLV1
                                   , "BdiBd21-3.2G0427700" # FEA2
                                   , "BdiBd21-3.2G0010000" # FEA3
                                   , "BdiBd21-3.1G0571300" # FEA4
                                   , "BdiBd21-3.1G0588300" # PIN1a
                                   , "BdiBd21-3.1G0135700" # KN1-like
                                   , "BdiBd21-3.1G0773000" # KNAT1-like1
                                   , "BdiBd21-3.1G0348000" # TMO6-like
                                   , "BdiBd21-3.1G1023300" # VND-like
                                   , "BdiBd21-3.5G0221500" # VND-like
                                   , "BdiBd21-3.2G0501500" # XCP1-like
                                   , "BdiBd21-3.3G0071200" # APL
                                   , "BdiBd21-3.2G0203400" # SLAH2
                                   , "BdiBd21-3.3G0204900" # PDF2-like 
                                   , "BdiBd21-3.2G0082000" # PDF1
                                   , "BdiBd21-3.2G0587500" # WOX9C-like1
                                   , "BdiBd21-3.4G0439300" # LHCA6
                                   , "BdiBd21-3.2G0749400" # STOMAGEN-1
),
cols = c("grey77", "black"), 
col.min = 0) +
  scale_x_discrete(labels=c("BdiBd21-3.5G0316500" = "BdWOX4", 
                            "BdiBd21-3.5G0168500" = "BdFCP1",
                            "BdiBd21-3.1G0402200" = "BdCLV1", 
                            "BdiBd21-3.2G0427700" = "BdFEA2",
                            "BdiBd21-3.2G0010000" = "BdFEA3", 
                            "BdiBd21-3.1G0571300" = "BdFEA4",
                            "BdiBd21-3.1G0588300" = "BdPIN1a",
                            "BdiBd21-3.1G0135700" = "BdKN1", 
                            "BdiBd21-3.1G0773000" = "BdKNAT1-like1",
                            "BdiBd21-3.1G0348000" = "BdTMO6-like",
                            "BdiBd21-3.1G1023300" = "BdVND-like1",
                            "BdiBd21-3.5G0221500" = "BdVND-like2",
                            "BdiBd21-3.2G0501500" = "BdXCP1-like",
                            "BdiBd21-3.3G0071200" = "BdAPL",
                            "BdiBd21-3.2G0203400" = "BdSLAH2",
                            "BdiBd21-3.3G0204900" = "BdPDF2-like",
                            "BdiBd21-3.2G0082000" = "BdPDF1",
                            "BdiBd21-3.2G0587500" = "BdWOX9C-like1",
                            "BdiBd21-3.2G0749400" = "BdSTOMAGEN-1",
                            "BdiBd21-3.4G0439300" = "BdLHCA6")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 8)) +
  labs(y = NULL, x = NULL)



###### 17. Rename clusters with new annotations ----------------------------------------------------------------------
Idents(wt_sid_nocca) <- "seurat_clusters"

### assign identities according to marker gene expression profiles
wt_sid_nocca <- RenameIdents(wt_sid_nocca, `0` = "Mesophyll", `1` = "Mesophyll", 
                             `2` = "Shoot apex", `3` = "Vasculature", 
                             `4` = "Epidermis", `5` = "Mesophyll", 
                             `6` = "Mesophyll", `7` = "Epidermis", 
                             `8` = "Mesophyll", `9` = "Mesophyll", 
                             `10` = "Epi.-Meso.-Vasc.", `11` = "Epidermis", 
                             `12` = "Meso.-Vasc.", `13` = "Epidermis", 
                             `14` = "Mesophyll", `15` = "Epi.-Meso.", 
                             `16` = "Vasculature", `17` = "Epidermis",
                             `18` = "Mesophyll", `19` = "Vasculature",
                             `20` = "Vasculature", `21` = "Vasculature",
                             `22` = "Vasculature", `23` = "Epidermis",
                             `24` = "Mesophyll", `25` = "Epidermis",
                             `26` = "Vasculature", `27` = "Epidermis",
                             `28` = "Mesophyll")

wt_sid_nocca$intermediate <- Idents(wt_sid_nocca)


### to merge the clusters with overlapping identities between tissues
wt_sid_nocca <- RenameIdents(wt_sid_nocca, "Epi.-Meso.-Vasc." = "Unknown",
                             "Epi.-Meso." = "Unknown",
                             "Meso.-Vasc." = "Unknown")

wt_sid_nocca$labels <- Idents(wt_sid_nocca) # save as metadata column

DimPlot(wt_sid_nocca, reduction = "umap", group.by = "labels")




###### 18. Identify markers ----------------------------------------------------------------------------------------
Idents(wt_sid_nocca) <- "labels"

all_markers <- FindAllMarkers(wt_sid_nocca)
all_markers <- filter(all_markers, p_val_adj < 0.5)
names(all_markers)[6] <- "tissue"


## to look for markers in a specific cluster
cluster.markers <- FindMarkers(data, ident.1 = "Epidermis", min.pct = 0.25) # markers for epidermal cluster that are expressed in at least 25% of the cells in that cluster
cluster.markers <- filter(cluster.markers, p_val_adj < 0.5) # filter for significant markers


write.table(all_markers, "E:/Scseq_brachy_dev/markers/all/all_markers.txt", row.names = F)








#---------------------------------------------------------------------------------------------------------------------

###### 19. Subset epidermal dataset ---------------------------------------------------------------------------------
Idents(wt_sid_nocca) <- "intermediate"
epidermis <- subset(wt_sid_nocca, idents = c("Epidermis", "Shoot apex", "Epi.-Meso.-Vasc.", "Epi.-Meso."))

dim(epidermis) # 24926 cells before iterative filtering


### 19.1 Find variable features
epidermis <- FindVariableFeatures(epidermis, selection.method = "vst", nfeatures = 2000)

### 19.2 Scale data
epidermis <- ScaleData(epidermis, features = rownames(epidermis))

### 19.3 run PCA
epidermis <- RunPCA(epidermis, features = VariableFeatures(object = epidermis), npcs = 30)

ElbowPlot(epidermis, ndims = 30)

### 19.4 find neighbours and clusters 
epidermis <- FindNeighbors(epidermis, dims = 1:25)
epidermis <- FindClusters(epidermis, resolution = 0.8)

### 19.5. run UMAP and visualise 
epidermis <- RunUMAP(epidermis, dims = 1:25)

DimPlot(epidermis)


### 19.6. take out clusters with mixed tissue identity (e.g. also mesophyllar or vascular markers) 
epidermis <- RenameIdents(epidermis, `0` = "Mixed tissues", `1` = "Stem cells", 
                          `2` = "Silica cells", `3` = "Unknown", 
                          `4` = "Prickle hair cells", `5` = "Stem cells", 
                          `6` = "Unknown", `7` = "Mixed tissues", 
                          `8` = "Mixed tissues", `9` = "Stem cells", 
                          `10` = "Stem cells", `11` = "Unknown", 
                          `12` = "Prickle hair cells", `13` = "Stomatal lineage",
                          `14` = "Stem cells", `15` = "Prickle hair cells",
                          `16` = "Stomatal lineage", `17` = "Stomatal lineage")
DimPlot(epidermis)

epidermis <- subset(epidermis, idents = c("Stem cells", "Prickle hair cells", "Stomatal lineage", "Unknown", "Silica cells"))
dim(epidermis) # now 18797 cells

epidermis <- FindVariableFeatures(epidermis, selection.method = "vst", nfeatures = 2000)
epidermis <- ScaleData(epidermis, features = rownames(epidermis))
epidermis <- RunPCA(epidermis, features = VariableFeatures(object = epidermis), npcs = 30)

epidermis <- FindNeighbors(epidermis, dims = 1:25)
epidermis <- FindClusters(epidermis, resolution = 0.8)
epidermis <- RunUMAP(epidermis, dims = 1:25)
DimPlot(epidermis)


### 19.7. additionally take out some of the shoot apical clusters
epidermis <- subset(epidermis, idents = c(0:2, 4:8, 10:17, 19))
dim(epidermis) # now 16108 cells

epidermis <- FindVariableFeatures(epidermis, selection.method = "vst", nfeatures = 2000)
epidermis <- ScaleData(epidermis, features = rownames(epidermis))
epidermis <- RunPCA(epidermis, features = VariableFeatures(object = epidermis), npcs = 30)

epidermis <- FindNeighbors(epidermis, dims = 1:25)
epidermis <- FindClusters(epidermis, resolution = 0.8)
epidermis <- RunUMAP(epidermis, dims = 1:25)

DimPlot(epidermis) # this is the final epidermal subset used in the publication


### 19.8. visualize
epidermis <- RenameIdents(epidermis, `0` = "Silica cell lineage", `1` = "Silica cell lineage", 
                          `2` = "Stomatal lineage", `3` = "Stem cells", 
                          `4` = "Unknown", `5` = "Stem cells", 
                          `6` = "Hair cell lineage", `7` = "Hair cell lineage", 
                          `8` = "Stomatal lineage", `9` = "Stomatal lineage", 
                          `10` = "Hair cell lineage", `11` = "Stage 0-1", 
                          `12` = "Hair cell lineage", `13` = "Hair cell lineage",
                          `14` = "Stomatal lineage", `15` = "Hair cell lineage",
                          `16` = "Stomatal lineage", `17` = "Stomatal lineage",
                          `18` = "Silica cell lineage", `19` = "Hair cell lineage")

epidermis$celltypes <- Idents(epidermis)

DimPlot(epidermis, label = F, shuffle = T, group.by = "celltypes") +
  scale_colour_manual(values = c("#574571", "#7fa074", "#c1d1aa", "lightgrey", "#2c4b27", "#b695bc")) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3))


DimPlot(epidermis, reduction = "umap", group.by = "seurat_clusters", label = T) + 
  scale_colour_manual(values = met.brewer("Hiroshige", n=20))


DimPlot(epidermis, reduction = "umap", group.by = "Phase") + theme(legend.position = "bottom") +
  scale_colour_manual(values = met.brewer("Derain", n=4))


### 19.9. mark epidermal cells on whole dataset
epi_cells <- WhichCells(epidermis)
DimPlot(wt_sid_nocca, cells.highlight = c(epi_cells), sizes.highlight = 0.3)


### 19.10. mark HC lineage on epidermal dataset
HCs <- WhichCells(epidermis, idents = c("Stage 0-1", "Hair cell lineage"))
DimPlot(epidermis, cells.highlight = c(HCs), sizes.highlight = 0.3)











#-------------------------------------------------------------------------------------------------

###### 20. Subset stomatal lineage subset ---------------------------------------------------------
Idents(epidermis) <- "seurat_clusters"

stomatal_files <- subset(epidermis, idents = c(2, 8, 9, 11, 14, 16, 17))
dim(stomatal_files) # 4675 cells

### 20.1. Find variable features 
stomatal_files <- FindVariableFeatures(stomatal_files, selection.method = "vst", nfeatures = 2000)

### 20.2. Scale data
stomatal_files <- ScaleData(stomatal_files, features = rownames(stomatal_files))

### 20.3. run PCA
stomatal_files <- RunPCA(stomatal_files, features = VariableFeatures(object = stomatal_files), npcs = 30)

ElbowPlot(stomatal_files, ndims = 30)

### 20.4. find neighbours and clusters
stomatal_files <- FindNeighbors(stomatal_files, dims = 1:25)
stomatal_files <- FindClusters(stomatal_files, resolution = 0.8)

### 20.5. run UMAP and visualise
stomatal_files <- RunUMAP(stomatal_files, dims = 1:25)

DimPlot(stomatal_files, reduction = "umap", group.by = "seurat_clusters", label = T) + 
  scale_colour_manual(values = met.brewer("Hiroshige", n=14))


## label datasets according to whether that year had a sid dataset (2021, 2022) or not (2023, 2024)
stomatal_files$year <- ifelse(test = stomatal_files$origin %in% "sid_2021", yes = "2021", 
                              no = ifelse(test = stomatal_files$origin %in% "WT_2021", yes = "2021", 
                                          no = ifelse(test = stomatal_files$origin %in% "sid_2022", yes = "2022",
                                                      no = ifelse(test = stomatal_files$origin %in% "WT_2022", yes = "2022",
                                                                  no = ifelse(test = stomatal_files$origin %in% "WT_2023-05-16", yes = "2023",
                                                                              no = ifelse(test = stomatal_files$origin %in% "WT_2023-05-17", yes = "2023", 
                                                                                          no = "2024"))))))

stomatal_files$sid_included <- ifelse(test = stomatal_files$year %in% "2021", yes = "sid included", 
                                      no = ifelse(test = stomatal_files$year %in% "2022", yes = "sid included", no = "sid not included"))
DimPlot(stomatal_files, reduction = "umap", split.by = "sid_included", group.by = "genotype", shuffle = T) +
  scale_colour_manual(values = rev(met.brewer("Hokusai3", n=2)))


DimPlot(stomatal_files, reduction = "umap", group.by = "Phase") + theme(legend.position = "bottom") +
  scale_colour_manual(values = met.brewer("Derain", n=4))


### 20.6. mark stomatal lineage cells on epidermal subset
stom <- WhichCells(stomatal_files)
DimPlot(epidermis, cells.highlight = c(stom), sizes.highlight = 0.3)












#---------------------------------------------------------------------------------------------------

###### 21. Subset guard cell lineage ----------------------------------------------------------------
Idents(stomatal_files) <- "seurat_clusters"
gc_lineage <- subset(stomatal_files, idents = c(2, 5, 6, 7, 9, 12))

dim(gc_lineage) # 1805

### 21.1. Find variable features ----------------------------------------------------------------------------------------------------------------------
gc_lineage <- FindVariableFeatures(gc_lineage, selection.method = "vst", nfeatures = 2000)

### 21.2. Scale data ----------------------------------------------------------------------------------------------------------------------------------
gc_lineage <- ScaleData(gc_lineage, features = rownames(gc_lineage))

### 21.3. run PCA ------------------------------------------------------------------------------------------------------------------------------------
gc_lineage <- RunPCA(gc_lineage, features = VariableFeatures(object = gc_lineage), npcs = 30)

ElbowPlot(gc_lineage, ndims = 30)

### 21.4. find neighbours and clusters ---------------------------------------------------------------------------------------------------------------
gc_lineage <- FindNeighbors(gc_lineage, dims = 1:25)
gc_lineage <- FindClusters(gc_lineage, resolution = 0.5) # 

### 21.5. run UMAP and visualise --------------------------------------------------------------------------------------------------------------------
gc_lineage <- RunUMAP(gc_lineage, dims = 1:25)

DimPlot(gc_lineage, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = met.brewer("Hiroshige", n=9))


gc_lineage <- RenameIdents(gc_lineage, "0" = "Stage 0-1", "1" = "Early GCs",
                           "2" = "Late GCs", "3" = "GMCs", 
                           "4" = "GMCs", "5" = "Dividing GMCs",
                           "6" = "Dividing GMCs")
gc_lineage$stages <- Idents(gc_lineage)

DimPlot(gc_lineage, reduction = "umap", group.by = "stages") +
  scale_colour_manual(values = c("#86C592", "#0A2E57", "#204877", "#3E709E", "#5D9EC1"))














#---------------------------------------------------------------------------------------------------

###### 22. Hair cell lineage ----------------------------------------------------------------
Idents(epidermis) <- "seurat_clusters"
hc_lineage <- subset(epidermis, idents = c(6, 7, 10, 12, 13, 15, 19, 11))

dim(hc_lineage) # 5032 cells

### 22.1. Find variable features ----------------------------------------------------------------------------------------------------------------------
hc_lineage <- FindVariableFeatures(hc_lineage, selection.method = "vst", nfeatures = 2000)

### 22.2. Scale data ----------------------------------------------------------------------------------------------------------------------------------
hc_lineage <- ScaleData(hc_lineage, features = rownames(hc_lineage))

### 22.3. run PCA ------------------------------------------------------------------------------------------------------------------------------------
hc_lineage <- RunPCA(hc_lineage, features = VariableFeatures(object = hc_lineage), npcs = 30)

ElbowPlot(hc_lineage, ndims = 30)

### 22.4. find neighbours and clusters ---------------------------------------------------------------------------------------------------------------
hc_lineage <- FindNeighbors(hc_lineage, dims = 1:25)
hc_lineage <- FindClusters(hc_lineage, resolution = 0.5) # 

### 22.5. run UMAP and visualise --------------------------------------------------------------------------------------------------------------------
hc_lineage <- RunUMAP(hc_lineage, dims = 1:25)

DimPlot(hc_lineage, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = met.brewer("Hiroshige", n=11))
