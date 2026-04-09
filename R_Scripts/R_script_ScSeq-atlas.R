### R script B. distachyon scSeq atlas (Berg et al. 2026)


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
whole_data <- merge(wt_str_nod, y=list(sid_str_nod, wt_strines_nod, sid_strines_nod, wt170523_str_nod, wt160523_str_nod,
                                         wt2024_1_str_nod, wt2024_2_str_nod, wt2024_3_str_nod))

whole_data <- FindVariableFeatures(whole_data, selection.method = "vst", nfeatures = 2000)

whole_data <- ScaleData(whole_data)

whole_data <- RunPCA(whole_data, features = VariableFeatures(object = whole_data), npcs = 100)

ElbowPlot(whole_data, ndims = 50)

whole_data <- FindNeighbors(whole_data, reduction = "pca", dims = 1:40)

### If you change the resolution, more or less clusters will be found but the underlying plot will look the same
whole_data <- FindClusters(whole_data, resolution = 0.8)

whole_data <- RunUMAP(whole_data, reduction = "pca", dims = 1:40)


###### 12. Visualize the dataset --------------------------------------------------------------------------------------
DimPlot(whole_data, reduction = "umap", group.by = "origin", shuffle = T)
DimPlot(whole_data, reduction = "umap", group.by = "seurat_clusters", label = T)

### add meta data that distinguishes between WT and sid (bdmute)
## for soupx-filtered data
whole_data$genotype <- ifelse(test = whole_data$origin %in% "sid_2021", yes = "bdmute", 
                                no = ifelse(test = whole_data$origin %in% "sid_2022", yes = "bdmute", no = "WT"))

DimPlot(whole_data, reduction = "umap", group.by = "genotype", shuffle = T)


### add meta data that distinguishes between different tissue datasets (2021+2022+2024 (leaf dev. zone) vs 2023 (vSAM + leaf primordia))
## for soupx-filtered data
whole_data$tissue <- ifelse(test = whole_data$origin %in% "WT_2023-05-16", yes = "SAM + leaf primordia", 
                              no = ifelse(test = whole_data$origin %in% "WT_2023-05-17", yes = "SAM + leaf primordia", 
                                          no = "leaf developmental zone"))

DimPlot(whole_data, reduction = "umap", group.by = "tissue", shuffle = T)


### Add metadata on cell cycle phases
## save g2/m cell cycle genes
g2m.genes <- c("BdiBd21-3.4G0363900", "BdiBd21-3.3G0530200","BdiBd21-3.2G0144200",
               "BdiBd21-3.2G0281100","BdiBd21-3.2G0675300","BdiBd21-3.1G0398300", 
               "BdiBd21-3.5G0241300")
s.genes <- c("BdiBd21-3.1G1030800", "BdiBd21-3.3G0027100", "BdiBd21-3.3G0090800", 
             "BdiBd21-3.3G0550700", "BdiBd21-3.4G0249100")

## Join layers to enable cell cycle scoring
whole_data <- JoinLayers(whole_data)

### give each cell a score based on cell cycle genes
whole_data <- CellCycleScoring(whole_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(whole_data, reduction = "umap", group.by = "Phase")


### save dataset
saveRDS(whole_data, "whole_data_all_genome_v1_2.rds")


### reload the data 
whole_data <- readRDS("whole_data_all_genome_v1_2.rds")




###### 13. Find out some numerical parameters -----------------------------------------------------------------------
### Number of cells
dim(whole_data)
# 69687 cells


### Mean number of UMIs
mean(whole_data$nCount_RNA) # 8205

### Median number of UMIs
median(whole_data$nCount_RNA) # 5673


### Mean number of features/genes
mean(whole_data$nFeature_RNA) # 2944

### Median number of features/genes
median(whole_data$nFeature_RNA) # 2660


### values for each of the included datasets
stats_scseq %>% group_by(ori) %>% summarise(mean_umi=mean(umi), median_umi=median(umi), 
                                     mean_fea=mean(fea), median_fea=median(fea), 
                                     cells=length(fea))



### How many genes are expressed in total?
## calculate number of cells expressing each individual gene
all_genes <- Matrix::rowSums(whole_data@assays$RNA$counts > 0)
all_genes <- as.data.frame(all_genes)
## filter out genes that are not expressed by any cell
expressed_genes <- all_genes %>% filter(all_genes>0) # 35584 out of 39068 genes are expressed in at least one cell
expressed_genes_strict <- all_genes %>% filter(all_genes>=5) # 31988 out of 39068 genes are expressed in at least 5 cells




###### 14. Visualize quality control plots -------------------------------------------------------------------------
FeaturePlot(whole_data, features = "nFeature_RNA") 
FeaturePlot(whole_data, features = "nCount_RNA") 
FeaturePlot(whole_data, features = "percent.mt") 
FeaturePlot(whole_data, features = "percent.chl")



###### 15. Tabulate cell statistics -----------------------------------------------------------------------------------
### How many cells are in each cluster?
table(whole_data$seurat_clusters)

### How many cells are in each replicate?
table(whole_data$origin)

### How does cluster membership vary by genotype?
table(Idents(whole_data), whole_data$genotype)




###### 16. Check marker gene expression on the dataset -------------------------------------------------------------
### Using a dot plot 
## if it says "-like" it means there were either several orthologues or it was not a 1:1 orthology
DotPlot(whole_data, features = c("BdiBd21-3.5G0316500" # WOX4
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
Idents(whole_data) <- "seurat_clusters"

### assign identities according to marker gene expression profiles
whole_data <- RenameIdents(whole_data, `0` = "Mesophyll", `1` = "Mesophyll", 
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

whole_data$intermediate <- Idents(whole_data)


### to merge the clusters with overlapping identities between tissues
whole_data <- RenameIdents(whole_data, "Epi.-Meso.-Vasc." = "Unknown",
                             "Epi.-Meso." = "Unknown",
                             "Meso.-Vasc." = "Unknown")

whole_data$labels <- Idents(whole_data) # save as metadata column

DimPlot(whole_data, reduction = "umap", group.by = "labels")




###### 19. Subset Vasculature and shoot apex to better annotate cell clusters ----------------------------------------
Idents(whole_data) <- "labels"

shoot_apex <- subset(whole_data, idents = c("Shoot apex", "Vasculature"))

dim(shoot_apex) # 16595 cells

### Find variable features
shoot_apex <- FindVariableFeatures(shoot_apex, selection.method = "vst", nfeatures = 2000)

### Scale data
shoot_apex <- ScaleData(shoot_apex, features = rownames(shoot_apex))

### run PCA
shoot_apex <- RunPCA(shoot_apex, features = VariableFeatures(object = shoot_apex), npcs = 30)

ElbowPlot(shoot_apex, ndims = 30)

### find neighbours and clusters
shoot_apex <- FindNeighbors(shoot_apex, dims = 1:25)
shoot_apex <- FindClusters(shoot_apex, resolution = 0.6)

### run UMAP and visualise
shoot_apex <- RunUMAP(shoot_apex, dims = 1:25)

DimPlot(shoot_apex, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = met.brewer("Hiroshige", n=16))

### check marker gene expression
Idents(shoot_apex) <- "seurat_clusters"
my_levels <- c(1, 3, 5, 7, 12, 15, 0, 2, 9, 11, 13, 6, 8, 10, 14, 4)
Idents(shoot_apex) <- factor(Idents(shoot_apex), levels = my_levels)

DotPlot(shoot_apex, features = c("BdiBd21-3.5G0316500" # WOX4
                                 , "BdiBd21-3.1G0402200" # FON1/CLV1
                                 , "BdiBd21-3.2G0427700" # FEA2
                                 , "BdiBd21-3.2G0010000" # FEA3
                                 , "BdiBd21-3.1G0571300" # FEA4
                                 , "BdiBd21-3.5G0168500" # FCP1
                                 , "BdiBd21-3.1G0135700" # KN1
                                 , "BdiBd21-3.1G0773000" # KNAT1-like1
                                 , "BdiBd21-3.1G0942300" # CRC
                                 , "BdiBd21-3.1G0588300" # PIN1a
                                 , "BdiBd21-3.2G0231400" # TED4-like2 procambium
                                 , "BdiBd21-3.1G0348000" # TMO6-like
                                 , "BdiBd21-3.1G1023300" # VND-like
                                 , "BdiBd21-3.2G0501500" # XCP1-like
                                 , "BdiBd21-3.3G0071200" # APL
                                 , "BdiBd21-3.2G0203400" # SLAH2
                                 , "BdiBd21-3.3G0204900" # PDF2-like 
                                 , "BdiBd21-3.2G0082000" # PDF1
                                 , "BdiBd21-3.2G0587500" # WOX9C-like1
                                 , "BdiBd21-3.4G0439300" # LHCA6
                                 , "BdiBd21-3.2G0749400" # STOMAGEN-1
), split.by = "annotations", # for old data, "labels"
cols = c("#5a5a83", "#8282aa", "#9d9dc7", "#c9c9dd", "#e3aba7", "#b1615c", "#d88782", "lightgrey")) +
  scale_x_discrete(labels=c("BdiBd21-3.5G0316500" = "BdWOX4",
                            "BdiBd21-3.1G0402200" = "BdCLV1", 
                            "BdiBd21-3.2G0427700" = "BdFEA2",
                            "BdiBd21-3.2G0010000" = "BdFEA3", 
                            "BdiBd21-3.1G0571300" = "BdFEA4", 
                            "BdiBd21-3.5G0168500" = "BdFCP1",
                            "BdiBd21-3.1G0135700" = "BdKN1", 
                            "BdiBd21-3.1G0773000" = "BdKNAT1-like1",
                            "BdiBd21-3.1G0942300" = "BdCRC",
                            "BdiBd21-3.1G0588300" = "BdPIN1a",
                            "BdiBd21-3.2G0231400" = "BdTED4-like2",
                            "BdiBd21-3.1G0348000" = "BdTMO6-like",
                            "BdiBd21-3.1G1023300" = "BdVND-like1",
                            "BdiBd21-3.2G0501500" = "BdXCP1-like",
                            "BdiBd21-3.3G0071200" = "BdAPL",
                            "BdiBd21-3.2G0203400" = "BdSLAH2",
                            "BdiBd21-3.3G0204900" = "BdPDF2-like",
                            "BdiBd21-3.2G0082000" = "BdPDF1",
                            "BdiBd21-3.2G0587500" = "BdWOX9C-like1",
                            "BdiBd21-3.4G0439300" = "BdLHCA6",
                            "BdiBd21-3.2G0749400" = "BdSTOMAGEN-1")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = 8)) +
  labs(y = NULL, x = NULL)

### Label clusters with tissue type
shoot_apex <- RenameIdents(shoot_apex, `0` = "Primordium", `1` = "SAM", 
                           `2` = "Primordium", `3` = "SAM", 
                           `4` = "Unclear", `5` = "SAM/Vasculature", 
                           `6` = "Vasculature", `7` = "SAM/Vasculature", 
                           `8` = "Vasculature", `9` = "Primordium", 
                           `10` = "Mesophyll", `11` = "Procambium", 
                           `12` = "SAM/Vasculature", `13` = "Procambium", 
                           `14` = "Epidermis", `15` = "SAM/Vasculature")

shoot_apex$annotations <- Idents(shoot_apex)

DimPlot(shoot_apex, group.by = "annotations") +
  scale_colour_manual(values = c("#9d9dc7", "#5a5a83", "lightgrey", "#8282aa", "#e3aba7", "#b1615c", "#c9c9dd", "#d88782"))


### Transfer cell type annotations to the whole dataset
SAM_cells <- WhichCells(shoot_apex, idents = c("SAM"))
SAM_vasc_cells <- WhichCells(shoot_apex, idents = c("SAM/Vasculature"))
vasc_cells <- WhichCells(shoot_apex, idents = c("Vasculature"))
epi_cells <- WhichCells(shoot_apex, idents = c("Epidermis"))
meso_cells <- WhichCells(shoot_apex, idents = c("Mesophyll"))
prim_cells <- WhichCells(shoot_apex, idents = c("Primordium"))
proc_cells <- WhichCells(shoot_apex, idents = c("Procambium"))
unc_cells <- WhichCells(shoot_apex, idents = c("Unclear"))

Idents(whole_data) <- "labels"

Idents(whole_data, WhichCells(object = whole_data, cells = SAM_cells, slot = 'data')) <-"SAM"
Idents(whole_data, WhichCells(object = whole_data, cells = SAM_vasc_cells, slot = 'data')) <-"SAM/Vasculature"
Idents(whole_data, WhichCells(object = whole_data, cells = vasc_cells, slot = 'data')) <-"Vasculature"
Idents(whole_data, WhichCells(object = whole_data, cells = epi_cells, slot = 'data')) <-"Epidermis"
Idents(whole_data, WhichCells(object = whole_data, cells = meso_cells, slot = 'data')) <-"Mesophyll"
Idents(whole_data, WhichCells(object = whole_data, cells = prim_cells, slot = 'data')) <-"Primordium"
Idents(whole_data, WhichCells(object = whole_data, cells = proc_cells, slot = 'data')) <-"Procambium"
Idents(whole_data, WhichCells(object = whole_data, cells = unc_cells, slot = 'data')) <-NA # labeling as NA is important so that STACAS does not recognize it as a separate cell type

whole_data <- RenameIdents(whole_data, "Unknown" = NA)

whole_data$newlabels <- Idents(whole_data)

DimPlot(whole_data) + scale_colour_manual(values = c("#c9c9dd", "#9d9dc7", "#b1615c", "#d88782", "#e3aba7", "#8282aa", "#5a5a83"))






###### 20. Semi-supervised integration with STACAS (https://pmc.ncbi.nlm.nih.gov/articles/PMC10825117/) ------------
DimPlot(whole_data, group.by = "newlabels")

whole_data <- SplitObject(whole_data, split.by = "origin")

whole_data <- Run.STACAS(whole_data, dims = 1:40, anchor.features = 2000, cell.labels = "newlabels")

whole_data <- RunUMAP(whole_data, dims = 1:40)

DefaultAssay(whole_data) <- "integrated"
whole_data <- FindNeighbors(whole_data, dims = 1:40)
whole_data <- FindClusters(whole_data, resolution = 1)

DimPlot(whole_data, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = c(
    "#F7AB58", "#E76254", "#ED8448", "#FFD47C", "#EB794C", #0,1,2,3,4
    "#3B6D98", "#F9B65F", "#5C9DB9", "#F2984E", "#F4A153", #5,6,7,8,9
    "#F9E5B9", "#E0E2C5", "#7CC2D7", "#4B84A6", "#F08F49", #10,11,12,13,14
    "#E96D50", "#5390AE", "#8DCBDA", "#65AAC5", "#FFE1A6", #15,16,17,18,19
    "#9DD4DD", "#AFDCDD", "#C8DFD1", "#FBC166", "#FFDA91", #20,21,22,23,24
    "#1E466E", "#2C5984", "#33628F", "#FECC6C", "#43799F", #25,26,27,28,29
    "#254F79", "#6EB7D1"  #30,31
  ))

DimPlot(whole_data, reduction = "umap", group.by = "origin", shuffle = T, alpha = 0.8) # how is data separated between original datasets

DimPlot(whole_data, reduction = "umap", group.by = "tissue", shuffle = T, alpha = 0.8) +
  theme(legend.position = "bottom") # how is data separated between tissue datasets

DimPlot(whole_data, reduction = "umap", group.by = "Phase")  + # show cell cycling phases
  theme(legend.position = "bottom") +
  scale_colour_manual(values = met.brewer("Derain", n = 4))


### assign tissue identity to STACAS integrated dataset
Idents(whole_data) <- "seurat_clusters"
whole_data <- RenameIdents(whole_data, `0` = "Mesophyll", `1` = "Mesophyll", 
                           `2` = "Mesophyll", `3` = "SAM/Primordium", 
                           `4` = "Mesophyll", `5` = "Epidermis", 
                           `6` = "Unclear", `7` = "Epidermis", 
                           `8` = "Mesophyll", `9` = "Mesophyll", 
                           `10` = "Vasculature", `11` = "SAM/Vasculature", 
                           `12` = "Epidermis", `13` = "Epidermis", 
                           `14` = "Mesophyll", `15` = "Mesophyll",
                           `16` = "Unclear", `17` = "Vasculature",
                           `18` = "Epidermis", `19` = "Vasculature", # 19 is Procambium
                           `20` = "Vasculature", `21` = "SAM/Vasculature",
                           `22` = "SAM/Vasculature", `23` = "Unclear",
                           `24` = "SAM/Primordium", `25` = "Epidermis",
                           `26` = "Epidermis", `27` = "Epidermis",
                           `28` = "SAM/Primordium", `29` = "Epidermis",
                           `30` = "Epidermis", `31` = "Epidermis")

whole_data$annotated <- Idents(whole_data)

Idents(whole_data) <- "annotated"

Idents(whole_data) <- factor(Idents(whole_data), levels = c("SAM/Vasculature", "SAM/Primordium", "Vasculature", "Epidermis", "Mesophyll", "Unclear"))
DimPlot(whole_data) + 
  scale_colour_manual(values = c("#454a74", "#97c684", "#efc86e", "#808fe1", "#6f9969", "lightgrey"))


###### assign tissue identity to STACAS integrated dataset with more details
Idents(whole_data) <- "seurat_clusters"
whole_data <- RenameIdents(whole_data, `0` = "Mesophyll", `1` = "Mesophyll", 
                           `2` = "Mesophyll", `3` = "SAM/Primordium", 
                           `4` = "Mesophyll", `5` = "Epidermis", 
                           `6` = "Unclear", `7` = "Epidermis", 
                           `8` = "Mesophyll", `9` = "Mesophyll", 
                           `10` = "Bundle sheath cells", `11` = "SAM/Xylem", 
                           `12` = "Epidermis", `13` = "Epidermis", 
                           `14` = "Mesophyll", `15` = "Mesophyll",
                           `16` = "Unclear", `17` = "Phloem",
                           `18` = "Epidermis", `19` = "Procambium",
                           `20` = "Phloem", `21` = "SAM/Vasculature",
                           `22` = "SAM/Vasculature", `23` = "Unclear",
                           `24` = "SAM/Primordium", `25` = "Epidermis",
                           `26` = "Epidermis", `27` = "Epidermis",
                           `28` = "SAM/Primordium", `29` = "Epidermis")

DimPlot(whole_data) +
  scale_colour_manual(values = c("#86C592", "#D8D97A", "#A8C971", "lightgrey", "#70C1C2", "#5D9EC1", "#3E709E", "#204877", "#0A2E57"))

### Quality control
par(mfrow=c(3,2))
p1 <- DimPlot(whole_data, reduction = "umap", group.by = "seurat_clusters")
p2 <- FeaturePlot(whole_data, features = "nFeature_RNA") 
p3 <- FeaturePlot(whole_data, features = "nCount_RNA") 
p4 <- FeaturePlot(whole_data, features = "percent.mt") 
p5 <- FeaturePlot(whole_data, features = "percent.chl")
p6 <- DimPlot(whole_data, reduction = "umap", group.by = "tissue", shuffle = T) + theme(legend.position = "bottom")
p1+p6+p2+p3+p4+p5

### Save object
saveRDS(whole_data, "E:/Scseq_brachy_dev/R_data_objects/whole_data_STACAS.rds")






###### 21. Link GRAS32 with mitotically active cells ------------------------------------------------------------
doi <- whole_data
doi <- AddModuleScore(doi, features = list("BdiBd21-3.1G0657800"), nbin = 2, name = "gene_expression")

### categorise cells into cells that have and don't have the genes of interest expressed
doi$goi_score <- ifelse(test = doi$gene_expression1 > 0, yes = "gene",
                        no = "not expressed")

doi$goi_score <- ifelse(test = doi$gene_expression1 > 1, yes = "gene",
                        no = "not expressed")

doi$goi_score <- ifelse(test = doi$gene_expression1 > 2, yes = "gene",
                        no = "not expressed")

doi$goi_score <- ifelse(test = doi$gene_expression1 > 2.5, yes = "gene",
                        no = "not expressed")

exp_table <- table(doi$goi_score, doi$Phase)
exp_table
## threshold expression > 0
#                 G1   G2M     S Undecided
#gene           1143  3870  1477        45
#not expressed 24793 13174 23504      1681

## threshold expression > 1
#                 G1   G2M     S Undecided
#gene            326  1721   499        40
#not expressed 25610 15323 24482      1686

# threshold expression > 2
#                 G1   G2M     S Undecided
#gene             23   314    53         7
#not expressed 25913 16730 24928      1719

#threshold expression > 2.5
#                 G1   G2M     S Undecided
#gene              5    83    10         1
#not expressed 25931 16961 24971      1725


## what percentage of cells expressing BdGRAS32 are mitotically active cells (G2M)
3870/(1143+1477+45+3870)
## --> 59% of BdGRAS-positive cells were categorised as G2M for expression threshold > 0

1721/(326+499+40+1721)
## --> 67% of BdGRAS-positive cells were categorised as G2M for expression threshold > 1

314/(23+53+7+314)
## --> 79% of BdGRAS-positive cells were categorised as G2M for expression threshold > 2

83/(5+10+1+83)
## --> 84% of BdGRAS-positive cells were categorised as G2M for expression threshold > 2.5








###### 22. Identify markers ----------------------------------------------------------------------------------------
## Differentially expressed genes in mitotically active cells (G2M phase)
Idents(whole_data) <- "Phase"
g2mphase.markers <- FindMarkers(whole_data, ident.1 = "G2M")

# filter for significant marker that are upregulated in G2M cells
g2mphase.markers <- filter(g2mphase.markers, p_val_adj < 0.5 & avg_log2FC > 0)

# filter for transcription factors
TFs <- read.table("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/TF_list_brachy_all.tsv", sep = "\t")
TFs <- as.vector(TFs$V1)

g2mphase_tfs <- filter(g2mphase.markers, rownames(g2mphase.markers) %in% TFs)

write.table(g2mphase_tfs, "E:/Scseq_brachy_dev/markers/revisions_paper/all_G2M_TFs.txt", sep = ";", col.names = T)

gras32 <- filter(g2mphase_tfs, rownames(g2mphase_tfs) == "BdiBd21-3.1G0657800")


## Differentially expressed genes in all clusters
setwd("E:/Scseq_brachy_dev/markers/revisions_paper/")
Idents(whole_data) <- "seurat_clusters"
all.markers <- FindAllMarkers(whole_data, min.pct = 0.25)


celltypes <- data.frame(cluster = c(0,1:31), celltype = c("Mesophyll", "Mesophyll", 
                                                          "Mesophyll", "SAM/Primordium", 
                                                          "Mesophyll", "Epidermis", 
                                                          "Unclear", "Epidermis", 
                                                          "Mesophyll", "Mesophyll", 
                                                          "Vasculature", "SAM/Vasculature", 
                                                          "Epidermis", "Epidermis", 
                                                          "Mesophyll", "Mesophyll",
                                                          "Unclear", "Vasculature",
                                                          "Epidermis", "Vasculature", 
                                                          "Vasculature", "SAM/Vasculature",
                                                          "SAM/Vasculature", "Unclear",
                                                          "SAM/Primordium", "Epidermis",
                                                          "Epidermis", "Epidermis",
                                                          "SAM/Primordium", "Epidermis",
                                                          "Epidermis", "Epidermis"))

all_markers_info <- merge(all.markers, celltypes, by="cluster", all.x=T)


all_markers_info_strict <- filter(all_markers_info, p_val_adj < 0.5)




write.table(all_markers_info_strict[,c(7,1,8,2:6)], "E:/Scseq_brachy_dev/markers/revisions_paper/markers_all_clusters.txt", sep=";", row.names = F)








#---------------------------------------------------------------------------------------------------------------------

###### 23. Subset epidermal dataset ---------------------------------------------------------------------------------
DefaultAssay(whole_data) <- "integrated"
Idents(whole_data) <- "annotated"

epidermis <- subset(whole_data, idents = c("Epidermis"))

dim(epidermis)
## contains 15034 cells

### Scale data
epidermis <- ScaleData(epidermis, features = rownames(epidermis))

### run PCA
epidermis <- RunPCA(epidermis, features = VariableFeatures(object = epidermis), npcs = 30)

ElbowPlot(epidermis, ndims = 30)

### find neighbours and clusters 
epidermis <- FindNeighbors(epidermis, dims = 1:25)
epidermis <- FindClusters(epidermis, resolution = 1)

### run UMAP and visualise 
epidermis <- RunUMAP(epidermis, dims = 1:25)

DimPlot(epidermis, reduction = "umap", group.by = "seurat_clusters", label = T) + 
  scale_colour_manual(values = met.brewer("Hiroshige", n=25))

DimPlot(epidermis, reduction = "umap", group.by = "Phase") + theme(legend.position = "bottom") +
  scale_colour_manual(values = met.brewer("Derain", n=4))

DimPlot(epidermis, reduction = "umap", group.by = "tissue", order = T) + theme(legend.position = "bottom")


### For paper: mark epidermal cells on whole dataset
epi_cells <- WhichCells(epidermis)

test <- whole_data
Idents(test, WhichCells(object = test, cells = epi_cells, slot = 'data')) <-"Epidermal subset"

DimPlot(test, cells.highlight = c(epi_cells), sizes.highlight = 0.3)

### mark HC lineage on epidermal dataset
HCs <- WhichCells(epidermis, idents = c("Stage 0-2", "Hair cell lineage"))

test <- epidermis
Idents(test, WhichCells(object = test, cells = HCs, slot = 'data')) <-"Hair cell lineage"

DimPlot(test, cells.highlight = c(HCs), sizes.highlight = 0.3)


### mark stomatal lineage on epidermal dataset
sto <- WhichCells(epidermis, idents = c("Stage 0-2", "Inter cells", "Stomatal lineage"))

test <- epidermis
Idents(test, WhichCells(object = test, cells = sto, slot = 'data')) <-"Stomatal lineage"

DimPlot(test, cells.highlight = c(sto), sizes.highlight = 0.3)


### Quality control
par(mfrow=c(3,2))
p1 <- DimPlot(epidermis, reduction = "umap")
p6 <- DimPlot(epidermis, reduction = "umap", group.by = "tissue") + theme(legend.position = "bottom")
p2 <- FeaturePlot(epidermis, features = "nFeature_RNA") 
p3 <- FeaturePlot(epidermis, features = "nCount_RNA") 
p4 <- FeaturePlot(epidermis, features = "percent.mt") 
p5 <- FeaturePlot(epidermis, features = "percent.chl")
p1+p6+p2+p3+p4+p5


### Annotate clusters
epidermis <- RenameIdents(epidermis, `0` = "Silica cell lineage", `1` = "Stage 0-2", 
                          `2` = "Hair cell lineage", `3` = "Stage 0-2", 
                          `4` = "Inter cells", `5` = "Stage 0-2", 
                          `6` = "Silica cell lineage", `7` = "Inter cells", 
                          `8` = "Hair cell lineage", `9` = "Silica cell lineage", 
                          `10` = "Stomatal lineage", `11` = "Stage 0-2", 
                          `12` = "Stage 0-2", `13` = "Stomatal lineage",
                          `14` = "Hair cell lineage", `15` = "Unclear",
                          `16` = "Unclear", `17` = "Inter cells",
                          `18` = "Hair cell lineage", `19` = "Stomatal lineage",
                          `20` = "Stomatal lineage", `21` = "Silica cell lineage",
                          `22` = "Stomatal lineage", `23` = "Hair cell lineage",
                          `24` = "Hair cell lineage")

epidermis$celltypes <- Idents(epidermis)

DimPlot(epidermis, label = F, shuffle = T) +
  scale_colour_manual(values = c("lightgrey", "#c1d1aa", "#7fa074", "#574571", "#2c4b27", "#b695bc")) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3))


### Find markers for every epidermis cluster compared to all remaining cells, report only the positive ones:
Idents(epidermis) <- "seurat_clusters"
epidermis.markers <- FindAllMarkers(epidermis, min.pct = 0.25)

celltypes <- data.frame(cluster = c(0,1:24), celltype = c("Silica cell lineage", "Stage 0-2", 
                                                          "Hair cell lineage", "Stage 0-2", 
                                                          "Stage 0-2", "Stage 0-2", 
                                                          "Silica cell lineage", "Inter-specialized cells", 
                                                          "Hair cell lineage", "Silica cell lineage", 
                                                          "Stomatal lineage", "Stage 0-2", 
                                                          "Stage 0-2", "Stomatal lineage",
                                                          "Hair cell lineage", "Unclear",
                                                          "Unclear", "Stage 0-2",
                                                          "Hair cell lineage", "Stomatal lineage",
                                                          "Stomatal lineage", "Silica cell lineage",
                                                          "Stomatal lineage", "Hair cell lineage",
                                                          "Hair cell lineage"))

epidermis_markers_info <- merge(epidermis.markers, celltypes, by="cluster", all.x=T)

epidermis_markers_info_strict <- filter(epidermis_markers_info, p_val_adj < 0.5)

write.table(epidermis_markers_info_strict[,c(7,1,8,2:6)], "E:/Scseq_brachy_dev/markers/revisions_paper/markers_epidermis.txt", sep=";", row.names = F)















#-------------------------------------------------------------------------------------------------

###### 24. Subset stomatal lineage subset ---------------------------------------------------------
Idents(epidermis) <- "celltypes"
stomatal_files <- subset(epidermis, idents = c("Stomatal lineage", "Inter cells", "Stage 0-2"))

dim(stomatal_files)
## 8406

DefaultAssay(stomatal_files) <- "integrated"

### Scale data
stomatal_files <- ScaleData(stomatal_files, features = rownames(stomatal_files))

### run PCA
stomatal_files <- RunPCA(stomatal_files, features = VariableFeatures(object = stomatal_files), npcs = 30)

ElbowPlot(stomatal_files, ndims = 30)

### find neighbours and clusters
stomatal_files <- FindNeighbors(stomatal_files, dims = 1:25)
stomatal_files <- FindClusters(stomatal_files, resolution = 0.8)

### run UMAP and visualise
stomatal_files <- RunUMAP(stomatal_files, dims = 1:25)

DimPlot(stomatal_files, reduction = "umap", group.by = "seurat_clusters", label = T) + 
  scale_colour_manual(values = met.brewer("Hiroshige", n=17))

### label datasets according to whether that year had a sid dataset (2021, 2022) or not (2023, 2024)
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

### mark stomatal lineage cells on epidermal subset
stom <- WhichCells(stomatal_files)
DimPlot(epidermis, cells.highlight = c(stom), sizes.highlight = 0.3)

### mark guard cell lineage on stomatal dataset
gc_cells <- WhichCells(gc_lineage)

DimPlot(stomatal_files, cells.highlight = c(gc_cells), sizes.highlight = 0.3)


### Annotate clusters
stomatal_files <- RenameIdents(stomatal_files, `0` = "Stage 0-2", `1` = "Stage 0-2", 
                               `2` = "Stage 0-2", `3` = "Inter-specialized cells", 
                               `4` = "Stage 0-2", `5` = "SCs", 
                               `6` = "SMCs", `7` = "Early GCs", 
                               `8` = "SMCs", `9` = "Early GMCs", 
                               `10` = "Inter-specialized cells", `11` = "Stage 0-2", 
                               `12` = "Late GCs", `13` = "Stage 0-2",
                               `14` = "Dividing GMCs", `15` = "Stage 0-2",
                               `16` = "Dividing GMCs")

stomatal_files$stages <- Idents(stomatal_files)

DimPlot(stomatal_files, label = F, shuffle = T, order = rev(c("Stage 0-2", "Inter-specialized cells",  
                                                              "SMCs", "SCs", "Early GMCs", 
                                                              "Dividing GMCs", "Early GCs", "Late GCs"))) +
  scale_colour_manual(values = c("#78C7B8", "#67AFC2", "#8CC483", 
                                 "#D8D97A", "#538EB9",
                                 "#356493", "#1D4573", "#0A2E57")) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 4))


## without interstomatal and SMC labels
DimPlot(stomatal_files, label = F, shuffle = T, order = rev(c("Stage 0-2", "Inter-specialized cells",  
                                                              "SMCs", "SCs", "Early GMCs", 
                                                              "Dividing GMCs", "Early GCs", "Late GCs"))) +
  scale_colour_manual(values = c("#78C7B8", "lightgrey", "lightgrey", 
                                 "#D8D97A", "#538EB9",
                                 "#356493", "#1D4573", "#0A2E57")) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 4))


### Quality control
par(mfrow=c(3,2))
p1 <- DimPlot(stomatal_files, reduction = "umap")
p6 <- DimPlot(stomatal_files, reduction = "umap", group.by = "tissue") + theme(legend.position = "bottom")
p2 <- FeaturePlot(stomatal_files, features = "nFeature_RNA") 
p3 <- FeaturePlot(stomatal_files, features = "nCount_RNA") 
p4 <- FeaturePlot(stomatal_files, features = "percent.mt") 
p5 <- FeaturePlot(stomatal_files, features = "percent.chl")
p1+p6+p2+p3+p4+p5


### Find markers
Idents(stomatal_files) <- "seurat_clusters"
stomatal.markers <- FindAllMarkers(stomatal_files, min.pct = 0.1)

celltypes <- data.frame(cluster = c(0,1:16), celltype = c("Stage 0-2", "Stage 0-2", 
                                                          "Stage 0-2", "Inter-specialized cells", 
                                                          "Stage 0-2", "SCs", 
                                                          "SMCs", "Early GCs", 
                                                          "SMCs", "Early GMCs", 
                                                          "Inter-specialized cells", "Stage 0-2", 
                                                          "Late GCs", "Stage 0-2",
                                                          "Dividing GMCs", "Stage 0-2",
                                                          "Dividing GMCs"))

stomatal_markers_info <- merge(stomatal.markers, celltypes, by="cluster", all.x=T)

write.table(stomatal_markers_info[,c(7,1,8,2:6)], "E:/Scseq_brachy_dev/markers/revisions_paper/markers_stomatal_files.txt", sep=";", row.names = F)










#---------------------------------------------------------------------------------------------------

###### 25. Subset guard cell lineage ----------------------------------------------------------------
Idents(stomatal_files) <- "seurat_clusters"
gc_lineage <- subset(stomatal_files, idents = c(2,11,13,9,14,16,7,12))

dim(gc_lineage)
## 2692

DefaultAssay(gc_lineage) <- "integrated"

### Scale data
gc_lineage <- ScaleData(gc_lineage, features = rownames(gc_lineage))

### run PCA
gc_lineage <- RunPCA(gc_lineage, features = VariableFeatures(object = gc_lineage), npcs = 30)

ElbowPlot(gc_lineage, ndims = 30)

### find neighbours and clusters
gc_lineage <- FindNeighbors(gc_lineage, dims = 1:25)
gc_lineage <- FindClusters(gc_lineage, resolution = 0.5) # 

### run UMAP and visualise
gc_lineage <- RunUMAP(gc_lineage, dims = 1:25)

DimPlot(gc_lineage, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = met.brewer("Hiroshige", n=9))


DimPlot(gc_lineage, reduction = "umap", group.by = "stages") +
  scale_colour_manual(values = c("#78C7B8", "#1D4573", "#538EB9", "#0A2E57", "#356493"))















#---------------------------------------------------------------------------------------------------

###### 26. Hair cell lineage ----------------------------------------------------------------
### HC lineage plus stage 0-2
Idents(epidermis) <- "celltypes"
hc_lineage <- subset(epidermis, idents = c("Stage 0-2", "Hair cell lineage"))

dim(hc_lineage)
## 8467

DefaultAssay(hc_lineage) <- "integrated"

### Scale data
hc_lineage <- ScaleData(hc_lineage, features = rownames(hc_lineage))

### run PCA
hc_lineage <- RunPCA(hc_lineage, features = VariableFeatures(object = hc_lineage), npcs = 30)

ElbowPlot(hc_lineage, ndims = 30)

### find neighbours and clusters
hc_lineage <- FindNeighbors(hc_lineage, dims = 1:25)
hc_lineage <- FindClusters(hc_lineage, resolution = 0.5) # 

### run UMAP and visualise
hc_lineage <- RunUMAP(hc_lineage, dims = 1:25)

DimPlot(hc_lineage, reduction = "umap", group.by = "seurat_clusters", label = T) +
  scale_colour_manual(values = met.brewer("Hiroshige", n=14))
