### Plot features onto single cell dataset UMAP
library(tidyverse)
setwd("E:/Scseq_brachy_dev/")

###### ggPlantmap files -----------------------------------------------------------------------------------------------
#devtools::install_github("leonardojo/ggPlantmap")

library(ggPlantmap)
library(tidyverse)
library(Seurat)
library(MetBrewer)

### Create Plantmap template
# file was created in icy following the tutorial on the ggplantmap github pages
dev.plantmap <- XML.to.ggPlantmap("E:/Scseq_brachy_dev/ggPlantmap/expression_browser1.xml")

dev.plantmap$ROI.name <- factor(dev.plantmap$ROI.name, levels = c("Shoot_apex", "Mesophyll",
                                                                  "Vasculature", "Epidermis",
                                                                  "Early_GMC", "Div_GMC",
                                                                  "Early_GC", "Late_GC",
                                                                  "Early_HC", "Middle_HC", 
                                                                  "Late_HC", "SC",
                                                                  "Silica_cell", "Interstomatal_cell"))

for(row in 1:nrow(dev.plantmap)) {
  if(dev.plantmap$ROI.name[row] == "Interstomatal_cell") {
    dev.plantmap$ROI.name[row] <- "Inter_specialized_cell"
  }
}


### Show the map
ggPlantmap.plot(dev.plantmap) + 
  guides(fill=guide_legend(title = "Cell type/stage")) +
  scale_fill_manual(values = met.brewer("Hiroshige", n = 15))

### save as R object
saveRDS(dev.plantmap, "leaf_dev_plantmap_revised.rds")




###### Generate average expression file for ggPlantmap ----------------------------------------------------------------
whole_data <- readRDS("E:/Scseq_brachy_dev/R_data_objects/whole_data_STACAS.rds")
epidermis <- readRDS("E:/Scseq_brachy_dev/R_data_objects/epidermis_withoutprim_STACAS.rds")
stomatal_files <- readRDS("E:/Scseq_brachy_dev/R_data_objects/stomatal_cell_files_withoutprim_STACAS")


### prepare data for tissue types (Shoot apex, Mesophyll, Vasculature, Epidermis)
DimPlot(whole_data)

Idents(whole_data) <- "seurat_clusters"
whole_data <- RenameIdents(whole_data, `0` = "Mesophyll", `1` = "Mesophyll", 
                           `2` = "Mesophyll", `3` = "Shoot_apex", 
                           `4` = "Mesophyll", `5` = "Epidermis", 
                           `6` = "Unclear", `7` = "Epidermis", 
                           `8` = "Mesophyll", `9` = "Mesophyll", 
                           `10` = "Vasculature", `11` = "Shoot_apex", 
                           `12` = "Epidermis", `13` = "Epidermis", 
                           `14` = "Mesophyll", `15` = "Mesophyll",
                           `16` = "Unclear", `17` = "Vasculature",
                           `18` = "Epidermis", `19` = "Vasculature",
                           `20` = "Vasculature", `21` = "Shoot_apex",
                           `22` = "Shoot_apex", `23` = "Unclear",
                           `24` = "Shoot_apex", `25` = "Epidermis",
                           `26` = "Epidermis", `27` = "Epidermis",
                           `28` = "Shoot_apex", `29` = "Epidermis",
                           `30` = "Epidermis", `31` = "Epidermis")

alltissues <- subset(whole_data, idents = c("Shoot_apex", "Mesophyll", "Vasculature", "Epidermis"))
DimPlot(alltissues)

## average gene expression per tissue identity
all2 <- as.data.frame(AverageExpression(alltissues, assays = "RNA"))
colnames(all2)
colnames(all2) <- c("Mesophyll", "Shoot_apex", "Epidermis", "Vasculature")
all2$gene <- rownames(all2)





### prepare data for epidermal cell types (Silica cells, early, middle, late HCs and inter-specialized cells)
DimPlot(epidermis, reduction = "umap", group.by = "seurat_clusters", label = T)

Idents(epidermis) <- "seurat_clusters"
epi <- subset(epidermis, idents = c(0, 2, 4, 6, 7, 8, 9, 14, 17, 21))

Idents(epi) <- "seurat_clusters"

epi <- RenameIdents(epi, `0` = "Silica_cell",
                    `2` = "Early_HC",
                    `4` = "Inter_specialized_cell",
                    `6` = "Silica_cell", 
                    `7` = "Inter_specialized_cell",
                    `8` = "Middle_HC", 
                    `9` = "Silica_cell",
                    `14` = "Late_HC",
                    `17` = "Inter_specialized_cell",
                    `21` = "Silica_cell")

DimPlot(epi)

## reorder cluster names to display legend in better order
epi$labels <- Idents(epi)
epi$labels <- factor(epi$labels, levels = c("Early_HC", "Middle_HC", "Late_HC", "Silica_cell", "Inter_specialized_cell"))

DimPlot(epi)
DimPlot(epi, group.by = "labels") + 
  scale_colour_manual(values = met.brewer("Hiroshige", n=5))

Idents(epi) <- "labels"

## average gene expression per cell identity
epi2 <- as.data.frame(AverageExpression(epi, assays = "RNA"))
colnames(epi2) <- c("Early_HC", "Middle_HC", "Late_HC", "Silica_cell", "Inter_specialized_cell")
epi2$gene <- rownames(epi2)





### prepare data for stomatal cell types and stages (Early_GMC, Div_GMC, Early_GC, Late_GC, SC)
DimPlot(stomatal_files)
Idents(stomatal_files) <- "seurat_clusters"
sto <- subset(stomatal_files, idents = c(5, 7, 9, 12, 14, 16))

sto <- RenameIdents(sto, `5` = "SC",
                    `7` = "Early_GC",  
                    `9` = "Early_GMC", 
                    `12` = "Late_GC", 
                    `14` = "Div_GMC",
                    `16` = "Div_GMC")

DimPlot(sto)

## average gene expression per cell identity
sto2 <- as.data.frame(AverageExpression(sto, assays = "RNA"))
colnames(sto2)
colnames(sto2) <- c("SC", "Early_GC", "Early_GMC", "Late_GC", "Div_GMC")
sto2$gene <- rownames(sto2)






### merge data and prepare final table
combined_data <- merge(all2, c(epi2, sto2), by = "gene")
combined_data <- select(combined_data, -gene.1)

combined_data2 <- combined_data %>%
  pivot_longer(-gene) ## making it tidy

write.table(combined_data2, "E:/Scseq_brachy_dev/ggPlantmap/averaged_expression_per_celltype_STACAS.csv")










###### Gene expression files -------------------------------------------------------------------------------------------------------
### Table to indicate file that needs to be loaded for the gene
features <- rownames(whole_data)
location_file <- data.frame(accession = features, all_cells1 = NA, all_cells2 = NA, epidermis = NA, stomatal_cell_files = NA)
rownames(location_file) <- location_file$accession

location_file[1:1197, "all_cells1"] <- "www/All_features_part1.arrow"
location_file[1198:2397, "all_cells1"] <- "www/All_features_part2.arrow"
location_file[2398:3597, "all_cells1"] <- "www/All_features_part3.arrow"
location_file[3598:4797, "all_cells1"] <- "www/All_features_part4.arrow"
location_file[4798:5997, "all_cells1"] <- "www/All_features_part5.arrow"
location_file[5998:7197, "all_cells1"] <- "www/All_features_part6.arrow"
location_file[7198:8397, "all_cells1"] <- "www/All_features_part7.arrow"
location_file[8398:9597, "all_cells1"] <- "www/All_features_part8.arrow"
location_file[9598:10797, "all_cells1"] <- "www/All_features_part9.arrow"
location_file[10798:11997, "all_cells1"] <- "www/All_features_part10.arrow"
location_file[11998:13197, "all_cells1"] <- "www/All_features_part11.arrow"
location_file[13198:14397, "all_cells1"] <- "www/All_features_part12.arrow"
location_file[14398:15597, "all_cells1"] <- "www/All_features_part13.arrow"
location_file[15598:16797, "all_cells1"] <- "www/All_features_part14.arrow"
location_file[16798:17997, "all_cells1"] <- "www/All_features_part15.arrow"
location_file[17998:19197, "all_cells1"] <- "www/All_features_part16.arrow"
location_file[19198:20397, "all_cells1"] <- "www/All_features_part17.arrow"
location_file[20398:21597, "all_cells1"] <- "www/All_features_part18.arrow"
location_file[21598:22797, "all_cells1"] <- "www/All_features_part19.arrow"
location_file[22798:23997, "all_cells1"] <- "www/All_features_part20.arrow"
location_file[23998:25197, "all_cells1"] <- "www/All_features_part21.arrow"
location_file[25198:26397, "all_cells1"] <- "www/All_features_part22.arrow"
location_file[26398:27597, "all_cells1"] <- "www/All_features_part23.arrow"
location_file[27598:28797, "all_cells1"] <- "www/All_features_part24.arrow"
location_file[28798:29997, "all_cells1"] <- "www/All_features_part25.arrow"
location_file[29998:32397, "all_cells1"] <- "www/All_features_part26.arrow"
location_file[32398:33997, "all_cells1"] <- "www/All_features_part27.arrow"
location_file[33998:35797, "all_cells1"] <- "www/All_features_part28.arrow"
location_file[35798:37597, "all_cells1"] <- "www/All_features_part29.arrow"
location_file[37598:39068, "all_cells1"] <- "www/All_features_part30.arrow"

location_file[1:1197, "all_cells2"] <- "www/All_features_part31.arrow"
location_file[1198:2397, "all_cells2"] <- "www/All_features_part32.arrow"
location_file[2398:3597, "all_cells2"] <- "www/All_features_part33.arrow"
location_file[3598:4797, "all_cells2"] <- "www/All_features_part34.arrow"
location_file[4798:5997, "all_cells2"] <- "www/All_features_part35.arrow"
location_file[5998:7197, "all_cells2"] <- "www/All_features_part36.arrow"
location_file[7198:8397, "all_cells2"] <- "www/All_features_part37.arrow"
location_file[8398:9597, "all_cells2"] <- "www/All_features_part38.arrow"
location_file[9598:10797, "all_cells2"] <- "www/All_features_part39.arrow"
location_file[10798:11997, "all_cells2"] <- "www/All_features_part40.arrow"
location_file[11998:13197, "all_cells2"] <- "www/All_features_part41.arrow"
location_file[13198:14397, "all_cells2"] <- "www/All_features_part42.arrow"
location_file[14398:15597, "all_cells2"] <- "www/All_features_part43.arrow"
location_file[15598:16797, "all_cells2"] <- "www/All_features_part44.arrow"
location_file[16798:17997, "all_cells2"] <- "www/All_features_part45.arrow"
location_file[17998:19197, "all_cells2"] <- "www/All_features_part46.arrow"
location_file[19198:20397, "all_cells2"] <- "www/All_features_part47.arrow"
location_file[20398:21597, "all_cells2"] <- "www/All_features_part48.arrow"
location_file[21598:22797, "all_cells2"] <- "www/All_features_part49.arrow"
location_file[22798:23997, "all_cells2"] <- "www/All_features_part50.arrow"
location_file[23998:25197, "all_cells2"] <- "www/All_features_part51.arrow"
location_file[25198:26397, "all_cells2"] <- "www/All_features_part52.arrow"
location_file[26398:27597, "all_cells2"] <- "www/All_features_part53.arrow"
location_file[27598:28797, "all_cells2"] <- "www/All_features_part54.arrow"
location_file[28798:29997, "all_cells2"] <- "www/All_features_part55.arrow"
location_file[29998:32397, "all_cells2"] <- "www/All_features_part56.arrow"
location_file[32398:33997, "all_cells2"] <- "www/All_features_part57.arrow"
location_file[33998:35797, "all_cells2"] <- "www/All_features_part58.arrow"
location_file[35798:37597, "all_cells2"] <- "www/All_features_part59.arrow"
location_file[37598:39068, "all_cells2"] <- "www/All_features_part60.arrow"

location_file[1:2897, "epidermis"] <- "www/Epi_part1.arrow"
location_file[2898:6597, "epidermis"] <- "www/Epi_part2.arrow"
location_file[6598:9997, "epidermis"] <- "www/Epi_part3.arrow"
location_file[9998:13397, "epidermis"] <- "www/Epi_part4.arrow"
location_file[13398:16997, "epidermis"] <- "www/Epi_part5.arrow"
location_file[16998:20197, "epidermis"] <- "www/Epi_part6.arrow"
location_file[20198:23897, "epidermis"] <- "www/Epi_part7.arrow"
location_file[23898:26997, "epidermis"] <- "www/Epi_part8.arrow"
location_file[26998:30497, "epidermis"] <- "www/Epi_part9.arrow"
location_file[30498:34797, "epidermis"] <- "www/Epi_part10.arrow"
location_file[34798:39068, "epidermis"] <- "www/Epi_part11.arrow"

location_file[1:6997, "stomatal_cell_files"] <- "www/Stom_part1.arrow"
location_file[6998:13997, "stomatal_cell_files"] <- "www/Stom_part2.arrow"
location_file[13998:20997, "stomatal_cell_files"] <- "www/Stom_part3.arrow"
location_file[20998:27997, "stomatal_cell_files"] <- "www/Stom_part4.arrow"
location_file[27998:33997, "stomatal_cell_files"] <- "www/Stom_part5.arrow"
location_file[33998:39068, "stomatal_cell_files"] <- "www/Stom_part6.arrow"


write_delim(location_file, "E:/Scseq_brachy_dev/Files for webpage/Location_file.txt", delim = "\t")
saveRDS(location_file, "E:/Scseq_brachy_dev/Files for webpage/Location_file_arrowfiles.rds")



#### save as arrow files instead of rds files -----------------------------------------------------------------
library(arrow)

### For all cells
features <- rownames(whole_data)
gene_info_1 <- FetchData(whole_data, vars = c(c("umap_1", "umap_2"), features), cells = 1:35000)
gene_info_2 <- FetchData(whole_data, vars = c(c("umap_1", "umap_2"), features), cells = 35001:69687)


gene_info_part1 <- gene_info_1[,1:1200]
gene_info_part1$cell <- rownames(gene_info_part1)
write_feather(gene_info_part1, sink = "D:/Github/website/www/All_features_part1.arrow")

gene_info_part2 <- gene_info_1[,1201:2400]
gene_info_part2$cell <- rownames(gene_info_part2)
write_feather(gene_info_part2, sink = "D:/Github/website/www/All_features_part2.arrow")

gene_info_part3 <- gene_info_1[,2401:3600]
gene_info_part3$cell <- rownames(gene_info_part3)
write_feather(gene_info_part3, sink = "D:/Github/website/www/All_features_part3.arrow")

gene_info_part4 <- gene_info_1[,3601:4800]
gene_info_part4$cell <- rownames(gene_info_part4)
write_feather(gene_info_part4, sink = "D:/Github/website/www/All_features_part4.arrow")

gene_info_part5 <- gene_info_1[,4801:6000]
gene_info_part5$cell <- rownames(gene_info_part5)
write_feather(gene_info_part5, sink = "D:/Github/website/www/All_features_part5.arrow")

gene_info_part6 <- gene_info_1[,6001:7200]
gene_info_part6$cell <- rownames(gene_info_part6)
write_feather(gene_info_part6, sink = "D:/Github/website/www/All_features_part6.arrow")

gene_info_part7 <- gene_info_1[,7201:8400]
gene_info_part7$cell <- rownames(gene_info_part7)
write_feather(gene_info_part7, sink = "D:/Github/website/www/All_features_part7.arrow")

gene_info_part8 <- gene_info_1[,8401:9600]
gene_info_part8$cell <- rownames(gene_info_part8)
write_feather(gene_info_part8, sink = "D:/Github/website/www/All_features_part8.arrow")

gene_info_part9 <- gene_info_1[,9601:10800]
gene_info_part9$cell <- rownames(gene_info_part9)
write_feather(gene_info_part9, sink = "D:/Github/website/www/All_features_part9.arrow")

gene_info_part10 <- gene_info_1[,10801:12000]
gene_info_part10$cell <- rownames(gene_info_part10)
write_feather(gene_info_part10, sink = "D:/Github/website/www/All_features_part10.arrow")

gene_info_part11 <- gene_info_1[,12001:13200]
gene_info_part11$cell <- rownames(gene_info_part11)
write_feather(gene_info_part11, sink = "D:/Github/website/www/All_features_part11.arrow")

gene_info_part12 <- gene_info_1[,13201:14400]
gene_info_part12$cell <- rownames(gene_info_part12)
write_feather(gene_info_part12, sink = "D:/Github/website/www/All_features_part12.arrow")

gene_info_part13 <- gene_info_1[,14401:15600]
gene_info_part13$cell <- rownames(gene_info_part13)
write_feather(gene_info_part13, sink = "D:/Github/website/www/All_features_part13.arrow")

gene_info_part14 <- gene_info_1[,15601:16800]
gene_info_part14$cell <- rownames(gene_info_part14)
write_feather(gene_info_part14, sink = "D:/Github/website/www/All_features_part14.arrow")

gene_info_part15 <- gene_info_1[,16801:18000]
gene_info_part15$cell <- rownames(gene_info_part15)
write_feather(gene_info_part15, sink = "D:/Github/website/www/All_features_part15.arrow")

gene_info_part16 <- gene_info_1[,18001:19200]
gene_info_part16$cell <- rownames(gene_info_part16)
write_feather(gene_info_part16, sink = "D:/Github/website/www/All_features_part16.arrow")

gene_info_part17 <- gene_info_1[,19201:20400]
gene_info_part17$cell <- rownames(gene_info_part17)
write_feather(gene_info_part17, sink = "D:/Github/website/www/All_features_part17.arrow")

gene_info_part18 <- gene_info_1[,20401:21600]
gene_info_part18$cell <- rownames(gene_info_part18)
write_feather(gene_info_part18, sink = "D:/Github/website/www/All_features_part18.arrow")

gene_info_part19 <- gene_info_1[,21601:22800]
gene_info_part19$cell <- rownames(gene_info_part19)
write_feather(gene_info_part19, sink = "D:/Github/website/www/All_features_part19.arrow")

gene_info_part20 <- gene_info_1[,22801:24000]
gene_info_part20$cell <- rownames(gene_info_part20)
write_feather(gene_info_part20, sink = "D:/Github/website/www/All_features_part20.arrow")

gene_info_part21 <- gene_info_1[,24001:25200]
gene_info_part21$cell <- rownames(gene_info_part21)
write_feather(gene_info_part21, sink = "D:/Github/website/www/All_features_part21.arrow")

gene_info_part22 <- gene_info_1[,25201:26400]
gene_info_part22$cell <- rownames(gene_info_part22)
write_feather(gene_info_part22, sink = "D:/Github/website/www/All_features_part22.arrow")

gene_info_part23 <- gene_info_1[,26401:27600]
gene_info_part23$cell <- rownames(gene_info_part23)
write_feather(gene_info_part23, sink = "D:/Github/website/www/All_features_part23.arrow")

gene_info_part24 <- gene_info_1[,27601:28800]
gene_info_part24$cell <- rownames(gene_info_part24)
write_feather(gene_info_part24, sink = "D:/Github/website/www/All_features_part24.arrow")

gene_info_part25 <- gene_info_1[,28801:30000]
gene_info_part25$cell <- rownames(gene_info_part25)
write_feather(gene_info_part25, sink = "D:/Github/website/www/All_features_part25.arrow")

gene_info_part26 <- gene_info_1[,30001:32400]
gene_info_part26$cell <- rownames(gene_info_part26)
write_feather(gene_info_part26, sink = "D:/Github/website/www/All_features_part26.arrow")

gene_info_part27 <- gene_info_1[,32401:34000]
gene_info_part27$cell <- rownames(gene_info_part27)
write_feather(gene_info_part27, sink = "D:/Github/website/www/All_features_part27.arrow")

gene_info_part28 <- gene_info_1[,34001:35800]
gene_info_part28$cell <- rownames(gene_info_part28)
write_feather(gene_info_part28, sink = "D:/Github/website/www/All_features_part28.arrow")

gene_info_part29 <- gene_info_1[,35801:37600]
gene_info_part29$cell <- rownames(gene_info_part29)
write_feather(gene_info_part29, sink = "D:/Github/website/www/All_features_part29.arrow")

gene_info_part30 <- gene_info_1[,37601:39070]
gene_info_part30$cell <- rownames(gene_info_part30)
write_feather(gene_info_part30, sink = "D:/Github/website/www/All_features_part30.arrow")

gene_info_part31 <- gene_info_2[,1:1200]
gene_info_part31$cell <- rownames(gene_info_part31)
write_feather(gene_info_part31, sink = "D:/Github/website/www/All_features_part31.arrow")

gene_info_part32 <- gene_info_2[,1201:2400]
gene_info_part32$cell <- rownames(gene_info_part32)
write_feather(gene_info_part32, sink = "D:/Github/website/www/All_features_part32.arrow")

gene_info_part33 <- gene_info_2[,2401:3600]
gene_info_part33$cell <- rownames(gene_info_part33)
write_feather(gene_info_part33, sink = "D:/Github/website/www/All_features_part33.arrow")

gene_info_part34 <- gene_info_2[,3601:4800]
gene_info_part34$cell <- rownames(gene_info_part34)
write_feather(gene_info_part34, sink = "D:/Github/website/www/All_features_part34.arrow")

gene_info_part35 <- gene_info_2[,4801:6000]
gene_info_part35$cell <- rownames(gene_info_part35)
write_feather(gene_info_part35, sink = "D:/Github/website/www/All_features_part35.arrow")

gene_info_part36 <- gene_info_2[,6001:7200]
gene_info_part36$cell <- rownames(gene_info_part36)
write_feather(gene_info_part36, sink = "D:/Github/website/www/All_features_part36.arrow")

gene_info_part37 <- gene_info_2[,7201:8400]
gene_info_part37$cell <- rownames(gene_info_part37)
write_feather(gene_info_part37, sink = "D:/Github/website/www/All_features_part37.arrow")

gene_info_part38 <- gene_info_2[,8401:9600]
gene_info_part38$cell <- rownames(gene_info_part38)
write_feather(gene_info_part38, sink = "D:/Github/website/www/All_features_part38.arrow")

gene_info_part39 <- gene_info_2[,9601:10800]
gene_info_part39$cell <- rownames(gene_info_part39)
write_feather(gene_info_part39, sink = "D:/Github/website/www/All_features_part39.arrow")

gene_info_part40 <- gene_info_2[,10801:12000]
gene_info_part40$cell <- rownames(gene_info_part40)
write_feather(gene_info_part40, sink = "D:/Github/website/www/All_features_part40.arrow")

gene_info_part41 <- gene_info_2[,12001:13200]
gene_info_part41$cell <- rownames(gene_info_part41)
write_feather(gene_info_part41, sink = "D:/Github/website/www/All_features_part41.arrow")

gene_info_part42 <- gene_info_2[,13201:14400]
gene_info_part42$cell <- rownames(gene_info_part42)
write_feather(gene_info_part42, sink = "D:/Github/website/www/All_features_part42.arrow")

gene_info_part43 <- gene_info_2[,14401:15600]
gene_info_part43$cell <- rownames(gene_info_part43)
write_feather(gene_info_part43, sink = "D:/Github/website/www/All_features_part43.arrow")

gene_info_part44 <- gene_info_2[,15601:16800]
gene_info_part44$cell <- rownames(gene_info_part44)
write_feather(gene_info_part44, sink = "D:/Github/website/www/All_features_part44.arrow")

gene_info_part45 <- gene_info_2[,16801:18000]
gene_info_part45$cell <- rownames(gene_info_part45)
write_feather(gene_info_part45, sink = "D:/Github/website/www/All_features_part45.arrow")

gene_info_part46 <- gene_info_2[,18001:19200]
gene_info_part46$cell <- rownames(gene_info_part46)
write_feather(gene_info_part46, sink = "D:/Github/website/www/All_features_part46.arrow")

gene_info_part47 <- gene_info_2[,19201:20400]
gene_info_part47$cell <- rownames(gene_info_part47)
write_feather(gene_info_part47, sink = "D:/Github/website/www/All_features_part47.arrow")

gene_info_part48 <- gene_info_2[,20401:21600]
gene_info_part48$cell <- rownames(gene_info_part48)
write_feather(gene_info_part48, sink = "D:/Github/website/www/All_features_part48.arrow")

gene_info_part49 <- gene_info_2[,21601:22800]
gene_info_part49$cell <- rownames(gene_info_part49)
write_feather(gene_info_part49, sink = "D:/Github/website/www/All_features_part49.arrow")

gene_info_part50 <- gene_info_2[,22801:24000]
gene_info_part50$cell <- rownames(gene_info_part50)
write_feather(gene_info_part50, sink = "D:/Github/website/www/All_features_part50.arrow")

gene_info_part51 <- gene_info_2[,24001:25200]
gene_info_part51$cell <- rownames(gene_info_part51)
write_feather(gene_info_part51, sink = "D:/Github/website/www/All_features_part51.arrow")

gene_info_part52 <- gene_info_2[,25201:26400]
gene_info_part52$cell <- rownames(gene_info_part52)
write_feather(gene_info_part52, sink = "D:/Github/website/www/All_features_part52.arrow")

gene_info_part53 <- gene_info_2[,26401:27600]
gene_info_part53$cell <- rownames(gene_info_part53)
write_feather(gene_info_part53, sink = "D:/Github/website/www/All_features_part53.arrow")

gene_info_part54 <- gene_info_2[,27601:28800]
gene_info_part54$cell <- rownames(gene_info_part54)
write_feather(gene_info_part54, sink = "D:/Github/website/www/All_features_part54.arrow")

gene_info_part55 <- gene_info_2[,28801:30000]
gene_info_part55$cell <- rownames(gene_info_part55)
write_feather(gene_info_part55, sink = "D:/Github/website/www/All_features_part55.arrow")

gene_info_part56 <- gene_info_2[,30001:32400]
gene_info_part56$cell <- rownames(gene_info_part56)
write_feather(gene_info_part56, sink = "D:/Github/website/www/All_features_part56.arrow")

gene_info_part57 <- gene_info_2[,32401:34000]
gene_info_part57$cell <- rownames(gene_info_part57)
write_feather(gene_info_part57, sink = "D:/Github/website/www/All_features_part57.arrow")

gene_info_part58 <- gene_info_2[,34001:35800]
gene_info_part58$cell <- rownames(gene_info_part58)
write_feather(gene_info_part58, sink = "D:/Github/website/www/All_features_part58.arrow")

gene_info_part59 <- gene_info_2[,35801:37600]
gene_info_part59$cell <- rownames(gene_info_part59)
write_feather(gene_info_part59, sink = "D:/Github/website/www/All_features_part59.arrow")

gene_info_part60 <- gene_info_2[,37601:39070]
gene_info_part60$cell <- rownames(gene_info_part60)
write_feather(gene_info_part60, sink = "D:/Github/website/www/All_features_part60.arrow")





### Epidermis
### The arrow files would exceed 50MB so I split the data into several files
features <- rownames(epidermis)
gene_info <- FetchData(epidermis, vars = c(c("umap_1", "umap_2"), "ident", features))

gene_info_part1 <- gene_info[,1:2900]
gene_info_part1$cell <- rownames(gene_info_part1)
write_feather(gene_info_part1, sink = "D:/Github/website/www/Epi_part1.arrow")

gene_info_part2 <- gene_info[,2901:6600]
gene_info_part2$cell <- rownames(gene_info_part2)
write_feather(gene_info_part2, sink = "D:/Github/website/www/Epi_part2.arrow")

gene_info_part3 <- gene_info[,6601:10000]
gene_info_part3$cell <- rownames(gene_info_part3)
write_feather(gene_info_part3, sink = "D:/Github/website/www/Epi_part3.arrow")

gene_info_part4 <- gene_info[,10001:13400]
gene_info_part4$cell <- rownames(gene_info_part4)
write_feather(gene_info_part4, sink = "D:/Github/website/www/Epi_part4.arrow")

gene_info_part5 <- gene_info[,13401:17000]
gene_info_part5$cell <- rownames(gene_info_part5)
write_feather(gene_info_part5, sink = "D:/Github/website/www/Epi_part5.arrow")

gene_info_part6 <- gene_info[,17001:20200]
gene_info_part6$cell <- rownames(gene_info_part6)
write_feather(gene_info_part6, sink = "D:/Github/website/www/Epi_part6.arrow")

gene_info_part7 <- gene_info[,20201:23900]
gene_info_part7$cell <- rownames(gene_info_part7)
write_feather(gene_info_part7, sink = "D:/Github/website/www/Epi_part7.arrow")

gene_info_part8 <- gene_info[,23901:27000]
gene_info_part8$cell <- rownames(gene_info_part8)
write_feather(gene_info_part8, sink = "D:/Github/website/www/Epi_part8.arrow")

gene_info_part9 <- gene_info[,27001:30500]
gene_info_part9$cell <- rownames(gene_info_part9)
write_feather(gene_info_part9, sink = "D:/Github/website/www/Epi_part9.arrow")

gene_info_part10 <- gene_info[,30501:34800]
gene_info_part10$cell <- rownames(gene_info_part10)
write_feather(gene_info_part10, sink = "D:/Github/website/www/Epi_part10.arrow")

gene_info_part11 <- gene_info[,34801:39070]
gene_info_part11$cell <- rownames(gene_info_part11)
write_feather(gene_info_part11, sink = "D:/Github/website/www/Epi_part11.arrow")






### Stomatal cell files
features <- rownames(stomatal_files)
gene_info <- FetchData(stomatal_files, vars = c(c("umap_1", "umap_2"), "ident", features))

gene_info_part1 <- gene_info[,1:7000]
gene_info_part1$cell <- rownames(gene_info_part1)
write_feather(gene_info_part1, sink = "D:/Github/website/www/Stom_part1.arrow")

gene_info_part2 <- gene_info[,7001:14000]
gene_info_part2$cell <- rownames(gene_info_part2)
write_feather(gene_info_part2, sink = "D:/Github/website/www/Stom_part2.arrow")

gene_info_part3 <- gene_info[,14001:21000]
gene_info_part3$cell <- rownames(gene_info_part3)
write_feather(gene_info_part3, sink = "D:/Github/website/www/Stom_part3.arrow")

gene_info_part4 <- gene_info[,21001:28000]
gene_info_part4$cell <- rownames(gene_info_part4)
write_feather(gene_info_part4, sink = "D:/Github/website/www/Stom_part4.arrow")

gene_info_part5 <- gene_info[,28001:34000]
gene_info_part5$cell <- rownames(gene_info_part5)
write_feather(gene_info_part5, sink = "D:/Github/website/www/Stom_part5.arrow")

gene_info_part6 <- gene_info[,34001:39071]
gene_info_part6$cell <- rownames(gene_info_part6)
write_feather(gene_info_part6, sink = "D:/Github/website/www/Stom_part6.arrow")

























# UMAP plots --------------------------------------------------------
library(tidyverse)
library(Seurat)
library(MetBrewer)

### get coordinates for full dataset
label_info_all <- FetchData(whole_data, vars = c(c("umap_1", "umap_2"), "annotated"))
saveRDS(label_info_all, "E:/Scseq_brachy_dev/Files for webpage/All_labelled_STACAS.rds")



### get coordinates for epidermal dataset
label_info_epi <- FetchData(epidermis, vars = c(c("umap_1", "umap_2"), "celltypes"))
saveRDS(label_info_epi, "E:/Scseq_brachy_dev/Files for webpage/Epidermis_labelled_STACAS.rds")



### get coordinates for stomatal dataset
label_info_stom <- FetchData(stomatal_files, vars = c(c("umap_1", "umap_2"), "stages"))
saveRDS(label_info_stom, "E:/Scseq_brachy_dev/Files for webpage/Stomata_labelled_STACAS.rds")
