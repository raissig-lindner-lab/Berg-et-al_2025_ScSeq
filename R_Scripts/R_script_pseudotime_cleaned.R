#### Pseudo-time analysis of ScSeq-data with Monocle

### Script based on http://cole-trapnell-lab.github.io/monocle-release/docs/#the-celldataset-class


### load original dataset
wt_sid_nocca <- readRDS("E:/Scseq_brachy_dev/R_data_objects/wt_sid_nocca_all_genome_v1_2.rds")
epidermis <- readRDS("E:/Scseq_brachy_dev/R_data_objects/epidermis_all_shorter_genome_v1_2.rds")
stomata <- readRDS("E:/Scseq_brachy_dev/R_data_objects/stomatal_lineage_2024.rds")
stomatal_files <- readRDS("E:/Scseq_brachy_dev/R_data_objects/stomatal_cell_files.rds")
gc_lineage <- readRDS("E:/Scseq_brachy_dev/R_data_objects/gc_lineage.rds")
hc_lineage <- readRDS("E:/Scseq_brachy_dev/R_data_objects/hc_lineage.rds")


### load libraries
library(tidyverse)
library(Seurat)
library(MetBrewer)
library(viridis)

############## MONOCLE -------------------------------------------------------------------------------------------

#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("satijalab/seurat-wrappers")
# if seuratwrappers can't install properly because "lazy loading failed for package "SeuratWrappers"", check if there is any dependencies missing

library(monocle3)
library(SeuratWrappers)

## based on https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.md and https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed 

### reupload a previous file
cds <- readRDS("E:/Scseq_brachy_dev/R_data_objects/gc_lineage_monocle.rds")
cds <- readRDS("E:/Scseq_brachy_dev/R_data_objects/hc_lineage_monocle.rds")


### or convert Seurat object into a CellDataSet object
data <- JoinLayers(wt_sid_nocca)
data <- JoinLayers(stomata)
data <- JoinLayers(stomatal_files)
data <- JoinLayers(gc_lineage)
data <- JoinLayers(hc_lineage)

cds <- as.cell_data_set(data)

cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(data[["RNA"]])
cds <- preprocess_cds(cds, num_dim = 50)

cds <- cluster_cells(cds, resolution=0.001) # GC lineage
cds <- cluster_cells(cds, resolution = 0.0001) # HC lineage

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1 + p2


cds <- learn_graph(cds, use_partition = F, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- estimate_size_factors(cds)


### colour cells by pseudotime
# check with featureplots which are the earliest cells, e.g. the root cluster
plot_cells(cds, genes = "BdiBd21-3.3G0074300", show_trajectory_graph = F, group_label_size = 5)


## define root cluster
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2])) # gc_lineage
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1])) # hc_lineage

  
## plot
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black")

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           cell_size = 0.1,
           cell_stroke = 1.5,
           trajectory_graph_color = NA)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = NA)

plot_cells(cds,
           color_cells_by = "stages",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = NA)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "white") + DarkTheme()

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = NA) + DarkTheme()


### export this information to the seurat object and visualise
gc_lineage$pseudotime <- pseudotime(cds)
FeaturePlot(gc_lineage, "pseudotime")
 
saveRDS(cds, "E:/Scseq_brachy_dev/R_data_objects/hc_lineage_monocle.rds")
saveRDS(cds, "E:/Scseq_brachy_dev/R_data_objects/gc_lineage_monocle.rds")





###### Expression patterns of selected genes over time (Monocle Pseudotime) -------------------------------------------------
### exctract info on genes of interest
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("BdiBd21-3.1G0523400", # SPCH1
                                                 "BdiBd21-3.1G0240400", # MUTE
                                                 "BdiBd21-3.2G0300000", # FAMA
                                                 "BdiBd21-3.1G0206300", # SCAP1
                                                 "BdiBd21-3.4G0234500" # MYB60
                          )))


genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("BdiBd21-3.2G0467800", # BdPOX
                                                 "BdiBd21-3.1G0421000", # BdSPL10
                                                 "BdiBd21-3.2G0477300", # ZmWOX3B
                                                 "BdiBd21-3.3G0074300", # BdSPL5
                                                 "BdiBd21-3.3G0655900", # 2cRR
                                                 "BdiBd21-3.3G0292700", # SLP1-like3
                                                 "BdiBd21-3.1G0921800" # SLP1-like5
                                                 ))) 



cds_subset <- cds[genes,]







### visualize gene expression over pseudotime gradient
plot_genes_in_pseudotime(cds_subset, color_cells_by = "stages",
                         panel_order = c("BdiBd21-3.1G0523400", # SPCH1
                                         "BdiBd21-3.1G0240400", # MUTE
                                         "BdiBd21-3.2G0300000", # FAMA
                                         "BdiBd21-3.1G0206300", # SCAP1
                                         "BdiBd21-3.4G0234500" # MYB60
                         ),
                         min_expr = 0.005,
                         cell_size = 1)


plot_genes_in_pseudotime(cds_subset, color_cells_by = "pseudotime",
                         panel_order = c("BdiBd21-3.3G0292700", # SLP1-like3
                                         "BdiBd21-3.2G0467800", # BdPOX
                                         "BdiBd21-3.1G0921800", # SLP1-like5
                                         "BdiBd21-3.3G0655900", # 2cRR
                                         "BdiBd21-3.1G0421000", # BdSPL10
                                         "BdiBd21-3.2G0477300", # ZmWOX3B
                                         "BdiBd21-3.3G0074300" # BdSPL5
                                         ),
                         min_expr = 0.005,
                         cell_size = 1)







### as a function to plot gene expression over pseudotime
gene_exp_over_pseudotime <- function(genes, colour_by = "seurat_clusters", 
                                     colours = met.brewer("Hiroshige"), min_exp = 0.005) {
  goi <- row.names(subset(fData(cds), gene_short_name %in% genes))
  
  cds_subset <- cds[goi,]
  
  plot_genes_in_pseudotime(cds_subset, 
                           color_cells_by = colour_by, 
                           panel_order = genes, 
                           min_expr = min_exp,
                           cell_size = 1) +
    scale_colour_manual(values = colours)
}




## SPCH2, MUTE, FAMA, SCAP1, MYB60 (SPCH2 at bottom) with the GC lineage subset
gene_exp_over_pseudotime(c("BdiBd21-3.4G0234500", "BdiBd21-3.1G0206300", "BdiBd21-3.2G0300000", "BdiBd21-3.1G0240400", "BdiBd21-3.1G0523400"),
                         colour_by = "stages",
                         colours = c("#78C7B8", "#1D4573", "#0A2E57", "#538EB9", "#356493"))