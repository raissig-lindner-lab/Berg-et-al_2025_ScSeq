### Code to plot featureplots with 3 colours
library(tidyverse)
library(arrow)
library(Seurat)

setwd("D:/Github/website/")

gene_info_gc <- FetchData(gc_lineage, vars = c(c("umap_1", "umap_2"), features = c("BdiBd21-3.3G0131200", # SPCH2
                                                                                   "BdiBd21-3.1G0240400", # MUTE
                                                                                   "BdiBd21-3.2G0300000", # FAMA
                                                                                   "BdiBd21-3.1G0206300", # SCAP1
                                                                                   "BdiBd21-3.4G0234500", # MYB60
                                                                                   "BdiBd21-3.2G0561000", # TMM
                                                                                   "BdiBd21-3.1G1009000", # SDD1-like
                                                                                   "BdiBd21-3.2G0614100", # LRR kinase
                                                                                   "BdiBd21-3.2G0049500" # GELP1
                                                                                   )))

gene_info_hc <- FetchData(hc_lineage, vars = c(c("umap_1", "umap_2"), features = c("BdiBd21-3.3G0292700", # SLP1-like3
                                                                                   "BdiBd21-3.2G0467800", # POX
                                                                                   "BdiBd21-3.1G0921800", # SLP1-like5
                                                                                   "BdiBd21-3.3G0655900", # 2cRR
                                                                                   "BdiBd21-3.1G0421000", # SPL10
                                                                                   "BdiBd21-3.2G0477300", # WOX3B
                                                                                   "BdiBd21-3.3G0074300" # SPL5
                                                                                   )))


### read files with data for labelled UMAP plots
label_info_all <- readRDS("D:/Github/website/www/All_labelled.rds")
label_info_epi <- readRDS("D:/Github/website/www/Epidermis_labelled.rds")
label_info_stom <- readRDS("D:/Github/website/www/Stomata_labelled.rds")

### Location file required for feature plots
location_file <- readRDS("D:/Github/website/www/Location_file_arrowfiles.rds") # this is the same file as for the web tool

featureplot_website <- function(dataset, gene=NULL) {
  if(gene == "x") {
    ggplot()+geom_point()+theme_classic()
  }
  
  else {
    ### check which files have the data for the selected gene and load the respective files
    file_name_all1 <- location_file[gene, "all_cells1"]
    file_name_all2 <- location_file[gene, "all_cells2"]
    file_name_epi <- location_file[gene, "epidermis"]
    file_name_stom <- location_file[gene, "stomatal_cell_files"]
    
    dt1 <- as.data.frame(read_feather(file_name_all1, col_select = c("cell", gene)))
    rownames(dt1) <- dt1$cell
    
    dt2 <- as.data.frame(read_feather(file_name_all2, col_select = c("cell", gene)))
    rownames(dt2) <- dt2$cell
    
    gene_info_all <- rbind(dt1, dt2)
    gene_info_all <- merge(gene_info_all, label_info_all[,1:2], by = c("row.names"))
    if(file_name_all1 == "www/All_features_part1.rds") {
      gene_info_all$umap_1 <- gene_info_all$umap_1.x
      gene_info_all$umap_2 <- gene_info_all$umap_2.x
    }
    
    gene_info_epi <- as.data.frame(read_feather(file_name_epi, col_select = c("cell", gene)))
    rownames(gene_info_epi) <- gene_info_epi$cell
    gene_info_epi <- merge(gene_info_epi, label_info_epi[,1:2], by = c("row.names"))
    if(file_name_epi == "www/Epi_part1.rds") {
      gene_info_epi$umap_1 <- gene_info_epi$umap_1.x
      gene_info_epi$umap_2 <- gene_info_epi$umap_2.x
    }
    
    gene_info_stom <- as.data.frame(read_feather(file_name_stom, col_select = c("cell", gene)))
    rownames(gene_info_stom) <- gene_info_stom$cell
    gene_info_stom <- merge(gene_info_stom, label_info_stom[,1:2], by = c("row.names"))
    if(file_name_stom == "www/Stom_part1.rds") {
      gene_info_stom$umap_1 <- gene_info_stom$umap_1.x
      gene_info_stom$umap_2 <- gene_info_stom$umap_2.x
    }
    
    
    if(dataset == "All cells") {
      feat_data <- gene_info_all[,c("umap_1", "umap_2", gene)]
      
      feat_data <- feat_data[order(feat_data[,3]),]
      
      ggplot(feat_data)+
        geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]), size = 0.5)+
        scale_colour_gradient2(low = "#FFFAA0", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
        theme_classic() +
        labs(colour = "Expression", title = gene)
    }
    
    else {
      if(dataset == "Epidermis") {
        feat_data <- gene_info_epi[,c("umap_1", "umap_2", gene)]
        
        feat_data <- feat_data[order(feat_data[,3]),]
        
        ggplot(feat_data)+
          geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]), size = 0.5)+
          scale_colour_gradient2(low = "#FFFAA0", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
          theme_classic()+
          labs(colour = "Expression", title = gene)
      }
      
      else {
        if(dataset == "Stomatal cell files") {
          feat_data <- gene_info_stom[,c("umap_1", "umap_2", gene)]
          
          feat_data <- feat_data[order(feat_data[,3]),]
          
          ggplot(feat_data)+
            geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]), size = 0.7)+
            scale_colour_gradient2(low = "#FFFAA0", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
            theme_classic() +
            labs(colour = "Expression", title = gene)
        }
        
        else {
          if(dataset == "GC lineage") {
            feat_data <- gene_info_gc[,c("umap_1", "umap_2", gene)]
            
            feat_data <- feat_data[order(feat_data[,3]),]
            
            ggplot(feat_data)+
              geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]), size = 0.7)+
              scale_colour_gradient2(low = "#FFFAA0", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
              theme_classic() +
              labs(colour = "Expression", title = gene)
          }
          
          else {
            if(dataset == "HC lineage") {
              feat_data <- gene_info_hc[,c("umap_1", "umap_2", gene)]
              
              feat_data <- feat_data[order(feat_data[,3]),]
              
              ggplot(feat_data)+
                geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]), size = 0.7)+
                scale_colour_gradient2(low = "#FFFAA0", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
                theme_classic() +
                labs(colour = "Expression", title = gene)
            }
          }
        }
      }
    }
  }
}


### Featureplots for paper (save as 3.5x3.5 PDF)
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.5G0316500") # WOX4
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.5G0168500") # FCP1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0402200") # CLV1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0427700") # FEA2
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0010000") # FEA3
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0571300") # FEA4
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0588300") # PIN1a
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0135700") # KN1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0773000") # KNAT1-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0348000") # TMO6-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G1023300") # VND-like1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.5G0221500") # VND-like2
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0501500") # XCP1-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.3G0071200") # APL
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0203400") # SLAH2
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.3G0204900") # PDF2-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.3G0192600") # HDG2-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0082000") # PDF1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0587500") # WOX9C-like1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0749400") # STOMAGEN-1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.4G0439300") # LHCA6
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0657800") # GRAS32 
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0164900") # Pavement cell GELP 
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.4G0052400") # AS1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.2G0089600") # TCP4-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0741600") # TSO1
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.5G0296400") # ML1-like
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.3G0783600") # PIN1b
featureplot_website(dataset = "All cells", gene = "BdiBd21-3.1G0942300") # CRC


featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0074300") # BdSPL5
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.2G0477300") # BdWOX3B
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0523400") # BdSPCH1
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0131200") # BdSPCH2
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0240400") # BdMUTE
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.2G0300000") # BdFAMA
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.4G0234500") # BdMYB60-like
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0206300") # BdSCAP1
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0662500") # BdERL1
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.4G0254800") # BdICE1
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.2G0762000") # BdSCRM2
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0655900") # 2cRR
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0657800") # GRAS32
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0292100") # SLP1-like 1
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0292000") # SLP1-like 2
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0292700") # SLP1-like 3
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.5G0061500") # SLP1-like 4
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0921800") # SLP1-like 5
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0921700") # SLP1-like 6
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0291900") # SLP1-like 7
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0292200") # SLP1-like 8
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0855400") # SLP1-like 9
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.2G0255500") # PME53-like
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.2G0561000") # TMM
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.1G0164900") # Pavement cell GELP
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.3G0181200") # L-gulonolactone oxidase
featureplot_website(dataset = "Epidermis", gene = "BdiBd21-3.5G0356800") # peroxidase


featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0131200") # SPCH2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0240400") # MUTE
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0300000") # FAMA
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0326500") # CST1 
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.4G0612300") # DUF567
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0255500") # PME53-like
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0614100") # LRR kinase
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0049500") # GELP1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0523400") # SPCH1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.4G0234500") # MYB60-like
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.4G0254800") # ICE1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0762000") # SCRM2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.5G0153600") # EPF2-1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.5G0306100") # EPF2-2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0609800") # ERECTA
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0662500") # ERL1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.5G0153700") # SERK1/2-like1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0620800") # SERK1/2-like2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0209300") # SERK1/2-like3
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.2G0561000") # TMM
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G1009000") # SDD1-like
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.5G0238000") # YDA1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0680900") # YDA2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0885800") # MPK3-like
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0650300") # MPK6-like
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0616400") # MKK4/5-like1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0709500") # MKK4/5-like2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0526300") # PAN1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0783500") # PAN2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0715200") # POLAR
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.4G0224600") # MLKS2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.3G0022600") # ROP9
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G1008100") # LPL2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0734900") # TAN1
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0233100") # WPRB2
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0164900") # Pavement cell GELP
featureplot_website(dataset = "Stomatal cell files", gene = "BdiBd21-3.1G0206300") # SCAP1

featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.3G0131200") # SPCH2
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.1G0240400") # MUTE
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.2G0300000") # FAMA
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.1G0206300") # SCAP1
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.4G0234500") # MYB60-like
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.2G0561000") # TMM
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.1G1009000") # SDD1-like
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.2G0614100") # LRR kinase
featureplot_website(dataset = "GC lineage", gene = "BdiBd21-3.2G0049500") # GELP1

featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.3G0292700") # SLP1-like3 
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.2G0467800") # POX
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.1G0921800") # SLP1-like5
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.3G0655900") # 2cRR
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.1G0421000") # SPL10
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.2G0477300") # WOX3B
featureplot_website(dataset = "HC lineage", gene = "BdiBd21-3.3G0074300") # SPL5














### To plot two genes against each other in a Feature Plot (as in Fig. 6C) -----------------------------------------------
doi <- stomatal_files
doi <- AddModuleScore(doi, features = "BdiBd21-3.1G0240400", nbin = 2, name = "gene_expression") # MUTE
doi <- AddModuleScore(doi, features = "BdiBd21-3.2G0300000", nbin = 2, name = "contra_expression") # FAMA

### categorise cells into cells that have and don't have the genes of interest expressed
doi$goi_score <- ifelse(test = doi$gene_expression1 > 0 & doi$contra_expression1 <= 0, yes = "gene",
                        no = ifelse(test = doi$gene_expression1 <= 0 & doi$contra_expression1 > 0, yes = "contra_gene",
                                    no = ifelse(test = doi$gene_expression1 > 0 & doi$contra_expression1 > 0, yes = "both",
                                                no = "neither")))

DimPlot(doi, reduction = "umap", group.by = "goi_score", 
        order = list("gene", "contra_gene", "both", "Neither"), 
        cols = c("#ECECEC", "#7e2e7b", "#ff00ff", "#71c837"))
