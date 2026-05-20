###### Gene regulatory analysis (GRN) using Mini-EX (https://github.com/VIB-PSB/MINI-EX)
library(Seurat)
library(tidyverse)
library(MetBrewer)

stomatal_files <- readRDS("E:/Scseq_brachy_dev/R_data_objects/stomatal_cell_files_withoutprim_STACAS.rds")

DimPlot(stomatal_files)
Idents(stomatal_files) <- "seurat_clusters"
gc_sc_lineage <- subset(stomatal_files, idents = c(3,5,6,7,8,9,10,11,12,13,14,16))
DimPlot(gc_sc_lineage)

### For paper: mark guard cell and subsidiary cell lineage on stomatal dataset
gc_sc_cells <- WhichCells(gc_sc_lineage)

DimPlot(stomatal_files, cells.highlight = c(gc_sc_cells), sizes.highlight = 0.3)


################# !!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!
### When uploading files to the HPC or any Linux system and before running the MINI-EX.sh script, use dos2unix on all the input files!


###### Prepare input files for MINI-EX -------------------------------------------------------------

expression.matrix <- as.data.frame(as.matrix(GetAssayData(object = gc_sc_lineage, assay = "RNA", layer = "counts")))
write.table(expression.matrix, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/stom_matrix_stages_sub.tsv", sep='\t', quote = FALSE)



### Find all differentially expressed markers (we only need upregulated ones)
## The markersOut points to the output obtained by Seurat FindAllMarkers using the command below:
# for subsetted GC and SC lineages
Idents(gc_sc_lineage) <- "stages" # by stage label
gc_sc_lineage <- RenameIdents(gc_sc_lineage, "Stage 0-2" = "0", 
                               "Early GMCs" = "1",
                               "Dividing GMCs" = "2",
                               "Early GCs" = "3",
                               "Late GCs" = "4",
                               "SMCs" = "5",
                               "SCs" = "6",
                               "Inter-specialized cells" = "7")
stom_markers_all <- FindAllMarkers(gc_sc_lineage, only.pos=TRUE)
write.table(stom_markers_all, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/stom_allMarkers_stages_sub.tsv", sep = "\t", quote=FALSE)



### The cellsToClusters points to a tab-separated file containing the cluster annotation for each cell in the expression matrix. It can be obtained using the Seurat command FetchData as shown below:
cells2clusters <-FetchData(gc_sc_lineage, vars = 'ident')

write.table(cells2clusters,"E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/stom_cells2clusters_stages_sub.tsv", sep="\t",quote=FALSE,col.names = FALSE)

# for seurat clustered analysis of the GC SC subset
clusterlabels <- data.frame(cluster=c("0", "1", 
                                      "2", "3", 
                                      "4", "5",             
                                      "6", "7"), 
                            label=c("Stage0to2", "EarlyGMCs", 
                                    "DividingGMCs", "EarlyGCs",
                                    "LateGCs", "SMCs", 
                                    "SCs", "Intercells"))

### The clustersToIdentities points to a tab-separated file containing the cell type annotation for each cluster. 
## Optionally, this file can contain a third column which specifies a cluster index. 
## This index is used to indicate the position of a cluster along a known developmental trajectory (see miniexExample_identities_with_idx.tsv in the INPUTS folder), 
## and translates to the column index in the regulator heatmaps (example). 
## If this column is not specified, an automatic index is created by sorting the cluster identities alphabetically.
clusters2Identities <- matrix(clusterlabels$label) 
rownames(clusters2Identities) <- 0:7 # for stage clustered input of the GC SC subset

## save for stages clustered input of the subset
write.table(clusters2Identities,"E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/stom_identities_stages_sub.tsv", sep="\t",quote=FALSE,col.names = FALSE)


### The TF_list is a list of TFs which is used in the GRNBoost2 run
tfs_all <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Bdi_TF_list.txt")

synonyms <- read.delim("E:/Scseq_brachy_dev/markers/BdistachyonBd21_3_460_v1.1.synonym.txt", header = F)
colnames(synonyms) <- c("BdiBd21_3", "Gene_ID")
synonyms$BdiBd21_3 <-substr(synonyms$BdiBd21_3,1,19) #remove transcript denominator (.1)
synonyms <- unique(synonyms)

tfs_all_2 <- merge(synonyms, tfs_all, by = "Gene_ID")

tfs_all_2 <- tfs_all_2[,c("BdiBd21_3", "Family")]
tfs_all_2 <- unique(tfs_all_2)

tf_all_list <- data.frame(tfs_all_2$BdiBd21_3)
write.table(tf_all_list, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/TF_list_brachy_all.tsv", sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)



### list with transcription factor motifs
tf2fam2motif <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Bdi_TF_binding_motifs_information.txt")

## merge with synonyms file to convert Bradi accession to Bd21-3
colnames(synonyms) <- c("BdiBd21_3", "Gene_id")

tf2fam2motif_2 <- merge(synonyms, tf2fam2motif, by = "Gene_id")

## this file now has 5 motifs less, this gives an error in Mini-EX because it can't find them
## which motifs?
df1 <- tf2fam2motif %>% select(Gene_id, Matrix_id)
df2 <- tf2fam2motif_2 %>% select(Gene_id, BdiBd21_3)

df3 <- merge(df1, df2, by = "Gene_id", all=T)
print(df3[c(105, 118, 146, 199, 224),])
#           Gene_id Matrix_id BdiBd21_3
#  105 Bradi2g36810   MP00434      <NA>
#  118 Bradi2g48057   MP00382      <NA>
#  146 Bradi3g01483   MP00229      <NA>
#  199 Bradi3g59320   MP00224      <NA>
#  224 Bradi4g32967   MP00219      <NA>

## these genes apparently do not exist in the Bd21-3 genome
## I will filter out these genes from the final FIMO output file used for MINI-EX ("Target_genes_all.out")
## done externally on the server in R by using the following lines of code:
#library(tidyverse)
#data <- read.delim("Target_genes_all.out", header=F)
#data2 <- data %>% filter(!V1 == "MP00434")
#data2 <- data2 %>% filter(!V1 == "MP00382")
#data2 <- data2 %>% filter(!V1 == "MP00229")
#data2 <- data2 %>% filter(!V1 == "MP00224")
#data2 <- data2 %>% filter(!V1 == "MP00219")
#write.table(data2, "Target_genes_all.out", sep = "\t", quote=F, col.names=F, row.names=F)



## add motif names for MUTE, FAMA, SPCH1/2, ICE1/SCRM2
extras <- data.frame("Gene_id" = c("Bradi1g18400", "Bradi2g22810", "Bradi1g38650", "Bradi3g09670", "Bradi4g17460", "Bradi2g59497"), 
                     "BdiBd21_3" = c("BdiBd21-3.1G0240400", "BdiBd21-3.2G0300000", "BdiBd21-3.1G0523400", "BdiBd21-3.3G0131200", "BdiBd21-3.4G0254800", "BdiBd21-3.2G0762000"),
                     "Family" = c("bHLH", "bHLH", "bHLH", "bHLH", "bHLH", "bHLH"),
                     "Matrix_id" = c("MUTE", "FAMA", "SPCH1", "SPCH2", "ICE1", "SCRM2"),
                     "Species" = "Brachypodium distachyon",
                     "Method" = "-",
                     "Datasource" = "-",
                     "Datasource_ID" = "-")
tf2fam2motif_2 <- rbind(tf2fam2motif_2, extras)

tf2fam2motif_2 <- tf2fam2motif_2[,c("BdiBd21_3", "Family", "Matrix_id")]



### add column specifying whether a TF has motif information available
### and reorder to have this order: TF, TF family, motif info available?, motif name
tf2fam2motif_2$Info <- "Y"
tf2fam2info2motif <- merge(tfs_all_2, tf2fam2motif_2, by = "BdiBd21_3", all = T)
tf2fam2info2motif <- tf2fam2info2motif %>% select(-Family.y)

colnames(tf2fam2info2motif)
tf2fam2info2motif <- tf2fam2info2motif[,c(1,2,4,3)]

for(row in 1:nrow(tf2fam2info2motif)) {
  if(is.na(tf2fam2info2motif$Info[row])) {
    tf2fam2info2motif$Info[row] <- "N"
  }
}

write.table(tf2fam2info2motif, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/new INPUTs/TF2fam2motif.tsv", quote = F, col.names = F, row.names = F, sep = "\t")



### ### transcription factor binding motifs
motif_info <- universalmotif::read_meme("E:/Scseq_brachy_dev/GRN_analysis/Bdi_TF_binding_motifs.meme")

## add bHLH motif for MUTE, FAMA, SPCH1/2, ICE1/SCRM2
mute <- universalmotif::create_motif("CANNTG", name = "Bradi1g18400", altname = "MUTE",
                                     alphabet = "DNA", type = "PPM", strand = "+-")
universalmotif::view_motifs(mute)

fama <- universalmotif::create_motif("CANNTG", name = "Bradi2g22810", altname = "FAMA",
                                     alphabet = "DNA", type = "PPM", strand = "+-")

spch1 <- universalmotif::create_motif("CANNTG", name = "Bradi1g38650", altname = "SPCH1",
                                      alphabet = "DNA", type = "PPM", strand = "+-")

spch2 <- universalmotif::create_motif("CANNTG", name = "Bradi3g09670", altname = "SPCH2",
                                      alphabet = "DNA", type = "PPM", strand = "+-")

ice1 <- universalmotif::create_motif("CANNTG", name = "Bradi4g17460", altname = "ICE1",
                                     alphabet = "DNA", type = "PPM", strand = "+-")

scrm2 <- universalmotif::create_motif("CANNTG", name = "Bradi2g59497", altname = "SCRM2",
                                      alphabet = "DNA", type = "PPM", strand = "+-")


motif_info <- append(motif_info, c(mute, fama, spch1, spch2, ice1, scrm2))


## write out file in .meme format and use in FIMO (https://meme-suite.org/meme/tools/fimo) to get TF binding sites (e.g. targets) --> use p-value threshold 0.01. It did not find matches for the stomatal TFs i added when using the default threshold setting (e-04)
universalmotif::write_meme(motif_info, "E:/Scseq_brachy_dev/GRN_analysis/Bdi_TF_binding_motifs_plus_stom_tfs.meme", overwrite = T)


### to split a fasta file into shorter segments (needed for the chromosomes that are too big for FIMO)
## I split chromosome 1 and 3 into three pieces each, chromosome 2 into two pieces and chromosome 4 and 5 were small enough to go as is
fasta <- Biostrings::readDNAStringSet("E:/Scseq_brachy_dev/JGI_files/genome_Bd21-3_v1.2/Shorter fasta files/Bd21_3_537_chr4.fa")
fasta <- (as.character(fasta$Bd1))
fasta <- (as.character(fasta$Bd2))
fasta <- (as.character(fasta$Bd3))
fasta <- (as.character(fasta$Bd4))

fasta2 <- as.list(substring(fasta, seq(1,nchar(fasta),20000000), seq(20000000,nchar(fasta),20000000)))
names(fasta2) <- 1:length(fasta2)

filenames <- names(fasta2)

# make sure to replace the number in the filename with the correct chromosome
for(i in filenames) {
  write.table(fasta2[i], paste("E:/Scseq_brachy_dev/JGI_files/genome_Bd21-3_v1.2/Shorter fasta files/Bd21_3_537_chr4_", i, ".fa"), 
              quote = F, col.names = F, row.names = F)
}

## For these files: Check whether the first line has ">info". 
## If this part is missing, FIMO won't recognize it as a fasta file so you need to manually add the first line



















###### FIMO, finding transcription factor binding sites (https://meme-suite.org/meme/tools/fimo) ----------------------------
### Use the above created .meme file and .fa files in the FIMO web tool to search for transcription factor binding sites
## i set "match p-value" to < 0.01 for chromosome 4 and 5 and <0.1 for the fragments of chromosome 1-3
## the rest is default FIMO settings






###### The next step with the FIMO files before using MINI-EX was done on the server as it takes much less time in a linux system due to possible parallel computing
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)

### for each FIMO file repeat the following (make sure to correct the Chr number in the script before executing)

### Load FIMO file
fimo_chr4 <- as.data.frame(rtracklayer::readGFF("FIMO_Chr4.gff"))
class(fimo_chr4$Alias) <- "character"
fimo_chr4 <- select(fimo_chr4, c(start, end, Alias))


### Load gene annotations
brachy_gff <- rtracklayer::readGFF("BdistachyonBd21_3_537_v1.2.gene_exons.gff3.gz")
bd_chr4 <- as.data.frame(brachy_gff) %>% filter(seqid == "Bd4" & type == "gene") %>% select(c(start, end, Name))

#print(tail(bd_chr1))


### split up fimo file into smaller pieces to speed up the process and now work on OnDemand version of R Studio from UBELIX 
### (otherwise the parallelization will not work because R on Windows only uses one core)
fimo_chr4_list <- split(fimo_chr4, ceiling(seq(nrow(fimo_chr4))/50000))

#test <- list(fimo_chr2_1[1:5000,], fimo_chr2_1[5001:10000,])

final_targets <- data.frame(motif_name = NA, target_gene = NA)
targetted <- data.frame(motif_name = NA, target_gene = NA)

map_to_gene2 <- function(x, bd_chr_object) {
  
  print("Start piece")
  
  for(gene in 1:nrow(bd_chr_object)) {
    print(bd_chr_object[gene, "Name"])
    
    for(line in 1:nrow(x)) {
      
      if(as.numeric(x$start[line]) >= as.numeric(bd_chr_object$start[gene])) {
        
        if(as.numeric(x$end[line]) <= as.numeric(bd_chr_object$end[gene])) {
          print(x$Alias[line])
          
          new_line <- data.frame(motif_name = x$Alias[line], target_gene = bd_chr_object[gene, "Name"])
          targetted <- na.omit(rbind(targetted, new_line))
        }
      }
    }
  }
  
  final_targets <- na.omit(rbind(final_targets, targetted))
  return(final_targets)
}


target_genes_chr4 <- data.frame(motif_name = NA, target_gene = NA)

target_genes_chr4 <- parallel::mclapply(fimo_chr4_list, map_to_gene2, bd_chr_object = bd_chr4, mc.cores = 16)

target_genes_chr4 <- as.data.frame(do.call(rbind, target_genes_chr4))

write.table(target_genes_chr4, "Target_genes_Chr4_allhits.txt", sep = "\t", quote = F, col.names = F, row.names = F)








###### Merge the target genes files to get the final input files for MINI-EX ---------------------------------------------------------------------------
### first merge the files per chromosome for the ones that were split before FIMO
dt1 <- read.delim("Target_genes_Chr1_1_allhits.txt", header = F)
dt2 <- read.delim("Target_genes_Chr1_2_allhits.txt", header = F)
dt3 <- read.delim("Target_genes_Chr1_3_allhits.txt", header = F)

dt <- rbind(dt1, dt2)
dt <- rbind(dt, dt3)

write.table(dt, "Target_genes_Chr1_allhits.txt", sep = "\t", quote = F, col.names = F, row.names = F)

### repeat for chromosomes 2 and 3

### now merge the chromosome files together and filter so that each motif X target gene combination only appears once
dt1 <- read.delim("Target_genes_Chr1_allhits.txt", header = F)
dt2 <- read.delim("Target_genes_Chr2_allhits.txt", header = F)
dt3 <- read.delim("Target_genes_Chr3_allhits.txt", header = F)
dt4 <- read.delim("Target_genes_Chr4_allhits.txt", header = F)
dt5 <- read.delim("Target_genes_Chr5_allhits.txt", header = F)

dt <- rbind(dt1, dt2)
dt <- rbind(dt, dt3)
dt <- rbind(dt, dt4)
dt <- rbind(dt, dt5)

dt <- unique(dt)

write.table(dt, "Target_genes_all.out", sep = "\t", quote = F, col.names = F, row.names = F)

### gzip the file using command line ("gzip Target_genes_all.out")



















###### MINI-EX --------------------------------------------------------------------------------------------------------------------

### Install Mini-EX from the github (https://github.com/VIB-PSB/MINI-EX), make sure you downloaded all the necessary files from there
### this means you need to install Nextflow and Singularity as described
### and you need to download the "bin" folder, as well as "Dockerfile", "miniex.config", "miniex.nf" and "requirements.txt"

### Then you need to configure the "miniex.config" file to refer to your input files (the ones we created above), this can be done in any text editor
### follow the instructions from the MINI-EX github: https://github.com/VIB-PSB/MINI-EX/blob/main/docs/configuration.md

### Make sure that you have all the files in the right position and then run MINI-EX












###### Work with the MINI-EX output files -------------------------------------------------------------------------------------------
### Which transcription factors are expressed?
tf_info <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_TF_info_file.tsv")

### How many are expressed?
table(tf_info$isTF_expressed)
# 24 are not expressed, 1738 are expressed



### In which clusters does a transcription factor act?
library(tidyverse)

tf_tg <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_regulons.tsv", header = F)

one_tf <- tf_tg %>% filter(V1 == "BdiBd21-3.1G0240400")
print(one_tf$V2)



### GSEA without Mini-EX
library(tidyverse)
library(GO.db)
library(AnnotationDbi)

##get GO term info
term2gene <- read_delim("E:/Scseq_brachy_dev/GO_terms_Brachy_UniProtDB_edited.txt") ### this is the list of GO terms for the Brachy genome from UniProt, I used my synonyms list to edit it for "Bd21-3"

bio_process <- data.frame(GO=unique(AnnotationDbi::as.data.frame(GOBPANCESTOR)[,1]))
#mol_function <- data.frame(GO=unique(AnnotationDbi::as.data.frame(GOMFANCESTOR)[,1]))
goterms <- Term(GOTERM)
term2name <- data.frame("GO"=names(goterms), "term"=goterms)
BP_term2name <- merge(bio_process, term2name, by="GO")
#MF_term2name <- merge(mol_function, term2name, by="GO")

go_info <- merge(term2gene, BP_term2name, by = "GO")
#go_info <- merge(term2gene, MF_term2name, by = "GO")


go_terms <- merge(term2gene, term2name, by = "GO")
write.table(go_terms, "E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis.csv")
go_terms <- read.table("E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis.csv")

write.table(go_info, "E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis_BP.csv") # biological process GO_terms

write.table(go_info, "E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis_MF.csv") # molecular function GO terms






#### Dot Plot ranked regulons/TFs (like regmap output of Mini-EX)
library(tidyverse)
regulons <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_rankedRegulons.tsv")
regulons <- as.data.frame(regulons)
regulons <- regulons %>% dplyr::select(-GOterm, -GOdescription, -hasTFrelevantGOterm)

#### add GO terms
go_terms <- read.table("E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis.csv")
go_terms$TF <- go_terms$Bd21_3_id

regulons <- merge(regulons, go_terms[, c(1,3,4)], by = "TF", all.x = T)

#### add gene info
info <- read_delim("E:/Scseq_brachy_dev/brachy_gene_protein_info_Bd21-3_v1.2.csv", delim=";")
info$TF <- info$Bd21_3_id
info <- info %>% dplyr::select(-GO, -Bd21_3_id)

regulons <- merge(regulons, info, by = "TF", all.x = T)

### write out to files
#write_delim(regulons, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_rankedRegulons_moreinfo.txt")

### read saved files
#regulons <- read_delim("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_rankedRegulons_moreinfo.txt")

### select only the columns i need for plotting
regulons <- regulons %>% dplyr::select(TF, cluster, isTF_DE, borda_clusterRank)
regulons <- unique(regulons)


clusterlabels <- data.frame(cluster=c("Stage0to2_Cluster_0", "EarlyGMCs_Cluster_1", 
                                      "DividingGMCs_Cluster_2", "EarlyGCs_Cluster_3", 
                                      "LateGCs_Cluster_4", "SMCs_Cluster_5",             
                                      "SCs_Cluster_6", "Intercells_Cluster_7"), 
                            label=c("Stage 0-2", "Early GMCs",
                                    "Dividing GMCs", "Early GCs", 
                                    "Late GCs", "SMCs",
                                    "SCs", "Inter-specialized cells"))
rownames(clusterlabels) <- clusterlabels$cluster


### filter for top x TFs per cluster
reg_topx <- data.frame(TF=NA, cluster=NA, isTF_DE=NA, borda_clusterRank=NA)
for(clusternr in unique(regulons$cluster)) {
  new_cluster <- regulons %>% filter(cluster==clusternr)
  new_cluster <- slice_max(new_cluster, order_by = dplyr::desc(borda_clusterRank), n = 7) # change n = for different numbers of top X transcription factors
  new_cluster$cluster <- as.character(clusterlabels[clusternr,]$label)
  
  reg_topx <- na.omit(rbind(reg_topx, new_cluster))
}



rownames(reg_topx) <- NULL

reg_topx <- reg_topx[order(factor(reg_topx$cluster, levels = c("Stage 0-2", "Early GMCs",
                                                               "Dividing GMCs", "Early GCs", 
                                                               "Late GCs", "SMCs",
                                                               "SCs", "Inter-specialized cells"))),]

reg_topx %>% mutate(cluster = factor(cluster, levels = c("Stage 0-2", "Early GMCs",
                                                         "Dividing GMCs", "Early GCs", 
                                                         "Late GCs", "SMCs",
                                                         "SCs", "Inter-specialized cells"))) %>% 
  arrange(reg_topx) %>%
  ggplot(aes(x = cluster, y = factor(TF, levels = unique(reg_topx$TF)), size = dplyr::desc(borda_clusterRank), colour = cluster))+
  geom_point()+
  guides(colour = "none",
         size = guide_legend(title = "Cluster-specific Borda rank"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        legend.position = "top")+
  scale_colour_manual(values = c("#78C7B8",
                                 "#538EB9", "#356493", 
                                 "#1D4573", "#0A2E57",
                                 "#8CC483", "#D8D97A", 
                                 "#67AFC2"))+
  labs(x=NULL, y=NULL)


# 8 cell stage/type annotations
"#78C7B8"  # Stage 0-2
"#67AFC2" # Interspecialized cells 
"#8CC483" # SMCs 
"#D8D97A" # SCs 
"#538EB9" # Early GMCs
"#356493" # Dividing GMCs 
"#1D4573" # Early GCs 
"#0A2E57" # Late GCs




### Compare targetomes
library(tidyverse)
edge_table <- read.delim("E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/stom_edgeTable.tsv")

## add gene info
info <- read_delim("E:/Scseq_brachy_dev/brachy_gene_protein_info_Bd21-3_v1.2.csv", delim=";")
info$TG <- info$Bd21_3_id
info <- info %>% dplyr::select(-protein_id, -protein_length, -Protein_info_uniprot)
info <- unique(info)



## For gc_sc_lineage, cluster annotations
#"Stage 0-2" = "0" 
#"Early GMCs" = "1"
#"Dividing GMCs" = "2"
#"Early GCs" = "3"
#"Late GCs" = "4"
#"SMCs" = "5"
#"SCs" = "6"
#"Inter-specialized cells" = "7"

## how many genes are potential targets for one or more of the TFs in a given cluster?
edge_table_genes <- edge_table %>% filter(cluster == "Cluster_2")
length(unique(edge_table_genes$TG))


### filter edge table for several TFs in specific clusters
edge_table_small <- edge_table %>% filter(TF %in% c("BdiBd21-3.1G0240400", # MUTE
                                                    "BdiBd21-3.2G0300000") & # FAMA
                                            cluster %in% c("Cluster_2")) #dividing GMCs


## or not cluster-specific
edge_table_small <- edge_table %>% filter(TF %in% c("BdiBd21-3.1G0240400", # MUTE
                                                    "BdiBd21-3.2G0300000")) # FAMA


edge_table_small <- edge_table_small %>% dplyr::select(-borda_rank, -borda_clusterRank)
edge_table_small <- merge(edge_table_small, info, by = "TG", all.x=T)


## filter for one TF
MUTE <- edge_table_small %>% filter(TF == "BdiBd21-3.1G0240400")
FAMA <- edge_table_small %>% filter(TF == "BdiBd21-3.2G0300000")


## compare overlapping genes and unique genes
overlap <- merge(MUTE, FAMA[,c(1,2,4)], by = "TG") # 262 genes

overlap$weight_ratio <- overlap$weight.x/overlap$weight.y
overlap$log_weight_ratio <- log(overlap$weight_ratio)

## select genes that are much stronger correlated to TF1 (MUTE) than to TF2 (FAMA) or vice versa (log < -3 and > 3) and plot
diff_corr <- overlap %>% filter(log_weight_ratio < -3 | log_weight_ratio > 3)

## plot
diff_corr %>%
  ggplot(mapping = aes(x=reorder(TG, log_weight_ratio), y=log_weight_ratio))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x = NULL, y="Log weight ratio")

MUTE_up_strict <- diff_corr %>% filter(log_weight_ratio>0) # 38 genes
FAMA_up_strict <- diff_corr %>% filter(log_weight_ratio<0) # 25 genes


## genes exclusive to one of the TFs
exclusive <- merge(MUTE, FAMA[,c(1,2,4)], by = "TG", all=T) %>% filter(is.na(TF.x) | is.na(TF.y))

table(exclusive$TF.x) # 13 genes exclusive to MUTE
table(exclusive$TF.y) # 6 genes exclusive to FAMA

MUTE_excl <- exclusive %>% filter(!is.na(TF.x))
FAMA_excl <- exclusive %>% filter(!is.na(TF.y))

## merge with the genes that are rather co-expressed with one of the TFs from the overlap set
MUTE_pref_strict <- rbind(MUTE_up_strict[,1:21], MUTE_excl) # filtered (log_weight_ratios < -3 and > 3)   --> 51 genes
FAMA_pref_strict <- rbind(FAMA_up_strict[,1:21], FAMA_excl) # filtered (log_weight_ratios < -3 and > 3)   --> 31 genes






## Show the top X most common GO terms
## GO term info
go_terms_x <- read.table("E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis.csv") # all GO terms
go_terms_x <- read.table("E:/Scseq_brachy_dev/GRN_analysis/GO_terms_GRN-analysis_BP.csv") # biological process GO terms


## write out Table S7 of the paper
MUTE_up_strict$group <- "MUTE preferred"
MUTE_excl$group <- "MUTE only"
MUTE_excl$weight_ratio <- NA
MUTE_excl$log_weight_ratio <- NA
FAMA_up_strict$group <- "FAMA preferred"
FAMA_excl$group <- "FAMA only"
FAMA_excl$weight_ratio <- NA
FAMA_excl$log_weight_ratio <- NA

table_s7 <- rbind(MUTE_up_strict, MUTE_excl)
table_s7 <- rbind(table_s7, FAMA_up_strict)
table_s7 <- rbind(table_s7, FAMA_excl)

write.table(table_s7, "E:/Scseq_brachy_dev/GRN_analysis/Mini-EX/Mini-EX_final_stomatal-files_paper-revisions/output/regulons/TableS7.txt", sep = ";")




plot_topGOs_comp <- function(dataset, name, n = 25, plot_only="no"){
  if(dataset == "overlap"){
    new <- overlap
  }
  
  else {
    if(dataset == "exclusive") {
      new <- exclusive
    }
    
    else {
      if(dataset == "MUTE_up_strict") {
        new <- MUTE_up_strict
      }
        
      else {
        if(dataset == "FAMA_up_strict") {
          new <- FAMA_up_strict
        }
            
        else {
          if(dataset == "MUTE_pref_strict") {
            new <- MUTE_pref_strict
          }
            
          else {
            if(dataset == "FAMA_pref_strict") {
              new <- FAMA_pref_strict
            }
              
            else {
              if(dataset == "MUTE_excl") {
                new <- MUTE_excl
              }
              
              else {
                if(dataset == "FAMA_excl") {
                  new <- FAMA_excl
                }
              }
            }
          }
        }
      }
    }
  }
  
  new$Bd21_3_id <- new$TG
  new <- new %>% dplyr::select(-GO)
  
  
  new_go <- merge(new, go_terms_x, all.x=T, by = "Bd21_3_id")
  
  new_go <- unique(new_go[,c("Bd21_3_id", "GO", "term")])
  
  sum_nr_go <- new_go %>% group_by(GO, term) %>% summarise(freq = length(GO))
  
  sum_nr_go <- na.omit(sum_nr_go[order(sum_nr_go$freq, decreasing = T),])
  
  top20 <- as.data.frame(sum_nr_go[1:n,])
  
  
  top20 <- top20 %>% dplyr::select(-GO)
  
  plot <- top20 %>% 
    mutate(term = fct_reorder(term, top20$freq[order(-top20$freq)])) %>%
    ggplot(mapping = aes(x = term, y = freq)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(y="Counts in Targetome", x=NULL, title = name)
  
  if(plot_only == "no") {
    print(plot)
    return(as.data.frame(top20)) 
  }
  
  else {
    return(plot)
  }
}

### save info into new data frame and plot
overlap_topGOs <- plot_topGOs_comp(dataset = "overlap", name = "MUTE and FAMA \noverlapping targetome", n = 20)
exclusive_topGOs <- plot_topGOs_comp(dataset = "exclusive", name = "MUTE and FAMA \nexclusive targetome", n = 20)
MUTE_up_strict_topGOs <- plot_topGOs_comp(dataset = "MUTE_up_strict", name = "Shared TGs \nup in MUTE cells vs. FAMA cells, strict", n = 20)
FAMA_up_strict_topGOs <- plot_topGOs_comp(dataset = "FAMA_up_strict", name = "Shared TGs \nup in FAMA cells vs. MUTE cells, strict", n = 20)
MUTE_pref_strict_topGOs <- plot_topGOs_comp(dataset = "MUTE_pref_strict", name = "TGs preferentially regulated \n by MUTE vs. FAMA, strict", n = 20)
FAMA_pref_strict_topGOs <- plot_topGOs_comp(dataset = "FAMA_pref_strict", name = "TGs preferentially regulated \n by FAMA vs. MUTE, strict", n = 20)
MUTE_pref_strict_topGOs <- plot_topGOs_comp(dataset = "MUTE_pref_strict", name = "TGs preferentially regulated \n by MUTE vs. FAMA, strict", n = 25)
FAMA_pref_strict_topGOs <- plot_topGOs_comp(dataset = "FAMA_pref_strict", name = "TGs preferentially regulated \n by FAMA vs. MUTE, strict", n = 25)
MUTE_excl_topGOs <- plot_topGOs_comp(dataset = "MUTE_excl", name = "MUTE \nexclusive targetome", n = 20)
FAMA_excl_topGOs <- plot_topGOs_comp(dataset = "FAMA_excl", name = "MUTE \nexclusive targetome", n = 20)


##### visualize
MUTE_pref_strict_topGOs$assay <- "MUTE_pref_strict"
FAMA_pref_strict_topGOs$assay <- "FAMA_pref_strict"

final_GOs <- rbind(MUTE_pref_strict_topGOs, FAMA_pref_strict_topGOs)

final_GOs %>% 
  mutate(assay = factor(assay, levels = c("MUTE_pref_strict", "FAMA_pref_strict"))) %>%
  mutate(term = fct_reorder(term, final_GOs$freq[order(-final_GOs$freq)])) %>%
  ggplot(mapping = aes(x = assay, y = term, size = freq)) +
  geom_point() +
  theme_bw() +
  labs(y= NULL, x= NULL)
