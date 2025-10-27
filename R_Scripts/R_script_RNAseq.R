###########
#RNASeq Reanalysis with Bd21_3
###########

#comparison between whole leaf zone and protoplasted leaf zone data set from Ines

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicAlignments")


library(DESeq2)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(edgeR)
library(readxl)

setwd("E:/Scseq_brachy_dev/RNA-bulk-seq/")


annot <- makeTxDbFromGFF("../JGI_files/genome_Bd21-3_v1.1/BdistachyonBd21_3_460_v1.1.gene_exons.gff3", format="gff3")

exons_by_genewholevsproto <- exonsBy(annot, by="gene") # all exons are grouped by 'gene'

#specify the bam files used for counting
bam_fileswholevsproto <- list.files("Y:/2025_Berg-et-al_ScSeq_GEO_data/bulk_RNASeq_data/", pattern="bam$", full=TRUE ) # put all the bam files in one folder and give here the path to that folder

bam_fileswholevsproto #this should list the number of files used for generating the list

bam_listwholevsproto <- BamFileList(bam_fileswholevsproto, yieldSize =100000)
bam_listwholevsproto

BiocParallel::register(BiocParallel::SerialParam())


### SumOver ---------------------------------------------------------------------------------------------------------

#this step takes longer
Totalread_countswholevsproto <- summarizeOverlaps(exons_by_genewholevsproto, bam_listwholevsproto,
                                                mode="Union",
                                                singleEnd=FALSE,
                                                ignore.strand=TRUE,
                                                fragments=TRUE ) #count the reads

Totalread_countswholevsproto
head(Totalread_countswholevsproto)
nrow(Totalread_countswholevsproto)
ncol(Totalread_countswholevsproto)
head(assay(Totalread_countswholevsproto))
colSums(assay(Totalread_countswholevsproto)) 


#save countdata as csv:
class(Totalread_countswholevsproto)
write.table(assays(Totalread_countswholevsproto)$counts,"wholevsprotoBd21_3_sumOver_raw.csv", sep=",", row.names=FALSE) 

#use assays(se)$counts to get count data out of multilayer SummarizedExperiment table
#this gives you a csv table of the name of your choosing (e.g., wholevsprotoBd21_3_sumOver_raw.csv)




### TPM normalization ---------------------------------------------------------------------

TPMwholevsproto <- assays(Totalread_countswholevsproto)$counts
head(TPMwholevsproto)
#attach gene length from gene file
#be aware that the file hast to be directly in the folder of the working directory, not in a subfolder
#!!!just change the name of the file from gff3.gz to gtf and then read it in
genewholevsproto<-read.table("E:/Scseq_brachy_dev/JGI_files/genome_Bd21-3_v1.1/BdistachyonBd21_3_460_v1.1.gene.gtf.gz", sep="\t", fill=T,quote = "", 
                          row.names = NULL, 
                          stringsAsFactors = FALSE)
head(genewholevsproto)

colnames(TPMwholevsproto) <- c("Protoplasted_1", "Protoplasted_2", "Protoplasted_3", "Non-protoplasted_1", "Non-protoplasted_2", "Non-protoplasted_3") 
#this line should contain the correct sample names and number of bamfiles

#extract locus name and attach as column
genewholevsproto$locus<-substr(genewholevsproto$V9,4,22) # !V9,4,22 specifies whole gene name from letter 4 to 22
#!!!check if your gene name has the same length

genewholevsproto$length<-abs(genewholevsproto$V5 - genewholevsproto$V4) #get gene length # abs makes an absolute number, so negative numbers don't exist
mRNAwholevsproto<-genewholevsproto[genewholevsproto$V3 == "mRNA", ] #keep only rows that say "mRNA". This is called subsetting
mRNAwholevsproto<-mRNAwholevsproto[ !duplicated(mRNAwholevsproto$locus), ] #keep only primary transcript (=1st row) in genes with alternative splicing

TPMwholevsproto<-as.data.frame((TPMwholevsproto)) 
TPMwholevsproto$locus <- rownames(TPMwholevsproto)

TPMwholevsproto$length<-mRNAwholevsproto[match(TPMwholevsproto$locus, mRNAwholevsproto$locus),11] #add length

#1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your per million scaling factor.
#3. Divide the RPK values by the per million scaling factor. This gives you TPM.


#wtprotoplasted 1
TPMwholevsproto$Protoplasted_1_RPK <- ((TPMwholevsproto$Protoplasted_1*1000) / TPMwholevsproto$length)
wholevsprotoSFwtproto1 <-(sum(TPMwholevsproto$Protoplasted_1_RPK, na.rm = TRUE))/1000000
TPMwholevsproto$Protoplasted_1_TPM <-TPMwholevsproto$Protoplasted_1_RPK/wholevsprotoSFwtproto1

#wtprotop2
TPMwholevsproto$Protoplasted_2_RPK <- ((TPMwholevsproto$Protoplasted_2*1000) / TPMwholevsproto$length)
wholevsprotoSFwtproto2 <-(sum(TPMwholevsproto$Protoplasted_2_RPK, na.rm = TRUE))/1000000
TPMwholevsproto$Protoplasted_2_TPM <-TPMwholevsproto$Protoplasted_2_RPK/wholevsprotoSFwtproto2

#wtprotop3
TPMwholevsproto$Protoplasted_3_RPK <- ((TPMwholevsproto$Protoplasted_3*1000) / TPMwholevsproto$length)
wholevsprotoSFwtproto3 <-(sum(TPMwholevsproto$Protoplasted_3_RPK, na.rm = TRUE))/1000000
TPMwholevsproto$Protoplasted_3_TPM <-TPMwholevsproto$Protoplasted_3_RPK/wholevsprotoSFwtproto3

#wtwholeLZ (non-protoplasted)
TPMwholevsproto$`Non-protoplasted_1_RPK` <- ((TPMwholevsproto$`Non-protoplasted_1`*1000) / TPMwholevsproto$length)
wholevsprotoSFwtwhole1 <-(sum(TPMwholevsproto$`Non-protoplasted_1_RPK`, na.rm = TRUE))/1000000
TPMwholevsproto$`Non-protoplasted_1_TPM` <-TPMwholevsproto$`Non-protoplasted_1_RPK`/wholevsprotoSFwtwhole1

#wtwholeLZ2
TPMwholevsproto$`Non-protoplasted_2_RPK` <- ((TPMwholevsproto$`Non-protoplasted_2`*1000) / TPMwholevsproto$length)
wholevsprotoSFwtwhole2 <-(sum(TPMwholevsproto$`Non-protoplasted_2_RPK`, na.rm = TRUE))/1000000
TPMwholevsproto$`Non-protoplasted_2_TPM` <-TPMwholevsproto$`Non-protoplasted_2_RPK`/wholevsprotoSFwtwhole2

#wtwholeLZ3
TPMwholevsproto$`Non-protoplasted_3_RPK` <- ((TPMwholevsproto$`Non-protoplasted_3`*1000) / TPMwholevsproto$length)
wholevsprotoSFwtwhole3 <-(sum(TPMwholevsproto$`Non-protoplasted_3_RPK`, na.rm = TRUE))/1000000
TPMwholevsproto$`Non-protoplasted_3_TPM` <-TPMwholevsproto$`Non-protoplasted_3_RPK`/wholevsprotoSFwtwhole3



#average TPM values, in the brackets, always keep the number of variables equal to the number of replicate of each sample
TPMwholevsproto$Protoplasted_TPM <-rowMeans(TPMwholevsproto [,c("Protoplasted_1_TPM","Protoplasted_2_TPM","Protoplasted_3_TPM")])
TPMwholevsproto$Non_protoplasted_TPM <-rowMeans(TPMwholevsproto [,c("Non-protoplasted_1_TPM","Non-protoplasted_2_TPM","Non-protoplasted_3_TPM")])


TPMwholevsproto2 <- TPMwholevsproto[,c(7, 1:6, 8:22)]

write.table(TPMwholevsproto2,'rawcounts_sumOver_TPM.csv', sep=',', row.names=FALSE) #creates sumOver_TPM.csv file

numZero <- colSums(TPMwholevsproto == 0, na.rm = T)
numZero
numNA <- colSums(is.na(TPMwholevsproto))
numNA
numnotZero <- colSums(TPMwholevsproto > 0, na.rm = T)
numnotZero 


sampleInfowholevsproto <- read.table ("wholevsproto.csv", sep=",", header=TRUE) 
#sample info shows the experimental setup #sample info says how many levels: different dev stages, different genotypes

colnames(Totalread_countswholevsproto) <- sampleInfowholevsproto[,1]
colData (Totalread_countswholevsproto)
dim(Totalread_countswholevsproto)
head(Totalread_countswholevsproto)
head(sampleInfowholevsproto)
dim(sampleInfowholevsproto)


##recreate SE file from TPM
View(TPMwholevsproto)
countsMTRwholevsprotop <- TPMwholevsproto[,1:6] #[,1:x], where x is the number of files used for making bamlist and total read counts
data.se_whole_vs_proto <- SummarizedExperiment(list(counts=as.matrix(countsMTRwholevsprotop)))
sampleInfowholevsproto<-DataFrame(sampleInfowholevsproto) 
sampleInfowholevsproto

Totalread_countswholevsprotoIdx_wholevsproto <- match(colnames(data.se_whole_vs_proto), sampleInfowholevsproto$run) #6 rows with 3 columns
head( cbind( colData(data.se_whole_vs_proto), sampleInfowholevsproto[ Totalread_countswholevsprotoIdx_wholevsproto, ] ) )
colData(data.se_whole_vs_proto) <- cbind( colData(data.se_whole_vs_proto), sampleInfowholevsproto)
colData(Totalread_countswholevsproto) 
###




### DESeq2 analysis -----------------------------------------------------------------------------------

library("DESeq2")

ddsFull_wholevsproto <- DESeqDataSet(data.se_whole_vs_proto, design= ~ StageII) 
ddsFull_wholevsproto 
head(ddsFull_wholevsproto)
dds_wholevsproto<-ddsFull_wholevsproto
dds_wholevsproto

dds_wholevsproto$StageII <- relevel(dds_wholevsproto$StageII, "leaf")


dds_wholevsproto
dim(dds_wholevsproto)

# Run the DESeq pipeline
design(dds_wholevsproto)
dds_wholevsproto <- DESeq(dds_wholevsproto)
dds_wholevsproto
res_wholevsproto <- results(dds_wholevsproto)
res_wholevsproto
mcols(res_wholevsproto, use.names=TRUE)
write.table(res_wholevsproto, "wholevsproto_Bd21_3_sumOver_raw_deseq2.csv", sep=",") 
sum( res_wholevsproto$padj< 0.1, na.rm=TRUE ) #gives you the number of all genes with a p adjusted values below 0.1



##########################################################################
#list of DESeq2 DEGs
deseq_dev_wholevsproto <- res_wholevsproto[ which(res_wholevsproto$padj < 0.1), ]
deseq_dev_wholevsproto$gene <- row.names(deseq_dev_wholevsproto)
deseq_dev_wholevsproto<-as.data.frame(deseq_dev_wholevsproto) #convert to data.frame
class(deseq_dev_wholevsproto)
colnames(deseq_dev_wholevsproto)[colnames(deseq_dev_wholevsproto) == 'gene'] <- 'Bd21_3'
View(deseq_dev_wholevsproto)
write.table(deseq_dev_wholevsproto, "wholevsproto_Bd21_3_sumOver_raw_deseq2_padj.csv", sep=",") 

deseq_dev_wholevsproto <- read.table("wholevsproto_Bd21_3_sumOver_raw_deseq2_padj.csv", sep =",")

#attach the info and synonym file to the Bd21_3_sumOver_raw_deseq2_dev_padj file
info_wholevsproto<-read.table("E:/Scseq_brachy_dev/JGI_files/genome_Bd21-3_v1.1/BdistachyonBd21_3_460_v1.1.annotation_info.txt", sep="\t", header = FALSE, fill=T,quote = "", 
                 row.names = NULL, 
                 stringsAsFactors = FALSE)

colnames(info_wholevsproto)[colnames(info_wholevsproto) == 'V2'] <- 'Bd21_3'
synonym_wholevsproto<-read.table("E:/Scseq_brachy_dev/JGI_files/genome_Bd21-3_v1.1/BdistachyonBd21_3_460_v1.1.synonym.txt", sep="\t", header=FALSE, fill=T,quote = "", 
                    row.names = NULL, 
                    stringsAsFactors = FALSE)
synonym_wholevsproto$V1<-substr(synonym_wholevsproto$V1,1,19) #remove .xyz transcript tag ##make sure you do not delete too little or too many
colnames(synonym_wholevsproto)[colnames(synonym_wholevsproto) == 'V1'] <- 'Bd21_3'
colnames(synonym_wholevsproto)[colnames(synonym_wholevsproto) == 'V2'] <- 'Bradi'
synonym_wholevsproto_withoutDupli<-synonym_wholevsproto[ !duplicated(synonym_wholevsproto$'Bradi'), ] #remove duplicates

info_all_wholevsproto<-merge(info_wholevsproto, synonym_wholevsproto, by='Bd21_3') #add Bradi info
info_all_wholevsproto_changed_order<-info_all_wholevsproto[c(1,17,2:16)] #change order
info_all_wholevsproto_withoutDupli<-info_all_wholevsproto_changed_order[ !duplicated(info_all_wholevsproto_changed_order$'Bradi'), ] #remove duplicates


wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge<-dplyr::left_join(deseq_dev_wholevsproto, info_all_wholevsproto, by = "Bd21_3") #this command merges the tables
wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge<-
  wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge[ !duplicated(wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge$'Bd21_3'), ] 
#the command above keeps only primary transcript (=1st row) in genes with alternative splicing

wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge<-wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge[c(7,8,1:6,9:23)] #change order of columns


write.table(wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge, 
            "wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge(info).csv", sep=",", row.names = FALSE)

write_excel_csv(wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge, "wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge(info).xlsx")

View(wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge)



### PLOTTING GRAPHS --------------------------------------------------------------------------------------------- 

wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge <- read_excel("wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge(info).xlsx")

### Volcano plot with ggplot
library(dplyr)
library(ggplot2)

#mark genes that are up- or downregulated
wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge <- wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge %>%
  mutate(gene_type = case_when(log2FoldChange >= 3 & padj < 0.05 ~ "up",
                               log2FoldChange <= -3 & padj < 0.05 ~ "down",
                               TRUE ~ "ns"))   
#plot
ggplot(wholevsproto_Bd21_3_sumOver_raw_deseq2_padj_merge, aes(x = log2FoldChange, y = -log10(padj), colour = gene_type)) + 
  geom_point(size=1, alpha=0.3) +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed") +
  xlim(c(-13.5,13.5)) +
  scale_colour_manual(values = c("forestgreen", "lightgrey", "purple")) +
  theme_classic()
