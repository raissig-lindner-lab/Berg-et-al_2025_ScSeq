library(readxl)
library(tidyverse)
library(MetBrewer)
library(agricolae)
#install.packages('devtools')
#devtools::install_github("lbmountain/licornetics")
library(licornetics)
library(ggpubr)



##### gene models
library(rtracklayer)

pme53_gff <- as.data.frame(readGFF("D:/OneDrive - Universitaet Bern/Scseq_brachy_dev/Gene_models/PME53 BdiBd21-3.2G0255500 v1.2.gff"))
gras32_gff <- as.data.frame(readGFF("D:/OneDrive - Universitaet Bern/Scseq_brachy_dev/Gene_models/GRAS32 BdiBd21-3.1G0657800 v1.2.gff"))
gff_file <- pme53_gff
gff_file <- gras32_gff

genebackbone <- gff_file %>% dplyr::filter(type %in% c("gene"))
gene3p <- gff_file %>% dplyr::filter(type %in% c("three_prime_UTR"))
gene5p <- gff_file %>% dplyr::filter(type %in% c("five_prime_UTR"))
exons <- gff_file %>% dplyr::filter(type %in% c("CDS"))
mutation <- gff_file %>% dplyr::filter(type %in% c("Polymorphism"))

ggplot() + geom_segment(genebackbone, mapping = aes(x = start, y = 0, xend = end, yend = 0), size = 1) + 
  geom_rect(gene3p, mapping = aes(xmin = start, ymin = -0.1, xmax = end, ymax = 0.1), fill = "white", colour = "black") + 
  geom_rect(gene5p, mapping = aes(xmin = start, ymin = -0.1, xmax = end, ymax = 0.1), fill = "white", colour = "black") + 
  geom_text(gene5p[1,], mapping = aes(x = start+100, y = -0.5, label = paste("5'-UTR"))) + 
  geom_rect(exons, mapping = aes(xmin = start, ymin = -0.1, xmax = end, ymax = 0.1), colour = "black", fill = "grey") + 
  geom_rect(mutation, mapping = aes(xmin=start, ymin=-0.15, xmax=end, ymax=0.15), colour = "black", fill = "black") +
  scale_y_continuous(limits = c(-1, 1)) + 
  scale_x_continuous(breaks = c(0, 100, 1000, 2000, 3000)) +
  theme_classic() + theme(strip.background = element_rect(colour = NA, fill = NA)) + 
  theme(axis.text.y = element_blank()) + 
  theme(axis.ticks.y = element_blank()) + 
  labs(x = "Gene position", y = NULL)




##### Physiology pme53-like ---------------------------------------------------------------------------------------------
#### Licor plots
p1 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"))

p2 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), 
                 type = "relgsw")

p3 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), 
                 type = "A")

p4 <- licorplots(c("wt", "pme53"), timeframe = 15:60, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), 
                 type = "WUE")

ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), legend = "bottom", common.legend = T)



### corrected by stomatal density
p1.2 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"),
                 stomden = c(79.5, 58.2))

p2.2 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), type = "relgsw",
                 stomden = c(79.5, 58.2))

p3.2 <- licorplots(c("wt", "pme53"), timeframe = 15:75, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), type = "A",
                 stomden = c(79.5, 58.2))

p4.2 <- licorplots(c("wt", "pme53"), timeframe = 15:60, timestamps = c(20,40,60), 
                 remove_outliers = "yes", colours = met.brewer("Hokusai3", n=3), legend_labels = c("WT", "pme53"), type = "WUE",
                 stomden = c(79.5, 58.2))

ggarrange(p1.2, p2.2, p3.2, p4.2, labels = c("A", "B", "C", "D"), legend = "bottom", common.legend = T)
# save as 6x6 pdf


##### Stomatal density pme53-like -------------------------------------------------------------------------------------

##load xlsx
stomden <- read_excel("pme53_stomatal_density.xlsx", sheet=1, col_names = T)


##calculate mean stomatal density for each individual
stomden1 <- stomden %>% group_by(line, identifier, individual) %>% summarise(mean_dens= mean(stomatal_density),
                                                                             sd_dens = sd(stomatal_density),
                                                                             stomata = sum(stomata))

stomden2 <- stomden1 %>% group_by(line) %>% summarise(mean_dens = mean(mean_dens),
                                                      mean_sd = mean(sd_dens))
stomden2
#  line  mean_dens mean_sd
#1 WT         79.5    16.5
#2 pme53      58.2    11.3


### test for significant differences
### all values
WT <- stomden %>% filter(line == "WT")
pme53 <- stomden %>% filter(line == "pme53")

t.test(WT$stomatal_density, pme53$stomatal_density) # significantly different!


### means of individuals
WT <- stomden1 %>% filter(line == "WT")
pme53 <- stomden1 %>% filter(line == "pme53")

t.test(WT$mean_dens, pme53$mean_dens) # significantly different!




##create data frame with significance levels
sigden <- data.frame(line=c("WT", "pme53"),
                     significance=c("a", "b"))


##reorder data
stomden1$line <- ordered(stomden1$line , levels=c("WT", "pme53"))

##plot stomatal density data
ggplot(stomden1, mapping=aes(x=line, y=mean_dens))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(sigden, mapping=aes(y=100, label=paste(significance)), show.legend = F)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  scale_y_continuous(limits = c(0, 100))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x=NULL, y=expression(paste("Stomatal density [stomata mm"^-2, "]")))

















##### MorphoGraphX gras32 -----------------------------------------------------------------------------------------

### load data (it needs to be in the current working directory!)
mgx_all <- data.frame(Genotype=NA, File=NA, Label=NA, area=NA)

identifier <- c("WT", "hom")

for (i in identifier) {
  print(i)
  
  files <- dir(pattern=i)
  for (onefile in files) {
    print(onefile)
    new_file <- read_delim(onefile, delim = ",", col_names = T)
    
    new_file$area <- new_file$`Area (µm²)`
    new_file$File <- print(onefile)
    new_file$Genotype <- print(i)
    
    new_file <- new_file %>% select(Genotype, File, Label, area)
    
    new_file <- na.omit(new_file)
    
    mgx_all <- na.omit(rbind(mgx_all, new_file))
    
    mgx_all <- mgx_all %>% filter(area > 8)
  }
}

#mgx <- read_delim("WT_1_cropped.csv", delim = ",") # to read in just a single file


for(row in 1:nrow(mgx_all)) {
  if(mgx_all$Genotype[row] == "hom") {
    mgx_all$Genotype[row] <- "bdgras32"
  }
}


### How many cells were counted per Individual?
summary(mgx_all) # min cell area was 9.74 µm2, max cell area was 234.51 µm2

table(mgx_all$Genotype) # 1268 cells were counted in total for gras32-like, 938 for WT
table(mgx_all$File) 
# gras32: between 203 cells and 298 cells
# WT: between 157 cells and 204 cells


##calculate cell counts and mean area for each individual
mgx_sum <- mgx_all %>% group_by(Genotype, File) %>% summarise(cells = length(area),
                                                              mean_area = mean(area))

### average cell counts and area per genotype
mgx_sum %>% group_by(Genotype) %>% summarise(mean_cells = mean(cells),
                                             mean_area_all = mean(mean_area))
#Genotype mean_cells mean_area_all
#WT             188.          73.1
#bdgras32       254.          52.1





### test for significant differences


## Two-sided Student's t-Test
# for means of individuals
WT <- mgx_all %>% filter(Genotype=="WT")
gras32 <- mgx_all %>% filter(Genotype=="bdgras32")

t.test(WT$area, gras32$area) # significantly different!


# for means of individuals
WT <- mgx_sum %>% filter(Genotype=="WT")
gras32 <- mgx_sum %>% filter(Genotype=="bdgras32")

t.test(WT$cells, gras32$cells) # significantly different!
t.test(WT$mean_area, gras32$mean_area) # significantly different!


##create data frame with significance levels
signif <- data.frame(Genotype=c("WT", "bdgras32"),
                     sig_cells=c("a", "b"),
                     sig_area=c("a", "b"))


### Reorder data
mgx_sum$Genotype <- ordered(mgx_sum$Genotype, levels=c("WT", "bdgras32"))

### Plot cell count per Genotype
ggplot(mgx_sum, mapping=aes(x=Genotype, y=cells))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(signif, mapping=aes(y=320, label=paste(sig_cells)), show.legend = F)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  scale_y_continuous(limits = c(0, 320))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x=NULL, y=expression(paste("Mean cell numbers")))


### Plot mean cell area per Genotype
ggplot(mgx_sum, mapping=aes(x=Genotype, y=mean_area))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(signif, mapping=aes(y=90, label=paste(sig_area)), show.legend = F)+
  geom_jitter(width=0.1, height = 0, alpha = 0.5)+
  scale_y_continuous(limits = c(0, 90))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x=NULL, y=expression(paste("Mean cell area [µm"^2,"]")))





### Reorder data
mgx_all$Genotype <- ordered(mgx_all$Genotype, levels=c("WT", "bdgras32"))

### Plot cell areas per Individual
ggplot(mgx_all, mapping=aes(x=Genotype, y=area))+
  geom_boxplot(outlier.shape = NA, colour = "black")+
  geom_text(signif, mapping=aes(y=250, label=paste(sig_area)), show.legend = F)+
  geom_jitter(width=0.1, height = 0, alpha = 0.3)+
  scale_y_continuous(limits = c(0, 250))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x=NULL, y=expression(paste("Cell area [µm"^2,"]")))
