### PARAGON ANALYSES - PCoA ###
### LAST UPDATED: 9/29/23 ###

### DATA WRANGLING / NORMALIZATION ###
# Packages
library(tidyverse)
library(reshape2)
library(edgeR)
library(randomcoloR)
library(patchwork)
library(vegan)
library(ape)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(clusterProfiler)
library(ggupset)
library(stringr)
library(stringi)

# Load in data
eukegg <- read.csv("eggnog_eukulele_fungi.csv",header=TRUE,row.names=1)
sal <- read.csv("salmon_wide.csv",header=TRUE,row.names=1)

# Wrangle data
eukegg$transcript_name <- str_remove_all(eukegg$transcript_name,".p1")
full <- left_join(sal,eukegg)
full <- full[c(1,17:19,2:16)]
full$CAZy <- NULL
full$KEGG_ko <- str_remove_all(full$KEGG_ko,"ko:")
full$KEGG_ko <- ifelse(is.na(full$KEGG_ko),"",full$KEGG_ko)
full$KEGG_ko <- str_split(full$KEGG_ko,",")
dfWide <- unnest(full,KEGG_ko)
colnames(dfWide)[1:3] <- c("Name","Taxonomy","KEGG")
dfWide <- subset(dfWide,grepl("Eukaryot",dfWide$Taxonomy))

# Separate Exp #1 and Exp #2 from Exp #3 (Experiment #3 was sequenced separately)
dfExp12 <- dfWide[,c("Name","Taxonomy","KEGG","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3 <- dfWide[,c("Name","Taxonomy","KEGG","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 
dfExp12b <- dfExp12
dfExp3b <- dfExp3

# Normalize Exp 1&2
dgeExp12 <- DGEList(counts=dfExp12b[4:13],genes=dfExp12b[1:3],group=c(rep("T0_Exp1",1),rep("T3_Exp1",1),rep("T6_Exp1",1),rep("T0_Exp2",1),rep("T3_Exp2",1),rep("T6_Exp2",1),rep("Water_Exp1",2),rep("Water_Exp2",2)))

dgeExp12$samples # Visualize groups

Exp12Norm <- calcNormFactors(dgeExp12,method = "TMM")
Exp12Cpm <- cpm(Exp12Norm, normalized.lib.sizes=TRUE, log=FALSE) 
Exp12Cpm <-as.data.frame(Exp12Cpm)  
Exp12Cpm<-data.frame(Exp12Norm$genes,Exp12Cpm)

# Normalize Exp 3
dgeExp3 <- DGEList(counts=dfExp3b[4:8],genes=dfExp3b[1:3],group=c(rep("T0_Exp3",1),rep("T3_Exp3",1),rep("T6_Exp3",1),rep("Water_Exp3",2)))

dgeExp3$samples # Visualize groups

Exp3Norm <- calcNormFactors(dgeExp3,method = "TMM")
Exp3Cpm <- cpm(Exp3Norm, normalized.lib.sizes=TRUE, log=FALSE) 
Exp3Cpm <-as.data.frame(Exp3Cpm)  
Exp3Cpm<-data.frame(Exp3Norm$genes,Exp3Cpm)

# Combine Normalized Exp 1&2 data with Normalized Exp 3 data
dfFin <- left_join(Exp12Cpm,Exp3Cpm)

### PCA Plot ###
dfPCA <- dfFin
dfPCA$Name <- NULL
dfPCA$Taxonomy <- NULL
dfPCA <- subset(dfPCA,KEGG!="")

dfPCA <- dfPCA %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()
rownames(dfPCA) <- dfPCA$KEGG
dfPCA$KEGG <- NULL
dfPCA <- as.data.frame(t(dfPCA))
rowz <- colsplit(rownames(dfPCA),"_",c("Sample","Exp"))
rowz$Sample <- ifelse(grepl("Water",rowz$Sample),"Water",rowz$Sample)
dfPCA <- cbind(dfPCA,rowz)
dfPCA <- dfPCA %>% group_by(Sample,Exp) %>% summarize_all(mean) %>% as.data.frame()
rownames(dfPCA) <- paste(dfPCA$Sample,dfPCA$Exp,sep="_")
dfPCA$Sample <- NULL
dfPCA$Exp <- NULL
dfPCA <- as.data.frame(t(dfPCA))
dfPCA <- subset(dfPCA,rowSums(dfPCA)!=0)
dfPCA <- as.data.frame(t(dfPCA))

distDf <- vegdist(dfPCA,  method = "bray")
pcoaOut <- pcoa(distDf)
pcoaOut$values
pcoaPlot <- data.frame(pcoaOut$vectors)
names <- colsplit(rownames(pcoaPlot),"_",c("Sample","Experiment"))
pcoaPlot <- cbind(pcoaPlot,names)
pcoaPlot$Sample <- ifelse(grepl("Water",pcoaPlot$Sample),"Water Column",pcoaPlot$Sample)
pcoaPlot$Sample <- ifelse(pcoaPlot$Sample=="RotT0","Net Trap T0",pcoaPlot$Sample)
pcoaPlot$Sample <- ifelse(pcoaPlot$Sample=="RotT3","Net Trap T3",pcoaPlot$Sample)
pcoaPlot$Sample <- ifelse(pcoaPlot$Sample=="RotT6","Net Trap T6",pcoaPlot$Sample)
pcoaPlot$Experiment <- ifelse(pcoaPlot$Experiment=="Exp1","Experiment #1 - 2021",pcoaPlot$Experiment)
pcoaPlot$Experiment <- ifelse(pcoaPlot$Experiment=="Exp2","Experiment #2 - 2021",pcoaPlot$Experiment)
pcoaPlot$Experiment <- ifelse(pcoaPlot$Experiment=="Exp3","Experiment #3 - 2022",pcoaPlot$Experiment)

keggPlot <- ggplot(pcoaPlot,aes(x=Axis.1,y=Axis.2))+geom_point(aes(fill=Sample,shape=Experiment),color="black",size=2)+scale_shape_manual(values=c(21,22,24))+scale_fill_manual(values=c("indianred","darkgoldenrod1","forestgreen","dodgerblue"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+xlab("36.0%")+ylab("19.2%")+ggtitle("")
keggPlot
ggsave("PCoA_KEGG.pdf")
