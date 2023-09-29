### PARAGON ANALYSES - Obiol Gene Groups ###
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
library(ggridges)

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
# dfWide <- subset(dfWide,grepl("Eukaryot",dfWide$Taxonomy))
#dfWide <- subset(dfWide, grepl("Fungi",dfWide$Taxonomy)|grepl("Ascomycot",dfWide$Taxonomy))
dfWide <- subset(dfWide, grepl("Discob",dfWide$Taxonomy))


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

# Fungi and protists
# dfFin <- subset(dfFin,grepl("Fungi",dfFin$Taxonomy)|grepl("Ascomycot",dfFin$Taxonomy))
#dfProt <- subset(dfFin,!grepl("Fungi",dfFin$Taxonomy)& !grepl("Ascomycot",dfFin$Taxonomy))
# dfFin <- subset(dfFin,grepl("Discoba",dfFin$Taxonomy))

# Load in biomarkers
ko <- read.csv("KEGG_KO_Obiol.csv",header=TRUE)
ko <- ko[c(1:3)]
colnames(ko)[1] <- "KEGG"
ko$KEGG <- str_remove_all(ko$KEGG,"ko:")

dfFin$Name <- NULL
dfFin$Taxonomy <- NULL
dfFin <- left_join(dfFin,ko)

dfFin <- subset(dfFin,!is.na(Category))
dfFin$KEGG <- NULL
dfFin$Gene <- NULL
dfFin <- dfFin %>% group_by(Category) %>% summarize_all(sum) %>% as.data.frame() 
dfFin <- melt(dfFin,id.vars="Category")

meanz <- dfFin %>% group_by(Category) %>% summarize(m=mean(value),s=sd(value))
dfFinA <- left_join(dfFin,meanz)
dfFinA$z <- (dfFinA$value-dfFinA$m)/dfFinA$s
dfFinA$value <- NULL
dfFinA$m <- NULL
dfFinA$s <- NULL
colz <- colsplit(dfFinA$variable,"_",c("Sample","Experiment"))
dfFinA <- cbind(dfFinA,colz)
dfFinA$variable <- NULL
dfFinA$Sample <- ifelse(grepl("Water",dfFinA$Sample),"Water Column",dfFinA$Sample)

dfFinA <- dfFinA %>% group_by(Sample,Category) %>% summarize(m=mean(z)) %>% as.data.frame()

dfFinA %>% ggplot(aes(x=Sample,y=Category,fill=m))+geom_tile(color="grey")+theme_classic(base_size=16)+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_fill_gradient2(high="red",low="blue",name="Z-score transformed\ntranscript abundance")+ggtitle("")+ylab("KEGG KO")


colz <- colsplit(dfFin$variable,"_",c("Sample","Experiment"))
dfFin <- cbind(dfFin,colz)
meanz <- dfFin %>% group_by(Category,Experiment) %>% summarize(m=mean(value),s=sd(value))
dfFinB <- left_join(dfFin,meanz)
dfFinB$z <- (dfFinB$value-dfFinB$m)/dfFinB$s
dfFinB$value <- NULL
dfFinB$m <- NULL
dfFinB$s <- NULL
dfFinB$variable <- NULL
dfFinB$Sample <- ifelse(grepl("Water",dfFinB$Sample),"Water Column",dfFinB$Sample)
dfFinB <- dfFinB %>% group_by(Sample,Category,Experiment) %>% summarize(m=mean(z)) %>% as.data.frame()

dfFinB %>% ggplot(aes(x=Sample,y=Category,fill=m))+geom_tile(color="grey")+theme_classic(base_size=16)+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_fill_gradient2(high="red",low="blue",name="Z-score transformed\ntranscript abundance")+ggtitle("")+ylab("KEGG KO")+facet_wrap(~Experiment)

ggsave("KEGG_HeatMap.pdf",width=8,height=6)
