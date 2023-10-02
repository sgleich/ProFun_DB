### PARAGON ANALYSES - FUNGI VS. NO FUNGI DB ###
### LAST UPDATED: 10/2/23 ###

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
# eukegg <- read.csv("eggnog_eukulele_fungi.csv",header=TRUE,row.names=1)
eukegg <- read.csv("eggnog_eukulele.csv",header=TRUE,row.names=1)
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

dfFin$Tax <- ifelse(grepl("Fungi",dfFin$Taxonomy)|grepl("Ascomycot",dfFin$Taxonomy),"Fungi","Other Eukaryote")
unique(dfFin$Tax)

dfFin$Name <- NULL
dfFin$Taxonomy <- NULL
dfFin$KEGG <- NULL

dfMelt <- melt(dfFin,id.vars="Tax")
dfMelt <- dfMelt %>% group_by(variable,Tax) %>% summarize(s=sum(value)) %>% as.data.frame()
colz <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
colz$Sample <- ifelse(grepl("Water",colz$Sample),"Water Column",colz$Sample)
dfMelt <- cbind(dfMelt,colz)
dfMelt <- dfMelt %>% group_by(Sample,Experiment,Tax) %>% summarize(m=mean(s)) %>% as.data.frame()

dfMelt$Sample <- ifelse(dfMelt$Sample=="RotT0","Net Trap T0",dfMelt$Sample)
dfMelt$Sample <- ifelse(dfMelt$Sample=="RotT3","Net Trap T3",dfMelt$Sample)
dfMelt$Sample <- ifelse(dfMelt$Sample=="RotT6","Net Trap T6",dfMelt$Sample)

dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp1","Experiment #1\n2021",dfMelt$Experiment)
dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp2","Experiment #2\n2021",dfMelt$Experiment)
dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp3","Experiment #3\n2022",dfMelt$Experiment)

dfMelt$Sample <- factor(dfMelt$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

woFun <- ggplot(dfMelt,aes(x=Sample, y=m,fill=Tax))+geom_bar(stat="identity",position="fill",color="black")+theme_classic()+facet_wrap(~Experiment)+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Relative Transcript Abundance")+scale_fill_manual(name="Taxonomic Group",values=c("black","grey"))+ggtitle("MMETSP")

woFun+wFun+plot_layout(guides="collect",nrow=2)
ggsave("MMETSP_Comparison.pdf",width=9,height=8)
