### PARAGON ANALYSES - Cathepsin Plot ###
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
library(ggtext)

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

# Cathepsin proteases
cath <- c("K01275","K01319","K01363","K01365","K01366","K01368","K01371","K01373","K01374","K01375","K01379","K01382","K08568","K08569","K09599","K09600","K09601","K09616","K13289","K24378","K24379","K26838")

dfFin$Name <- NULL
dfFin$Taxonomy <- NULL
dfFin <- subset(dfFin,KEGG %in% cath)

# Find most abundant cathepsins
find <- dfFin %>% group_by(KEGG) %>% summarize_all(sum)
find$sum <- rowSums(find[2:16])
find <- find %>% arrange(desc(sum))
cath2 <- find$KEGG[1:12]

dfMelt <- melt(dfFin,id.vars="KEGG")
`%ni%` <- Negate(`%in%`)
dfMelt$KEGG <- ifelse(dfMelt$KEGG %ni% cath2,"Other Cathepsin",dfMelt$KEGG)
dfMelt <- dfMelt %>% group_by(KEGG,variable) %>% summarise(s=sum(value)) %>% as.data.frame()

# Redo labels
colz <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
colz$Sample <- ifelse(grepl("Water",colz$Sample),"Water Column",colz$Sample)

colz$Experiment <- ifelse(colz$Experiment=="Exp1","Experiment #1 - 2021",colz$Experiment)

colz$Experiment <- ifelse(colz$Experiment=="Exp2","Experiment #2 - 2021",colz$Experiment)

colz$Experiment <- ifelse(colz$Experiment=="Exp3","Experiment #3 - 2022",colz$Experiment)

dfMelt <- cbind(dfMelt,colz)
dfMelt$variable <- NULL
dfMelt <- dfMelt %>% group_by(KEGG,Sample,Experiment) %>% summarize(m=mean(s)) %>% as.data.frame()

dfMelt$Sample2 <- as.numeric(as.factor(dfMelt$Sample))
colrs <- distinctColorPalette(length(cath2)+1)

dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01275","Cathepsin C",NA)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01363","Cathepsin B",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01365","Cathepsin L",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01366","Cathepsin H",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01368","Cathepsin S",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01371","Cathepsin K",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01373","Cathepsin F",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01379","Cathepsin D",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K01382","Cathepsin E",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K08568","Cathepsin X",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K08569","Cathepsin W",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="K13289","Cathepsin A",dfMelt$KEGG2)
dfMelt$KEGG2 <- ifelse(dfMelt$KEGG=="Other Cathepsin","Other Cathepsin Protease",dfMelt$KEGG2)

# Plots (Rotting exp, water column)
r3 <- dfMelt %>% filter(Sample!="Water Column" & Experiment=="Experiment #3 - 2022") %>% ggplot(aes(x=Sample2, y=m, fill = KEGG2))+geom_area()+theme_bw()+scale_fill_manual(name="Cathepsin Protease",values=c(colrs))+scale_x_continuous(breaks=c(1,2,3),labels = c("T0","T3","T6"))+xlab("Incubation Time Point")+ylab("")+ylim(0,120000)


w3 <- dfMelt %>% filter(Sample=="Water Column" & Experiment=="Experiment #3 - 2022") %>% ggplot(aes(x=Sample2, y=m, fill = KEGG2))+geom_bar(stat="identity")+theme_bw()+scale_fill_manual(name="Cathepsin Protease",values=c(colrs))+xlab("Water Column")+ylab("Normalized Expression (TPM)")+scale_x_continuous(breaks=c(4),labels=c(""))+ylim(0,120000)

exp <- w1+r1+w2+r2+w3+r3+plot_layout(widths=c(5,20),guides = "collect",nrow=3)
exp
ggsave("Cathepsin_Plot.pdf",width=7,height=8)
