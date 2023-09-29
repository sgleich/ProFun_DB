### PARAGON ANALYSES - KEGG Heatmap ###
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

# Wrangle data
dfFin$Taxonomy <- NULL
dfFin$Name <- NULL
dfFin <- dfFin %>% group_by(KEGG) %>% summarize_all(sum)%>% as.data.frame()
dfFin <- dfFin[rowSums(dfFin[2:16] == 0)==0, ]

dfMelt <- melt(dfFin,id.vars="KEGG")
dfMelt <- subset(dfMelt,KEGG!="")
dfMelt <- subset(dfMelt,!is.na(KEGG))

colz <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
colz$Sample <- ifelse(grepl("Water",colz$Sample),"Water Column",colz$Sample)
colz$Sample <- ifelse(grepl("RotT0",colz$Sample),"Net Trap T0",colz$Sample)
colz$Sample <- ifelse(grepl("RotT3",colz$Sample),"Net Trap T3",colz$Sample)
colz$Sample <- ifelse(grepl("RotT6",colz$Sample),"Net Trap T6",colz$Sample)

dfMelt <- cbind(dfMelt,colz)

dfMeans <- dfMelt %>% group_by(KEGG,Experiment)%>%summarize(mz=mean(value),sdz=sd(value))

dfMelt <- left_join(dfMelt,dfMeans)
dfMelt$z <- (dfMelt$value-dfMelt$mz)/dfMelt$sdz

dfMelt$variable <- NULL
dfMelt$value <- NULL
dfMelt$mz <- NULL
dfMelt$sdz <- NULL

# Summarize Zscore
out <- dfMelt %>% group_by(Experiment,Sample,KEGG)%>%summarize(m=mean(z))
out <- subset(out,!is.na(m))

out$Sample <- factor(out$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

outF<- NULL
koz <- unique(out$KEGG)
for (i in 1:length(koz)){
  t <- subset(out,KEGG==koz[i])
  exp1 <- subset(t,Experiment=="Exp1")
  exp2 <- subset(t,Experiment=="Exp2")
  exp3 <- subset(t,Experiment=="Exp3")
  row1 <- exp1[(which.max(exp1$m)),2]
  row2 <- exp2[(which.max(exp2$m)),2]
  row3 <- exp3[(which.max(exp3$m)),2]
  if (row1$Sample==row2$Sample){
    if(row2$Sample==row3$Sample){
      sampl <- data.frame(sampl=row2$Sample,kegg=koz[i])
      
      outF <- rbind(outF,sampl)
      
    }
  }
  
}

outFinWater <- subset(outF, sampl=="Water Column")
outFinRot0<- subset(outF, sampl=="Net Trap T0")
outFinRot3 <- subset(outF, sampl=="Net Trap T3")
outFinRot6 <- subset(outF, sampl=="Net Trap T6")

outFinWater <- subset(out,KEGG %in% outFinWater$kegg)
outFinRot0 <- subset(out,KEGG %in% outFinRot0$kegg)
outFinRot3 <- subset(out,KEGG %in% outFinRot3$kegg)
outFinRot6 <- subset(out,KEGG %in% outFinRot6$kegg)
outFinWater$color <- "black"
outFinRot0$color <- "indianred"
outFinRot3$color <- "dodgerblue"
outFinRot6$color <- "forestgreen"

outFin <- rbind(outFinWater,outFinRot0,outFinRot3,outFinRot6)
outFin$color <- factor(outFin$color,levels=c("forestgreen","dodgerblue","indianred","black"))

outFin$Experiment <- ifelse(outFin$Experiment=="Exp1","Experiment #1\n2021",outFin$Experiment)
outFin$Experiment <- ifelse(outFin$Experiment=="Exp2","Experiment #2\n2021",outFin$Experiment)
outFin$Experiment <- ifelse(outFin$Experiment=="Exp3","Experiment #3\n2022",outFin$Experiment)

outFin$High2 <- as.numeric(outFin$color)

colrs <- c(rep("black",length(unique(outFinWater$KEGG))),rep("indianred",length(unique(outFinRot0$KEGG))),rep("dodgerblue",length(unique(outFinRot3$KEGG))),rep("forestgreen",length(unique(outFinRot6$KEGG))))


ggplot(outFin,aes(x=Sample,y=reorder(KEGG,High2),fill=m))+geom_tile(color="grey")+facet_grid(~Experiment)+theme_classic(base_size=16)+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_fill_gradient2(high="red",low="blue",name="Z-score transformed\ntranscript abundance")+ggtitle("")+ylab("KEGG KO")+theme(axis.text.y = element_markdown(colour = rev(colrs))) 
ggsave("KEGG_HeatMap.pdf",width=8,height=6)
