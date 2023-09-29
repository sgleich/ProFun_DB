### PARAGON ANALYSES - CAZy Upset Plot ###
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
full$KEGG_ko <- NULL
full$CAZy <- ifelse(is.na(full$CAZy),"",full$CAZy)
full$CAZy <- str_split(full$CAZy,",")
dfWide <- unnest(full,CAZy)
colnames(dfWide)[1:3] <- c("Name","Taxonomy","CAZy")
dfWide <- subset(dfWide,grepl("Eukaryot",dfWide$Taxonomy))

# Separate Exp #1 and Exp #2 from Exp #3 (Experiment #3 was sequenced separately)
dfExp12 <- dfWide[,c("Name","Taxonomy","CAZy","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3 <- dfWide[,c("Name","Taxonomy","CAZy","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 
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

# Upset plot
dfFin$Name <- NULL
dfFin$Taxonomy <- NULL
dfFin <- dfFin %>% group_by(CAZy) %>% summarize_all(sum)%>%as.data.frame()

dfMelt <- melt(dfFin,id.vars="CAZy")
namez <- colsplit(dfMelt$variable,"_",c("Sample","Exp"))
namez$Sample <- ifelse(grepl("Water",namez$Sample),"Water Column",namez$Sample)
dfMelt <- cbind(dfMelt,namez)
dfMelt$variable <- NULL
dfMelt <- dfMelt %>% group_by(Exp,Sample,CAZy) %>% summarize(m=mean(value))
dfMelt <- subset(dfMelt,CAZy!="")
dfMelt$SampleExp <- paste(dfMelt$Sample,dfMelt$Exp,sep="_")
dfMelt$Exp <- NULL
dfMelt$Sample <- NULL

dfMelt <- dfMelt %>% pivot_wider(id_cols="CAZy",names_from="SampleExp",values_from="m")%>% as.data.frame()
rownames(dfMelt) <- dfMelt$CAZy
dfMelt$CAZy <- NULL

dfMeltNew <- select(dfMelt,contains("Exp3"))

colnames(dfMeltNew)
colnames(dfMeltNew) <- c("Net Trap T0","Net Trap T3","Net Trap T6","Water Column")

dfMeltBin <- ifelse(dfMeltNew > 0, 1,0)
dfMeltBin <- as.data.frame(dfMeltBin)
dfMeltBin$Intersect <- apply(dfMeltBin > 0, 1, function(x){toString(names(dfMeltBin)[x])})
dfMeltBin$Intersect <- stri_replace_all_fixed(dfMeltBin$Intersect, " ", "")
dfMeltBin$Intersect <- as.list(strsplit(dfMeltBin$Intersect, ","))

keep <- list(c("NetTrapT0"),c("NetTrapT3"),c("NetTrapT6"),c("WaterColumn"),c("NetTrapT3","NetTrapT6"),c("NetTrapT0","NetTrapT3"),c("NetTrapT0","NetTrapT3","NetTrapT6"),c("NetTrapT0","NetTrapT3","NetTrapT6","WaterColumn"))

keep_df <- subset(dfMeltBin,Intersect %in% keep)
keep_df$C <- ifelse(grepl("AA",rownames(keep_df)),"AA",NA)
keep_df$C <- ifelse(grepl("CBM",rownames(keep_df)),"CBM",keep_df$C)
keep_df$C <- ifelse(grepl("GH",rownames(keep_df)),"GH",keep_df$C)
keep_df$C <- ifelse(grepl("GT",rownames(keep_df)),"GT",keep_df$C)
keep_df$C <- ifelse(grepl("CE",rownames(keep_df)),"CE",keep_df$C)
keep_df$C <- ifelse(grepl("PL",rownames(keep_df)),"PL",keep_df$C)

colrs <- c("GH"="dodgerblue",
           "GT"="indianred",
           "CE"="darkgoldenrod3",
           "CBM"="forestgreen",
           "AA"="salmon","PL"="purple")

exp1 <- keep_df%>% ggplot(aes(x=Intersect,fill=C)) + geom_bar(stat = "count", position="stack") + scale_x_upset(sets=c("WaterColumn","NetTrapT0","NetTrapT3","NetTrapT6"))+theme_classic(base_size = 12)+xlab("")+ylab("Number of Shared/Unique CAZymes")+theme(plot.margin = margin(10, 10, 10, 100))+ggtitle("Experiment #1\n2021")+scale_fill_manual(name=c("CAZyme Class"),values=c(colrs))+ylim(0,85)
exp1

exp1+exp2+exp3+plot_layout(guides="collect",nrow=1)#+plot_annotation(tag_levels="a")
ggsave("CAZy_Upset.pdf",width=18,height=6)
