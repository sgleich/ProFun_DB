### PARAGON ANALYSES - TAXA BARPLOT - FUNGI ###
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
dfWide <- subset(dfWide,grepl("Fungi",dfWide$Taxonomy)|grepl("Ascomycota",dfWide$Taxonomy))

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

##### TAXA BARPLOT #####
dfTax <- dfFin[c(2,4:18)]

taxz <- colsplit(dfTax$Taxonomy,";",c("k","p","c","o","f","g","s"))
taxz2 <- NULL
for(i in 1:nrow(taxz)){
  x <- taxz[i,]
  if(x$c==" Ascomycota"){
    x$s <- x$g
    x$g <- x$f
    x$f <- x$o
    x$o <- x$c
    x$c <- " Fungi"
  }
  taxz2 <- rbind(taxz2,x)
}

taxz2$fin <- ifelse(taxz2$f==" Sordariomycetes","Sordariomycetes",NA)
taxz2$fin <- ifelse(taxz2$f==" Leotiomycetes","Leotiomycetes",taxz2$fin)
taxz2$fin <- ifelse(taxz2$f==" Saccharomycotina","Saccharomycotina",taxz2$fin)
taxz2$fin <- ifelse(is.na(taxz2$fin) & taxz2$o==" Ascomycota","Unknown Ascomycota",taxz2$fin)
taxz2$fin <- ifelse(is.na(taxz2$fin) & taxz2$o=="","Unknown Fungi",taxz2$fin)
unique(taxz2$fin)

dfTax$tax <- taxz2$fin
dfTax$Taxonomy <- NULL
dfTaxMelt <- melt(dfTax,id.vars="tax")
dfTaxMelt <- dfTaxMelt %>% group_by(variable,tax)%>% summarize(s=sum(value))

varSplit <- colsplit(dfTaxMelt$variable,"_",c("Sample","Exp"))
dfTaxMelt <- cbind(dfTaxMelt,varSplit)

colrs <- c("Unknown Ascomycota"="#DA9774","Leotiomycetes"="#AEB2D7","Unknown Fungi"="#ACDF68","Saccharomycotina"="#C164CD", "Sordariomycetes"="#A8DCBF")

taxPltFxn <- function(df,exp,title){
  subs <- subset(df, Exp==paste(exp))
  subsWater <- subset(subs,grepl("Water",subs$Sample))
  subsWater <- subsWater %>% group_by(tax) %>% summarize(s=mean(s))
  subsWater$Sample <- "Water Column"
  subs <- subset(subs,!grepl("Water",subs$Sample))
  subs$variable <- NULL
  subs$Exp <- NULL
  subs$Sample <- ifelse(subs$Sample=="RotT0",0,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT3",3,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT6",6,subs$Sample)
  
  
  waterPlot <- ggplot(subsWater,aes(x=Sample,y=s,fill=tax))+geom_bar(stat="identity",color="grey20",position="fill")+scale_fill_manual(values=c(colrs))+theme_classic(base_size=16)+xlab("")+ylab("Relative Transcript Abundance")+theme(axis.text.x = element_text(angle = 45, hjust =1))+theme(legend.position = "none")
  
  rotPlot <- ggplot(subs,aes(x=as.numeric(Sample),y=s,fill=tax))+geom_area(alpha = 0.6, position = 'fill')+geom_col(width = 1.5, color = 'gray20', position = 'fill')+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic(base_size = 16)+xlab("")+ylab("Relative Transcript Abundance")+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_x_continuous(breaks=c(0,3,6),labels=c("Net Trap\nDay 0","Net Trap\nDay 3","Net Trap\nDay 6"))+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  out <- waterPlot+rotPlot+plot_layout(widths=c(5,20),guides = "collect",nrow=1)+plot_annotation(title=paste(title))
  return(out)}

pExp1 <- taxPltFxn(dfTaxMelt,"Exp1","Experiment #1\n2021")
pExp2 <- taxPltFxn(dfTaxMelt,"Exp2","Experiment #2\n2021")
pExp3 <- taxPltFxn(dfTaxMelt,"Exp3","Experiment #3\n2022")
# ggarrange(pExp1,pExp2,pExp3,common.legend = TRUE,nrow=1,ncol=3,labels=c("a","b","c"),font.label = list(size=12,color="black",face="plain"),hjust=-0.1,vjust=0.5)
ggarrange(pExp1,pExp2,pExp3,common.legend = TRUE,nrow=1,ncol=3)

ggsave("TaxaBarplot_Fungi.pdf",width=13,height=6)
