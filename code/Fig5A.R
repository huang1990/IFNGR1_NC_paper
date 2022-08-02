library("apeglm")
library( "DESeq2" )
library(ggplot2)
library(tidyr)
library(dplyr)
countData <- read.csv("B6merge.txt", sep = "\t", header=T, row.names=1)
head(countData)
countData=apply(countData,c(1,2),as.integer)
countData[is.na(countData)] <- 0
###GROUP1 VS GROOUP2
countData <- countData[,c(4,5,6,1,2,3)]
head(countData)
###input prepare
metaData =read.table("samplelist.txt", sep = '\t', header = T)
head(metaData)
metaData <- metaData[c(4,5,6,1,2,3),]
head(metaData)
#metaData<-metaData[order(metaData$condition),]
rownames(metaData)=metaData$sampleID
###check all rowname in countdata
all(rownames(metaData) %in% colnames(countData))
metaData$condition <- factor(metaData$condition)
#metaData=metaData[,c('condition','response','patient','combined')]

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~condition)
####pre-filter rows have at least 11 reads total.
keep <- rowSums(counts(dds)) >= 6
dds <- dds[keep,]
##on factor level

dds$condition <- factor(dds$condition, levels = c("IFNGR1KOclone5", "scr"))## the second one is ref level 
##differential expression analysis
dds <- DESeq(dds)
####
res <- results(dds)
res <- results(dds, contrast=c("condition","IFNGR1KOclone5", "scr")) ###the second one if ref level
res
#plotMA(resOrdered)
resOrdered <- res[order(res$pvalue),]
###save expression
ntd <- normTransform(dds)
filteredAssay=(assay(ntd))[rownames((assay(ntd)))%in%rownames(resOrdered),]
#write.csv(as.data.frame(filteredAssay),
 #         file="./scrVSKO5/scrVSKO5_all_normalized_expression.csv")

##organize dataframe
rownames(resOrdered) <- sub("\\..*","\\", rownames(resOrdered))
resOrdered$eensembleID=rownames(resOrdered) 
rownames(resOrdered) <- NULL
library(org.Mm.eg.db)
library(tidyr)
ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
                                    key=resOrdered$eensembleID,
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
out <- merge(as.data.frame(resOrdered), ens2symbol, by.x="eensembleID",by.y="ENSEMBL")
out<-out[!duplicated(out$SYMBOL),]
library(tidyr)
out<-out %>% drop_na()
###volcano plot 
out$Diffexpressed<-"Not Significant"
out$Diffexpressed[out$log2FoldChange > 1 & out$pvalue < 0.05] <- "Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
out$Diffexpressed[out$log2FoldChange < -1 & out$pvalue < 0.05] <- "Down"
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
head(out)
out$delabel <- NA

##label name 
#delabel<- c("EGF","LAMA4","COL1A1","EFNA1","COL1A2","COL2A1","ERBB3","DDIT4","EIF4EBP1","THEM4","ITGA7","GNG8","ATF4",
#"PCK2","FGF21","FZD4","SESN2","DDIT4","EIF4EBP1","SLC3A2","WNT9A")
legend_title<-""
#out$delabel[out$diffexpressed != "NO"] <- out$SYMBOL[out$diffexpressed != "NO"]
library(ggrepel)
ggplot(data=out, aes(x=log2FoldChange, y=-log10(pvalue), col=Diffexpressed)) + 
  geom_point() +
  theme_minimal() +
 # geom_text_repel(
 # data = out[tolower(out$SYMBOL)%in%tolower(delabel),], #subset(out, -log10(out$padj)>15)
 #   aes(label = SYMBOL), size=5, color="red4", box.padding = unit(0.35, "lines"), 
#    point.padding = unit(0.4, "lines")) + 
  #geom_text(aes(label=ifelse(log10(padj)>15,as.character(SYMBOL),'')),hjust=0,vjust=0 )
  ylim(0,40) + 
  xlab("Log2FoldChange") +  
  ylab("P value (-log10)") +
  scale_color_manual(legend_title, values = c("blue", "black","red")) + 
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text = element_text(size=14, color = "black"), axis.title = element_text(size=14), 
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        legend.position = c(0.15,0.90)) 
  





