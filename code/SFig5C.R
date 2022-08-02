library(dplyr)
library(tidyr)
##SKCM high vs low 
df=read.csv("SKCM80_Myc_STAT3_exp.txt",sep = "\t")

df <- subset(df, select = -c(CCND2.894, EIF4E.1977, HSPD1.3329, PPAT.5471) )
##mediun if IFNGR1
median(df$IFNGR1.3459)
# creating a column to dataframe based on values in other columns:
df <- df %>% 
  mutate(type = if_else(IFNGR1.3459 >= 900.816, "High", "Low"))

##transfer data matrix
df_new <- tidyr::gather(df,key='symbol',value='expression',2:9,-1) 
df_new$symbol <- gsub("\\.[0-9]+$","",df_new$symbol)
df_new <- df_new[df_new$symbol!="IFNGR1",]

##plot 
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#df_new$symbol <- gsub("\\.[0-9]+$","",df_new$symbol)
df_new$expression <- log10(df_new$expression)
source("./Function_for_violin_plot.R")
Data_summary <- summarySE(df_new, measurevar="expression", groupvars=c("type","symbol"))
head(Data_summary)

df_new$symbol<- factor(df_new$symbol,
                   levels = c("ENO1", "FASN","FKBP4", "ODC1", "GADD45A", "JUNB", "VEGFA"),ordered = TRUE)

##plot
if(T){
  mytheme <- theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   axis.line = element_line(color = 'black'),
                   plot.title = element_text(size = 24,color="black",hjust = 0.5),
                   axis.title = element_text(size = 22,color ="black"),
                   axis.text = element_text(size= 18,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(size=18, color = "black"),
                   panel.grid=element_blank(),
                   legend.position = c(0.9,0.85),
                   legend.text = element_text(size= 18),
                   legend.title= element_text(size= 18)
  )
}

gene_split_violin <- ggplot(df_new,aes(x= symbol,y= expression,fill= type))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x= symbol, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci),
                width= 0.05,
                position= position_dodge(0.5),
                color="black",
                alpha = 0.8,
                size= 0.5) +
  ylim(-1,7) +
  scale_fill_manual(
    labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)),
                    values = c("blue", "green"))+
  labs(fill="",y=("Expression (log10)"),x=NULL,title = "TCGA_SKCM_deconvolution") +
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = type),
                     label = "p.signif",
                     method = "wilcox.test",
                     label.y =6,
                     size=7,
                     hide.ns = T)
gene_split_violin


ggsave(gene_split_violin,
       filename = "./Output/gene_split_violin.pdf",
       height = 10,width = 16,units = "cm")
