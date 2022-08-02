#setwd("/Users/emmahuang/Downloads/Chong's_lab/Lewis/TCGA_UVM_SKCM/final_figure/selectFigure")
library(ggplot2)
library(reshape2)
library(plyr)
df <- read.table("GSE78220_60_Myc_STAT3_exp.txt", sep = "\t", header = T)
df$expression <- log10(df$expression)
source("./Function_for_violin_plot.R")
Data_summary <- summarySE(df, measurevar="expression", groupvars=c("type","symbol"))

df$symbol<- factor(df$symbol,
                         levels = c("ENO1", "FASN","FKBP4", "ODC1", "GADD45A", "JUNB", "VEGFA"),ordered = TRUE)
df <- df %>%
  arrange(expression) %>%
  mutate(type=factor(type, levels = c("res", "non")))

VEGFA <- df[df$symbol=="VEGFA",]
non<-VEGFA[VEGFA$type=="non",]
res <- VEGFA[VEGFA$type=="res",]
wilcox.test(non$expression,res$expression)
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

gene_split_violin <- ggplot(df,aes(x= symbol,y= expression,fill= type))+
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
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(fill="",y=("Expression (log10)"),x=NULL,title = "GSE78220") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = type),
                     label = "p.signif",
                     method = "wilcox.test",
                     label.y =4,
                     size=7,
                     hide.ns = T)
gene_split_violin

ggsave(gene_split_violin,
                         filename = "./Output/gene_split_violin.pdf",
                         height = 10,width = 16,units = "cm")
