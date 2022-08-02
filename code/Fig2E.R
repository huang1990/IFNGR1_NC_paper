library(dplyr)
library(tidyr)
library(ggplot2)
df=read.csv("SKCM80.immuneGene.txt",sep = "\t")
head(df)
##mediun if IFNGR1
mean(df$IFNGR1.3459)
# creating a column to dataframe based on values in other columns:
df <- df %>% 
  mutate(type = if_else(IFNGR1.3459 >= 888.9253, "High", "Low"))#80SKCM:795.157,UVM:376.01#90SKCM:677.72 UVM:400.34  #totalSKCM 900.816, UVM: 451.418

##transfer data matrix
df_new <- tidyr::gather(df,key='symbol',value='expression',2:13,-1) ##
head(df_new)
#df_new$expression <- log10(df_new$expression)
df_new$symbol <- gsub("\\.[0-9]+$","",df_new$symbol)
df_new$symbol<-gsub("\\.", "-",df_new$symbol)
head(df_new)
res <- df_new %>% 
  dplyr::filter(!is.infinite(expression)) %>%
  group_by(type,symbol) %>%
  dplyr::summarise(value = list(expression)) %>%
  spread(type,value) %>%
  group_by(symbol) %>%
  mutate(high_mean=mean(unlist(High))) %>%
  mutate(low_mean=mean(unlist(Low))) %>%
  mutate(p_value=wilcox.test(unlist(High),unlist(Low))$p.value)

new <- res[,c(1,4,5,6)]
new

##plot 
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)

df_new$expression <- log10(df_new$expression)
UVM_SKCM <- df_new[,-1]

##reorder x axis
UVM_SKCM$symbol<- factor(df_new$symbol,
                        levels = c("IFNGR1", "CD3D","CD4","CD8A", "HLA-A", "HLA-B", 
                "HLA-C", "HLA-DRA", "GZMB", "IFNG", "PRF1", "TNF"),ordered = TRUE)


p1<-ggplot(data = UVM_SKCM, aes(x=symbol, y=expression)) + 
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
  ylim(-1,7) +
  labs(fill = "IFNGR1") +
  ylab("Expression (log10)") +
  xlab("Gene")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(
    labels = c(expression(IFNGR1^High),expression(IFNGR1^Low)), 
    values = c("blue", "green")) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=15, color = "black"),
        #    axis.text.x = element_text(size=10, color = "black", angle = 45, hjust=1), axis.title = element_text(size=20, face = "bold"), 
        axis.text.x = element_text(size=15, color = "black"), axis.title = element_text(size=22), 
        axis.title.x = element_blank(),plot.title = element_text(size = 27, hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        #legend.title = element_text(size=13),
        legend.position = c(0.9,0.7)) #c(0.9,0.8)
p1 + stat_compare_means(aes(group = type),method="wilcox.test", label = "p.signif",
                        label.y = 6.5, size =6, hide.ns = TRUE) + ggtitle("TCGA_SKCM_deconvolution")

