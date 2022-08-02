df=read.csv("SKCM_UVM08_Immune_expression.txt",sep = "\t")
head(df)
##transfer data matrix
df_new <- tidyr::gather(df,key='symbol',value='expression',2:13,-1)
head(df_new)
#df_new$expression <- log10(df_new$expression)
df_new$symbol <- gsub("\\.[0-9]+$","",df_new$symbol)
df_new$symbol<-gsub("\\.", "-",df_new$symbol)
library(tidyr)
head(df_new)
res <- df_new %>% 
  dplyr::filter(!is.infinite(expression)) %>%
  group_by(type,symbol) %>%
  dplyr::summarise(value = list(expression)) %>%
  spread(type,value) %>%
  group_by(symbol) %>%
  mutate(high_mean=mean(unlist(SKCM))) %>%
  mutate(low_mean=mean(unlist(UVM))) %>%
  mutate(p_value=wilcox.test(unlist(SKCM),unlist(UVM))$p.value)

new <- res[,c(1,4,5,6)]

#write.csv(new, "UVM_SKCM_STAT3Target.csv")


##plot 
library(ggsignif)
library(ggplot2)
library(tidyverse)
library(ggpubr)
#df_new$symbol <- gsub("\\.[0-9]+$","",df_new$symbol)
df_new$expression <- log10(df_new$expression)
UVM_SKCM <- df_new[,-1]

UVM_SKCM$symbol<- factor(UVM_SKCM$symbol, 
                         levels = c("IFNGR1", "CD3D","CD4","CD8A", "HLA-A", "HLA-B", 
                                    "HLA-C", "HLA-DRA", "GZMB", "IFNG", "PRF1", "TNF"),ordered = TRUE)

p1<-ggplot(data = UVM_SKCM, aes(x=symbol, y=expression)) + 
  geom_boxplot(aes(fill=type), outlier.size = 0.1) +
  ylim(-1,7) +
  labs(fill = "") +
  ylab("Expression (log10)") +
  xlab("Gene")+
  theme_bw() +
  #eliminates background, gridlines, and chart border
  scale_fill_manual(values=c("coral3", "cyan3"))+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black')) +
  theme(axis.text.y = element_text(size=15, color = "black", family = "Arial" ),
        axis.text.x = element_text(size=15, color = "black", family = "Arial"), 
        axis.title = element_text(size=20, family = "Arial"), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 27, hjust = 0.5, family = "Arial"),
        legend.text = element_text(size=14, family = "Arial"),
        
        #legend.title = element_text(size=16),
        legend.position = c(0.92,0.75))   ## "coral3", "cyan3" #c(0.92,0.8)
p1 + stat_compare_means(aes(group = type),method="t.test", label = "p.signif",
                        label.y = 6.5, size =6, hide.ns = TRUE) + ggtitle("SKCM vs UVM (after deconvolution)")

