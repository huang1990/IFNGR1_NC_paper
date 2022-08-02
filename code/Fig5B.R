####KEGG plot 
a=read.csv('KEGG.pathway.txt',sep = '\t')
#colnames(a)=c("Category","Term","Count","Percentage","PValue","Genes","List Total","Pop Hits","Pop Total"	,"Fold Enrichment","Bonferroni","Benjamini","FDR")
a=subset(a, a$P.value<0.05)
a$log2p=-log2(a$P.value)
a <- a[!(a$Term=="Prostate cancer" | a$Term=="Breast cancer" | a$Term=="Colorectal cancer" | a$Term=="Melanoma" | a$Term=="Gastric cancer" | a$Term=="Basal cell carcinoma"), ]

ggplot(a, aes(x=reorder(Term,log2p),y=log2p))+ geom_bar(stat="identity",color="violetred4",fill="violetred4",width=0.8) + 
  coord_flip() + ylab("P value (-log2)") + xlab("Enriched Pathways")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",size=1),
        axis.title=element_text(size=14),axis.text.x=element_text(size=14,color="black"),
        axis.text.y=element_text(size=12,color="black"),
        plot.title = element_text(hjust = 0.5))


