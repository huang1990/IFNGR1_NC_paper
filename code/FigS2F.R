library(survival)
library(survminer)
library(Cairo)
library(dplyr)
clinical=read.table("SKCM.clinical.txt", sep = "\t", header=T) 
expression=read.table("SKCM.IFNGR1.txt",sep = "\t", header =T)

head(clinical)
head(expression)
#rownames(ids)=ids$external_gene_name
rownames(expression)=expression$gene
Name=c()
Pvalue=c()
Symbol=c()

symbol="IFNGR1|3459"
name="IFNGR1|3459"

#gene=data.frame("expression"=t(expression[name,])[-1,]) 
gene=data.frame(t(expression[name,]))
gene$number=as.numeric(as.character(gene[,1]))

qgene=subset(gene, !is.na(gene$number))
gene$submitter_id= substr(row.names(gene),1,12)
gene=gene[substr(rownames(gene),14,14)=="0",]
gene$rank= ifelse(gene$number> median(gene$number),"High","Low")
clinical$submitter_id=gsub("-",'.',clinical$submitter_id)

gene_clinical=merge(clinical,gene,by.x="submitter_id",by.y="submitter_id")
gene_clinical$vital_status

gene_clinical$status[gene_clinical$vital_status=="Alive"]=0
gene_clinical$status[gene_clinical$vital_status=="Dead"]=1
gene_clinical$days=ifelse(gene_clinical$vital_status=="Alive",gene_clinical$days_to_last_follow_up,gene_clinical$days_to_death)
gene_clinical=gene_clinical[!is.na(gene_clinical$days),]
gene_clinical[,c('days','status','rank')]
sfit= survfit(Surv(days,status)~rank,data=gene_clinical)
Name=c(Name,name)
Symbol=c(Symbol,symbol)
Pvalue=c(Pvalue,surv_pvalue(sfit)$pval)
if(surv_pvalue(sfit)$pval<=0.5){
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(
        plot.title=element_text(hjust=0.5)
      )
  }
  ###SKCM set xlim(0,4000)
  ggsurvplot(sfit,pval=TRUE,
             size=1,
             palette = c("blue", "green"),
             xlim=c(0,4000),break.x.by = 1000, 
             title = "SKCM",ggtheme=custom_theme(),
             legend.title = "IFNGR1",
             legend.labs = c("High", "Low"),
             pval.coord = c(3200, 0.85),
             pval.size = 6,
             xlab="Time / Days",
             font.x = c(20,face="bold"),
             font.y = c(20, face = "bold"),
             font.tickslab = c(18),
             font.legend = c(18),
             font.title=c(22, "bold"))
}


