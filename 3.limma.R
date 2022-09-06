if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

library("BiocManager")

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
BiocManager::install("WGCNA")
BiocManager::install("impute",force =TRUE)
BiocManager::install("IRanges",force =TRUE)
BiocManager::install("S4Vectors",force =TRUE)
BiocManager::install("stats",force =TRUE)

library("WGCNA")
BiocManager::install("GSEABase")
BiocManager::install("XML",force =TRUE)
library(GSEABase)
BiocManager::install("GSVA",force = TRUE)
library(GSVA)
library(readxl)

rt=read_excel("mainAll_time.xls")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp1=rt[,2:ncol(rt)]
exp=exp1[,2:ncol(exp1)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)



#Determine if the original data went to the log
max(rt)
if(max(rt)>30) rt=log2(rt+1)     #If the maximum value of rt is greater than 30, take log

#Use normalizeBetweenArrays for correction, and assign the value of rt1 after correction
rt1=normalizeBetweenArrays(as.matrix(rt))

#Not standardized
cols=rainbow(ncol(rt)) 
par(cex = 0.7)
if(ncol(rt)>40) par(cex = 0.5)  

boxplot(rt,las=2,col =cols ) 
#dev.off()

# standardized
cols=rainbow(ncol(rt1)) 
par(cex = 0.5)
if(ncol(rt1)>40) par(cex = 0.5)  
pdf(file = "limma_nor.pdf",width=5,height = 4.5)
boxplot(rt1,las=2,col =cols )
dev.off()

#save 
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file="limma_norexp.txt",sep="\t",quote=F,col.names = F)





####DE 
data=rt1

h0Data=data[,as.vector(colnames(data)[1:2])]
h2Data=data[,as.vector(colnames(data)[3:4])]
h8Data=data[,as.vector(colnames(data)[5:6])]
h16Data=data[,as.vector(colnames(data)[7:8])]
rt=cbind(h0Data,h2Data,h8Data,h16Data)

h0Num=ncol(h0Data)
h2Num=ncol(h2Data)
h8Num=ncol(h8Data)
h16Num=ncol(h16Data)

#limma all 
Type=c(rep("h0",h0Num),rep("h2",h2Num),rep("h8",h8Num),rep("h16",h16Num))
#define
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("h0","h2","h8","h16")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(h0-h2,h0-h8,h0-h16,levels=design)
help(make.names)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
head(fit2)
Diff1=topTable(fit2,adjust ="fdr",number=length(rownames(data)))
Diff<-cbind(Diff1,rt1)



#Save differential results for all genes
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="limma_DIFF_all.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (adj.P.Val < 0.02 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="limma_DIFF_af.xls",sep="\t",quote=F,col.names=F)
dim(diffSigOut)





#Heat map showing the top 30 most diverse genes
Diff=Diff[order(as.numeric(as.vector(Diff$adj.P.Val))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(100)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]

#####
install.packages("pheatmap")
library(pheatmap)

pheatmap(afExp)













#Grouping Tags
Type=c(rep("h0",h0Num),rep("h2",h2Num),rep("h8",h8Num),rep("h16",h16Num))

names(Type)=colnames(rt)
Type=as.data.frame(Type)
#Annotation colors for grouped labels
ann_colors=list(gene_class=c(h0='#CC6666',h2='#3366FF',h8='#FDDCA9',h16="#FF6600"))
pdf(file="limma_DIFF_heatmap1.pdf",height=7,width=10)
pheatmap(afExp,                                                                     
         annotation=Type,                                                           
         color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     
         #not sure 
    
         cluster_cols =F,                                                           
         show_colnames = F,                                                         
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=ann_colors
)
dev.off()


