

rm(list=ls())
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db",force = TRUE)
BiocManager::install("impute",force = TRUE)
BiocManager::install("preprocessCore",force = TRUE)
install.packages("WGCNA")
BiocManager::install("stats",force = TRUE)
library(WGCNA)
library(WGCNA)
library(data.table)
library(stringr)
BiocManager::install("openxlsx",force = TRUE)
library(openxlsx)

library(readxl)

#####
femData = read.csv("limma_DIFF_af (no).csv",stringsAsFactors = FALSE,check.names=F)
# Take a quick look at what is in the data set:
dim(femData)
View(femData[1:6,1:9])
datExpr0 = as.data.frame(t(femData[, -1]))


names(datExpr0) = femData$GeneName  
rownames(datExpr0) = names(femData)[-1]
View(datExpr0[1:5,1:5])
dataExpr <-datExpr0
dim(dataExpr)

#
gsg = goodSamplesGenes(dataExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
head(dataExpr)[,1:8]


#.Soft Threshold(Automatic or manual) 
powers <- c(c(1:15), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers)

# auto   if none  then picture it to chose
sft$powerEstimate #best beta

# Plot the results:
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h =0.8, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

# WGCNA  step by step
WGCNA_matrix<-dataExpr
adjacency <- adjacency(WGCNA_matrix, power = 20)#

TOM <-TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM),method = "average")


# module 
minModuleSize <-30 
dynamicMods<- cutreeDynamic(
  dendro= geneTree,
  distM= dissTOM,
  deepSplit= 4,
  pamStage= TRUE,
  pamRespectsDendro= FALSE,
  minClusterSize= minModuleSize)
 


# pictures

dynamicColors =labels2colors(dynamicMods)
tiff(file="wgcna_limma_dynamicColors_plotDendroAndColors.tiff",width=17,height=15,units="cm", compression="lzw", res=600)
p1<-plotDendroAndColors(
  geneTree,
  dynamicColors,
  "DynamicTree Cut",
  dendroLabels= FALSE, hang = 0.03,
  addGuide= TRUE,
  guideHang= 0.05,
  main= "Gene dendrogram and module colors")
print(p1)
dev.off()

#.between modules
MEList<- moduleEigengenes(WGCNA_matrix,colors = dynamicColors)
MEs<- MEList$eigengenes 
MEDiss<- 1-cor(MEs)
METree<- hclust(as.dist(MEDiss), method = "average")
tiff(file="wgcna_limma_METree.tiff",width=34,height=10,units="cm", compression="lzw", res=600)#宽度设置比较长，方便观察
p1<-plot(METree,
         main= "Clustering of module eigengenes",
         xlab= "", sub = "")
print(p1)
dev.off()

# Merge module
MEDissThres<- 0.1
merge<- mergeCloseModules(WGCNA_matrix,dynamicColors,cutHeight= MEDissThres,verbose = 3)
mergedColors <- merge$colors
table(mergedColors)
mergedMEs<- merge$newMEs 
tiff(file="wgcna_limma_mergedColors_plotDendroAndColors.tiff",width=17,height=17,units="cm", compression="lzw", res=600)
p1<-plotDendroAndColors(
  geneTree,
  cbind(dynamicColors,mergedColors),
  c("DynamicTree Cut", "Merged dynamic"),
  dendroLabels= FALSE,
  hang= 0.03,
  addGuide= TRUE,
  guideHang= 0.05)
print(p1)
dev.off()
table(mergedColors)

 
MEDiss<- 1-cor(mergedMEs) 
METree <- hclust(as.dist(MEDiss),method = "average")
tiff(file="wgcna_limma_METree.mergedColors.tiff",width=17,height=10,units="cm", compression="lzw", res=600)
p1<-plot(METree,
         main= "Clustering of module eigengenes",
         xlab= "", sub = "")
print(p1)
dev.off()


 
# modules heatmap
tiff(file="wgcna_limma_mergedColors.module_heatmap.tiff",width=17,height=17,units="cm", compression="lzw", res=600)
plotTOM<- dissTOM^7 
p1<-TOMplot(plotTOM, geneTree, mergedColors, main="Network heatmapplot, all genes")
print(p1)
dev.off()



tiff(file="wgcna_limma_mergedColors.Eigengeneadjacency heatmap.tiff",width=17,height=21,units="cm",compression="lzw", res=600)
p2<-plotEigengeneNetworks(mergedMEs,"Eigengene adjacencyheatmap", marHeatmap = c(2,2,2,2),plotDendrograms = T,xLabelsAngle = 90) 
print(p2)
dev.off()


# print out module gene
datME=moduleEigengenes(WGCNA_matrix,mergedColors)[[1]]

color1=as.character(mergedColors)
class(WGCNA_matrix)
class(datME)
datKME=signedKME(WGCNA_matrix[1], datME, outputColumnName="kME")
dataExpr1=as.data.frame(t(WGCNA_matrix));
datSummary=rownames(dataExpr1)
datout=data.frame(datSummary,colorNEW=color1,datKME )
write.table(datout, "wgcna_limma_all_gene_module.xls", sep="\t", row.names=F,quote=F)
dim(datout)



