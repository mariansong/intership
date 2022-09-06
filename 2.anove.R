
#ANOVA analysis

#######import the table 
library(readxl)
rt=read_excel("mainAll_time.xls")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp1=rt[,2:ncol(rt)]
exp=exp1[,2:ncol(exp1)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)



b<-data.frame(t(rt))
###
#row(b)
#row.names(b) <- c("","0a","0b","2a","2b","8a","8b","16a","16b")
write.csv(b,file = "anova_b1.csv")
##a litle change & write  b1_change
dat2<-data.frame(read.csv("anova_b1_change.csv"))

#
#Then, the time variable is set to the categorical variable.

dat2$time<-factor(dat2$time,levels = c("0","2","8","16"))
dat2$time



#Then.do ANOVA.

baseformula <- " ~ time"
#n m Counting, initialization
n<-0
m<-0
#b  Store the result of the for loop
b<-array(data = NA,dim=length(sample))
b1<-b
b2<-b
for (i in 2:ncol(dat2)) {
  formula <- paste(colnames(dat2)[i], baseformula, sep="")
  f <-summary(aov(as.formula(formula), data=dat2))[[1]][["F value"]][1]
  p <- summary(aov(as.formula(formula), data=dat2))[[1]][["Pr(>F)"]][1]
  if(f!=0){
    n=n+1
    print(paste(formula, ": f=", f, sep=""))
    b1[n]<-paste(formula, ": f=", f, sep="")
  }
  
  write.table(b1,file="anova_bf.txt")
  
  
  if(p!=0){
    m=m+1
    print(paste(formula, ": p=", p, sep=""))
    b2[m]<-paste(formula, ": p=", p, sep="")
  }
  
  write.table(b2,file="anova_bp.txt")
  
}
print(n)
print(m)
  


###ANOVAP&F anovaP&F.xls
ANOVAPf<-read_excel("anova_P&F.xls")
anova_data<-data.frame(ANOVAPf)
##Correcting the p-values with BH (fdr) method
#library(qvalue)
param <- names(anova_data)
p <- anova_data[,3]
p_fdr <- p.adjust(p, method = "fdr")
#p_fdr<-qvalue(anova2way_data$ANOVAP)
#Storing the corrected p-value back in the ANOVA object, in a new column
colN<-length(anova_data)+1
anova_data[,colN]<-p_fdr


names(anova_data)[colN]<-paste0(param,"FDR")

names(anova_data)[colN]<-"FDR"

#anova_data$GeneNamesFDR<0.02
##CHOSE
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

FDRcal<-data.frame(rt)
a1<-FDRcal

outall<-cbind(anova_data$GeneName,a1,anova_data$Fvalue,anova_data$Pvalue,anova_data$FDR)
colnames(outall) <- c("GeneName","0a","0b","2a","2b","8a","8b","16a","16b","ANOVAF","ANOVAP","FDR")
write.table(outall,file="anovaoutall.xls",quote=FALSE,sep="\t",row.names=FALSE)

re1 = outall %>% filter(FDR<0.02)
head(re1)
##630
write.table(re1,file="anova_diff_results.xls",quote=FALSE,sep="\t",row.names=FALSE)
dim(re1)



















# In this paper, we just need FDR,so ignore logFC

########can  ignore#######LOG FC 
##########################
#TRY calculate FC (FC for two groups)  
FDRcal<-data.frame(mainAll_time)
a1<-FDRcal


#Pre-generate 2 all-0 vectors of the same length as the number of lines in the input file, which will be used to store the p value and the difference multiplier (log2FC)
log2_FC<-c(rep(0,nrow(a1)))
FC<-c(rep(0,nrow(a1)))


for(i in 1:nrow(a1)){
  
  if(sum(a1[i,2:3])==0&&sum(a1[i,4:9])==0){
    
    log2_FC[i]<- "NA"
    FC[i] <- "NA"
    
  }else{
    
    log2_FC[i]<-log2((mean(as.numeric(a1[i,2:3]))+0.001)/(mean(as.numeric(a1[i,4:9]))+0.001))
    FC[i]<-(mean(as.numeric(a1[i,2:3]))+0.001)/(mean(as.numeric(a1[i,4:9]))+0.001)
    
  }
  
}


# Add log2FC, p value and FDR, in 3 columns, to the end of the original file.

out<-cbind(a1,log2_FC,FC)

write.table(out,file="anova.out.xls",quote=FALSE,sep="\t",row.names=FALSE)
anova_data
out
outall<-cbind(out,anova_data$ANOVA.F,anova_data$ANOVA.P,anova_data$GeneNameFDR)
colnames(outall) <- c("GeneName","0a","0b","2a","2b","8a","8b","16a","16b","log2_FC","FC","ANOVAF","ANOVAP","FDR")
write.table(outall,file="anovaoutall.xls",quote=FALSE,sep="\t",row.names=FALSE)






