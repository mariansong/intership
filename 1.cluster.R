

library(readxl)
#mainAll_time <- read_excel("mainAll_time.xls")
sample0  <- read_excel("mainAll_time.xls")
a<-t(sample0)
write.table(a,"sample.csv",sep=",") 

sample<- data.frame(read_excel("sample_lable.xls"))
#change 


#Aggregate function for average

dim(sample)
#[1]    8 256
sample1<-aggregate(sample[,2:256],by=list(sample$Lable),mean,na.rm= TRUE)

#Set new line name
row.names(sample1)<-sample1[,1]
sample1<-data.frame(t(sample1[,-1]))
sample1<-sample1[-1,]

#Next, start the analysis, install the Mfuzz package 

BiocManager::install("Mfuzz")
library(Mfuzz)

#Build the object
sample1<-as.matrix(sample1)
sample1<- ExpressionSet(assayData = sample1)
#Handling missing values and outliers
sample1 <- filter.NA(sample1, thres = 0.25)#Exclusion of genes with more than 25% of measured deletions
sample1 <- fill.NA(sample1, mode = 'mean')
sample1 <- filter.std(sample1, min.std = 0)
#Standardized processing
sample1 <- standardise(sample1)
#Set random seeds, set the number of clusters to be displayed, and then cluster
set.seed(123)
cluster_num <- 10
sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))

#plot picture
mfuzz.plot2(sample1, cl = sample1_cluster, mfrow = c(2, 5),
            time.labels = colnames(sample1),centre=TRUE,x11=F)


