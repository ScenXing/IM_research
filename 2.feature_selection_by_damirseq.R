library(DESeq2) # 1.12.4
library(glmnet) # 2.0-5
library(pROC)   # ROC curve
library(pheatmap) # heatmap
library(DaMiRseq) # feature selections
library("biomaRt") # Retrive information for Ensembl
library(hash) # Perl hash
library("reshape2") # melt function
library(ggplot2) # plotting

###########################################
############------DESeq2------#############
###########################################
#---Count data---#
colnames(countData)<-grep('CSF',unlist(strsplit(colnames(countData),'_')),value=T)
#---Annotation information---#
colData$class <-factor(colData$class,levels = c(0,1,2))
colData<-colData[,c('batch','gender','age','collection_days','class')]
colData$gender<-as.factor(colData$gender)
colData$batch<-as.factor(colData$batch)
colData<-colData[intersect(colnames(countData),row.names(colData)),]
countData<-countData[,row.names(colData)]
#-----feature selections-----#
SE<-DaMiR.makeSE(countData,colData) #assay(SE), colData(SE)
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,hyper = "yes", th.cv=3)
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.2)
sv <- DaMiR.SV(data_filt)
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=4)

###########################################
#####--------Feature selection---------####
###########################################

data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.3)
data_reduced <- DaMiR.FReduct(data_reduced$data)
DaMiR.MDSplot(data_reduced, df,type="pearson")

df.importance <- DaMiR.FSort(data_reduced, df)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,n.pred = 53)
DaMiR.Clustplot(selected_features$data, df)
