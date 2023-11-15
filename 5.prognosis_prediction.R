library(DESeq2) ### 1.12.4
library(glmnet) ### 2.0-5
library(pROC)   ### ROC curve
library(pheatmap) ### heatmap
library(DaMiRseq) ### feature selections

###########################################
###--------------Obtain DEGs------------###
###########################################
library("AnnotationDbi")
library("org.Hs.eg.db")
ensg2sym <- mapIds(org.Hs.eg.db,
                     keys=row.names(countData),
                     column=c("ENSEMBL","SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals="first")
#-----feature selections-----#
SE<-DaMiR.makeSE(countData,colData) #assay(SE), colData(SE)
#-3.3.3 normalization
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,hyper = "yes", th.cv=3) 
#-3.3.4 sample filtering
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.2)
sv <- DaMiR.SV(data_filt)
data_adjust<-DaMiR.SVadjust(data_filt,sv,n.sv=4)

#-----feature selection-----#
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.2) 
data_reduced <- DaMiR.FReduct(data_reduced$data) 
DaMiR.MDSplot(data_reduced, df,type="pearson")

df.importance <- DaMiR.FSort(data_reduced, df)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,n.pred = 10)
DaMiR.Clustplot(selected_features$data, df) #heatmap

################################################
#####-----Prediction model for Prognosis----####
################################################
ydat<-as.matrix(data_reduced)
cvfit <- cv.glmnet(x = ydat[-idx_test,], y = colData$class[-idx_test],family = "binomial", 
                   type.multinomial="grouped",
                   type.measure = "class",nfolds =51-length(idx_test))
tscores<-cbind(colData[-idx_test,],
               predict(cvfit, newx=ydat[-idx_test,], s="lambda.1se", type="class"),
               predict(cvfit, newx=ydat[-idx_test,], s="lambda.1se", type="response")
)
colnames(tscores)[7:8]<-c('pre_class','pre_score')
tscores$class<-as.numeric(as.character(tscores$class))

#################################################
#####------------ROC of the model------------####
#################################################
gscores<-cbind(colData[idx_test,],
               predict(cvfit, newx=ydat[idx_test,], s="lambda.1se", type="class"),
               predict(cvfit, newx=ydat[idx_test,], s="lambda.1se", type="response")
)
colnames(gscores)[7:8]<-c('pre_class','pre_score')
gscores$class<-as.numeric(as.character(gscores$class))

library(pROC) #calculate ROC
auc(gscores$class,gscores$pre_score)
auc(tscores$class,tscores$pre_score)

library(plotROC) #geom_roc package
rocplot<-ggplot(ascores, aes(d=class, m =pre_score,color=group))+
  geom_roc(n.cuts=0)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA))+
  ggtitle("Performence of prediction model")+  
  annotate("text", x = .75, y = .25,
           label = paste("AUC =", round(0.775, 2))) +
  scale_x_continuous("Specificity",breaks = seq(0, 1, by = .25),labels = sort(seq(0,1,by = 0.25),decreasing = T))+
  scale_y_continuous("Sensitivity",breaks = seq(0, 1, by = .25),labels = seq(0,1,by = 0.25))

#################################################
####------genes consisting of the model------####
#################################################
coef(cvfit,cvfit$lambda.1se)
lambda.1se.num<-which(cvfit$glmnet.fit$lambda==cvfit$lambda.1se)
beta.1se<-cvfit$glmnet.fit$beta[,lambda.1se.num]
beta.1se<-beta.1se[beta.1se!=0]
exp.1se<-ydat[,names(beta.1se)]

exp.1se<-as.data.frame(cbind(exp.1se,class=as.numeric(as.character(colData$class))))
exp.1se<-melt(data = exp.1se,id.vars=c("class"),variable.name="gene",value.name="expression")
exp.1se$class<-as.factor(exp.1se$class)

pdf("beta_exp_diff.pdf",height = 2,width = 4)
 ggplot(exp.1se, aes(x=class, y=expression)) + 
  geom_boxplot(aes(fill = class), alpha = 0.5,show.legend = T)+
  facet_grid(~gene,scales = "free",space = "free")
dev.off()

#################################################
####------------Heatmap plots----------------####
#################################################
library(pheatmap)
pheatmap(t(ydat), scale="none",clustering_distance_row="correlation",annotation=colData,
         show_colnames=F,cellwidth=6,cellheight =12,
         fontsize=6, fontsize_row=5)
#-important genes-#
pheatmap(t(ydat[,names(beta.1se)]), scale="none",clustering_distance_row="correlation",annotation=colData,
         show_colnames=F,cellwidth=6,cellheight =18,
         fontsize=6, fontsize_row=5)
