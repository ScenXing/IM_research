library(glmnet)
library(DaMiRseq)

library("biomaRt") #retrive information for ensembl
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
genetmp<-
  getBM(attributes = c("external_gene_name","chromosome_name", "ensembl_gene_id","gene_biotype"),
        filters = "ensembl_gene_id",
        values = row.names(countData),
        mart=mart)
genetmp<-genetmp[!grepl('_',genetmp$chromosome_name),]
pgene<-genetmp$ensembl_gene_id[genetmp$gene_biotype %in% 'protein_coding']

###########################################
#-------BM/CTRL feature selection---------#
###########################################
colBac<-colData[colData$class %in% c(0,1),]
colBac$class <-factor(colBac$class,levels = c(0,1)) 
countBac<-countData[,row.names(colBac)]
SE<-DaMiR.makeSE(countBac,colBac)
tSE <- SE[, tcohort, drop=FALSE]
vSE <- SE[, vcohort, drop=FALSE]
#-normalization
data_norm <- DaMiR.normalization(tSE, minCounts=10, fSample=0.5,hyper = "yes", th.cv=3,type="rlog")
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.5)
sv <- DaMiR.SV(data_filt)
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=8)
df<-colData(data_adjust)
DaMiR.corrplot(sv,df)

#-----feature selection-----#
data_clean<-DaMiR.transpose(assay(data_adjust))
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.75,th.VIP =3) # 
data_reduced <- DaMiR.FReduct(data_reduced$data,th.corr = 0.9) # remove highly correlated features

#-PCA figure
DaMiR.MDSplot(data_reduced, df,type="pearson")

df.importance <- DaMiR.FSort(data_reduced, df)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,n.pred =9) #ALL FEATURES
DaMiR.Clustplot(selected_features$data, df) #heatmap

###########################################
#-----Machine learning(glmnet package)----#
###########################################
tdat<-as.matrix(selected_features$data) 
#focus on high expressed genes in BM
#--generate model by cross validation--#
cvfit <- cv.glmnet(x = tdat, y = colBac$class[tcohort],family = "binomial",
                   type.measure = "class",nfolds =length(tcohort),
                   intercept=FALSE)
tscores<-cbind(colBac[tcohort,],
               predict(cvfit, newx=tdat, s="lambda.1se", type="class"),
               predict(cvfit, newx=tdat, s="lambda.1se", type="response"))

colnames(tscores)[6:7]<-c('pred_class','pre_score')

#--------prediction on the validation cohorts---------#
vdat<-t(v_adjust[colnames(tdat),])
vscores<-cbind(colBac[vcohort,],
               predict(cvfit, newx=vdat, s="lambda.1se", type="class"),
               predict(cvfit, newx=vdat, s="lambda.1se", type="response"))
colnames(vscores)[6:7]<-c('pred_class','pre_score')

###########################################
#------------plot ROC curve---------------#
###########################################
#-calculate ROC
library(pROC) 
auc(tscores$class,tscores$pre_score) #training
table(tscores$class,tscores$pred_class)
auc(vscores$class,vscores$pre_score) #test
table(vscores$class,vscores$pred_class)

#find the used beta list
coef(cvfit,cvfit$lambda.1se)

#---ROC for training and test cohort---#
library(plotROC) #geom_roc package
ascores<-rbind(tscores,vscores)
ascores$group<-rep(c('training','test'),c(nrow(tscores),nrow(vscores)))
ascores$class<-as.numeric(as.character(ascores$class))
rocplot<-
  ggplot(ascores, aes(d=class, m =pre_score,color=group))+
  geom_roc(n.cuts=0)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA))+
  theme(panel.grid=element_blank())+
  ggtitle("Performence of BM model")+  
  annotate("text", x = .75, y = .25) +
  scale_x_continuous("Specificity",breaks = seq(0, 1, by = .25),labels = sort(seq(0,1,by = 0.25),decreasing = T))+
  scale_y_continuous("Sensitivity",breaks = seq(0, 1, by = .25),labels = seq(0,1,by = 0.25))

###########################################
#-----------Suspected samples-------------#
###########################################
colnames(countPol)<-grep('CSF',unlist(strsplit(colnames(countPol),'_')),value=T) #change sample names
row.names(countPol)<-grep('ENSG',unlist(strsplit(rownames(countPol),"\\.")),value=T) #change gene names
countPol<-countPol[,row.names(colPol)]
colPol<-colPol[,c('batch','gender','age','collection_days','class')]
colPol$class<-'pol'
###########################################
#----Prediction of suspected samples------#
###########################################
#-prediction
preBM<-function(countPre,colPre){
  countPre2<-cbind(countPre,countBac[,vcohort])
  colPre2<-rbind(colPre,colBac[vcohort,])
  pSE<-DaMiR.makeSE(countPre2,colPre2)
  p_norm <- DaMiR.iTSnorm(tSE[rownames(data_norm)],pSE[rownames(data_norm)],normtype = "rlog",method = "precise")
  p_adjust <- DaMiR.iTSadjust(data_adjust,p_norm)
  pdat<-t(p_adjust[colnames(tdat),row.names(colPre)])
  pscores<-cbind(colPre,
                 predict(cvfit, newx=pdat, s="lambda.1se", type="class"),
                 predict(cvfit, newx=pdat, s="lambda.1se", type="response"))
  colnames(pscores)[6:7]<-c('pred_class','pre_score')  
  pscores
}
pscores<-preBM(countPol,colPol)