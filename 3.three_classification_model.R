library(DESeq2) # 1.12.4
library(glmnet) # 2.0-5
library(pROC)   # ROC curve
library(pheatmap) # heatmap
library(DaMiRseq) # feature selections
library("biomaRt") # Retrive information for ensembl
library(hash) # Perl hash
library("reshape2") # melt function
library(ggplot2) # plotting

###########################################
######--------machine learning-------######
###########################################
rdat<-as.matrix(data_reduced)
#--generate model by cross validation--#
cvfit <- cv.glmnet(x = rdat, y = colData$class,family = "multinomial",
                   type.measure = "class",nfolds =30)
ascores<-cbind(colData,
               predict(cvfit, newx=rdat, s="lambda.1se", type="class"),
               predict(cvfit, newx=rdat, s="lambda.1se", type="response"))
colnames(ascores)[6:9]<-c('pred_class',0,1,2)

#----plot ROC curve by pROC package----#
#------plot ROC curve----#
library(multiROC)
#----data.frame for multiROC package---#
nclass<-matrix(nrow = nrow(ascores),ncol = 3)
for (i in 1:nrow(ascores)){
  x<-as.numeric(as.character(ascores[i,5]))
  if(x==0){
    nclass[i,]<-c(1,0,0)
  }else if(x==1){
    nclass[i,]<-c(0,1,0)
  }else if(x==2){
    nclass[i,]<-c(0,0,1)
  }
}
nclass<-cbind(nclass,ascores[,7:9])
colnames(nclass)<-c('ctrl_true','bacteria_true','virus_true','ctrl_pred_m1','bacteria_pred_m1','virus_pred_m1')
res <- multi_roc(nclass, force_diag=T)

#---------plot ROC curve-----------#
#-Change the format of results to a ggplot2 friendly format.
n_method <- length(unique(res$Methods))
n_group <- length(unique(res$Groups))
res_df <- data.frame(Specificity= numeric(0), Sensitivity= numeric(0), Group = character(0), AUC = numeric(0), Method = character(0))
for (i in 1:n_method) {
  for (j in 1:n_group) {
    temp_data_1 <- data.frame(Specificity=res$Specificity[[i]][j],
                              Sensitivity=res$Sensitivity[[i]][j],
                              Group=unique(res$Groups)[j],
                              AUC=res$AUC[[i]][j],
                              Method = unique(res$Methods)[i])
    colnames(temp_data_1) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
    res_df <- rbind(res_df, temp_data_1)
    
  }
  temp_data_2 <- data.frame(Specificity=res$Specificity[[i]][n_group+1],
                            Sensitivity=res$Sensitivity[[i]][n_group+1],
                            Group= "Macro",
                            AUC=res$AUC[[i]][n_group+1],
                            Method = unique(res$Methods)[i])
  temp_data_3 <- data.frame(Specificity=res$Specificity[[i]][n_group+2],
                            Sensitivity=res$Sensitivity[[i]][n_group+2],
                            Group= "Micro",
                            AUC=res$AUC[[i]][n_group+2],
                            Method = unique(res$Methods)[i])
  colnames(temp_data_2) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
  colnames(temp_data_3) <- c("Specificity", "Sensitivity", "Group", "AUC", "Method")
  res_df <- rbind(res_df, temp_data_2)
  res_df <- rbind(res_df, temp_data_3)
}
#-4.2 Plot
library(ggplot2)
fig.roc<-
  ggplot(res_df[res_df$Group %in% c('bacteria','virus','ctrl'),],aes(x = 1-Specificity, y=Sensitivity)) + 
  geom_path(aes(color = Group)) + 
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.justification=c(1, 0), 
        legend.position=c(.95, .05), legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, linetype="solid", colour ="black"))


#################################################
###---Obtain genes consisting of the model--#####
#################################################
lambda.1se.num<-which(cvfit$glmnet.fit$lambda==cvfit$lambda.1se) #lambda.1se(used lambda)

beta.1se.t1<-cvfit$glmnet.fit$beta$`1`[,lambda.1se.num] #find the used beta list 
beta.1se.t1<-beta.1se.t1[beta.1se.t1!=0] 
exp.1se.t1<-t(rdat)[names(beta.1se.t1),]

beta.1se.t2<-cvfit$glmnet.fit$beta$`2`[,lambda.1se.num] #find the used beta list 
beta.1se.t2<-beta.1se.t2[beta.1se.t2!=0] 
exp.1se.t2<-t(rdat)[names(beta.1se.t2),]

#--------genes in the glmnet model---------#
#-convert ensg id to gene symbol
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='http://grch37.ensembl.org')
genetmp<-
  getBM(attributes = c("external_gene_name","chromosome_name", "start_position","ensembl_gene_id","gene_biotype"),
        filters = "ensembl_gene_id",
        values = beta.1se.gene,
        mart=mart)
genetmp<-genetmp[!grepl('_',genetmp$chromosome_name),]
setdiff(beta.1se.gene,genetmp$ensembl_gene_id) #unmatched ID
hubgenes<-genetmp
ensg2sym <- hash( keys=hubgenes$ensembl_gene_id, values=hubgenes$external_gene_name )

#--genes in class(0,1,2) models--#
t1<-c()
for (i in 1:nrow(exp.1se.t1)){
  t1<-c(t1,which(hubgenes$ensembl_gene_id==row.names(exp.1se.t1)[i]))
}
row.names(exp.1se.t1)<-hubgenes$external_gene_name[t1]
exp.1se.t1<-as.data.frame(cbind(t(exp.1se.t1),class=as.numeric(as.character(colData$class))))
exp.1se.t1<-melt(data = exp.1se.t1,id.vars=c("class"),variable.name="gene",value.name="expression")
exp.1se.t1$class<-as.factor(exp.1se.t1$class)
pdf("t1_exp.pdf",height = 4,width = 12)
ggplot(exp.1se.t1, aes(x=class, y=expression)) + 
  geom_boxplot(aes(fill = class), alpha = 0.5,show.legend = T)+
  facet_grid(~gene,scales = "free",space = "free")
dev.off()


t2<-c()
for (i in 1:nrow(exp.1se.t2)){
  t2<-c(t2,which(hubgenes$ensembl_gene_id==row.names(exp.1se.t2)[i]))
}
row.names(exp.1se.t2)<-hubgenes$external_gene_name[t2]
exp.1se.t2<-as.data.frame(cbind(t(exp.1se.t2),class=as.numeric(as.character(colData$class))))
exp.1se.t2<-melt(data = exp.1se.t2,id.vars=c("class"),variable.name="gene",value.name="expression")
exp.1se.t2$class<-as.factor(exp.1se.t2$class)

#####################################################
###----------Obtain genes in the model-----------####
#####################################################
pcg.t1<-cvfit$glmnet.fit$beta$`1`[,lambda.1se.num]
pcg.t1<-pcg.t1[pcg.t1>0 & names(pcg.t1) %in% pcg]
exp.pcg.t1<-t(rdat)[names(pcg.t1),]

pcg.t2<-cvfit$glmnet.fit$beta$`2`[,lambda.1se.num]
pcg.t2<-pcg.t2[pcg.t2>0 & names(pcg.t2) %in% pcg]
exp.pcg.t2<-t(rdat)[names(pcg.t2),]

pcg.gene<-unique(c(names(pcg.t1),names(pcg.t2)))

#--genes in class(0,1,2) models--#
t1<-c()
for (i in 1:nrow(exp.pcg.t1)){
  t1<-c(t1,which(hubgenes$ensembl_gene_id==row.names(exp.pcg.t1)[i]))
}
row.names(exp.pcg.t1)<-hubgenes$external_gene_name[t1]
exp.pcg.t1<-as.data.frame(cbind(t(exp.pcg.t1),class=as.numeric(as.character(colData$class))))
exp.pcg.t1<-melt(data = exp.pcg.t1,id.vars=c("class"),variable.name="gene",value.name="expression")
exp.pcg.t1$class<-as.factor(exp.pcg.t1$class)

#-Three important genes-#
exp.pcg.t1<-exp.pcg.t1[exp.pcg.t1$gene %in% c('ASRGL1','NR2F6','OLFML3'),]
pdf("t1_exp.pcg.3genes.pdf",height = 2,width = 5)
ggplot(exp.pcg.t1, aes(x=class, y=expression)) + 
  geom_boxplot(aes(fill = class), alpha = 0.5,show.legend = T,outlier.shape = NA)+
  facet_grid(~gene,scales = "free",space = "free")+
  ylim(3.3,6)+theme_bw() + theme(panel.grid=element_blank())
dev.off()

t2<-c()
for (i in 1:nrow(exp.pcg.t2)){
  t2<-c(t2,which(hubgenes$ensembl_gene_id==row.names(exp.pcg.t2)[i]))
}
row.names(exp.pcg.t2)<-hubgenes$external_gene_name[t2]
exp.pcg.t2<-as.data.frame(cbind(t(exp.pcg.t2),class=as.numeric(as.character(colData$class))))
exp.pcg.t2<-melt(data = exp.pcg.t2,id.vars=c("class"),variable.name="gene",value.name="expression")
exp.pcg.t2$class<-as.factor(exp.pcg.t2$class)
#--3 genes--#
exp.pcg.t2<-exp.pcg.t2[exp.pcg.t2$gene %in% c('STIP1','PGAM5','AKAP8'),]
pdf("t2_exp_pcg_3_genes.pdf",height = 2,width = 5)
ggplot(exp.pcg.t2, aes(x=class, y=expression)) + 
  geom_boxplot(aes(fill = class), alpha = 0.5,show.legend = T,outlier.shape = NA)+
  facet_grid(~gene,scales = "free",space = "free")+
  ylim(3.5,6.3)+theme_bw() + theme(panel.grid=element_blank())
dev.off()

#-all genes-#
ggplot(exp.pcg.t2, aes(x=class, y=expression)) + 
  geom_boxplot(aes(fill = class), alpha = 0.5,show.legend = T)+
  facet_grid(~gene,scales = "free",space = "free")