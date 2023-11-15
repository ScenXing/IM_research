library(DESeq2)
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(pheatmap)

#---function: change ENSEMBL to ENTREZID---#
ensg2entrez<-
  function(fc2){
    fc2$ENSEMBL<-grep('ENSG',unlist(strsplit(rownames(fc2),"\\.")),value=T)
    fc2<-fc2[,c('ENSEMBL','log2FoldChange')]
    fc2.id <- bitr(fc2$ENSEMBL,fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
    fc2<-merge(fc2.id,fc2)
    fc2<-arrange(fc2,desc(log2FoldChange))
    fc2sort<-fc2$log2FoldChange
    names(fc2sort)<-fc2$ENTREZID
    fc2sort<-fc2sort[!is.na(fc2sort)]  #remove genes without foldchange value
    rm(fc2,fc2.id)
    fc2sort}

#------IM vs Control------#
row.names(countData)<-grep('ENSG',unlist(strsplit(rownames(countData),"\\.")),value=T)
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,
                              design=~gender+batch+age+class)
dds$class<- relevel(dds$class, ref='0')
dds <- DESeq(dds)
resBac <- results(dds,contrast=c("class","1","0"))
resBac <- resBac[order(resBac$padj),]
fc2bac=as.data.frame(resBac)
degBac<-fc2bac[abs(fc2bac$log2FoldChange)>0.5 & fc2bac$padj<0.05 &!is.na(fc2bac$padj),]
write.table(degBac,"BM_vs_CTRL_DEGs.tsv",sep = "\t",quote = F) 

resVir <- results(dds,contrast=c("class","2","0"))
resVir <- resVir[order(resVir$padj),]
fc2vir=as.data.frame(resVir)
degVir<-fc2vir[abs(fc2vir$log2FoldChange)>0.5 & fc2vir$padj<0.05 &!is.na(fc2vir$padj),]
write.table(degVir,"VM_vs_CTRL_DEGs.tsv",sep = "\t",quote = F)

#--retrieve normalized gene counts-#
genenor <- assay(vst(dds, blind=FALSE))

#------------Bacteria vs Ctrl---------------#
#top DEGs
bactopgene<-row.names(fc2bac[abs(fc2bac$log2FoldChange)>0.5 & fc2bac$padj<0.05 &!is.na(fc2bac$padj),])[1:20]
bactop<-genenor[bactopgene,colData$class %in% c(0,1)]
hp1bac<-pheatmap(bactop,show_colnames = F,annotation_col = colBac[,c(1,2,5)],fontsize_row=5,
                 treeheight_col=20,treeheight_row=10)

#---GO KEGG enichments---#
bacgenes<-ensg2entrez(fc2bac)
bacGO<-gseGO(bacgenes,ont = 'BP',keyType = "ENTREZID",OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05)
bacGO<-setReadable(bacGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #readable
bacKEGG<-gseKEGG(bacgenes,organism = 'hsa',pvalueCutoff = 0.05)
bacKEGG<-setReadable(bacKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #readable
View(as.data.frame(bacKEGG)) #48 terms
#----plot of KEGG terms----#
p1bac<-dotplot(bacKEGG,x='p.adjust',size='Count',color='p.adjust',font.size=8)


#----------------Virus vs Ctrl--------------------#
#top DEGs
virtopgene<-row.names(fc2vir[abs(fc2vir$log2FoldChange)>0.5 & fc2vir$pvalue<0.05 & !is.na(fc2vir$padj),])[1:16]
virtop<-genenor[virtopgene,colData$class %in% c(0,2)]
hp1vir<-pheatmap(virtop,show_colnames = F,annotation_col = colVir[,c(1,2,5)],fontsize_row=5,
                 treeheight_col=20,treeheight_row=10)

#---GO_KEGG_enrichment---#
virgenes<-ensg2entrez(fc2vir)
virGO<-gseGO(virgenes,ont = 'BP',keyType = "ENTREZID",OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05)
virGO<-setReadable(virGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #readable
virKEGG<-gseKEGG(virgenes,organism = 'hsa',pvalueCutoff = 0.05)
virKEGG<-setReadable(virKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #readable
View(as.data.frame(virKEGG)) #40 terms

#----Plot of KEGG terms----#
p1vir<-dotplot(virKEGG,x='p.adjust',size='Count',color='p.adjust',font.size=8)
