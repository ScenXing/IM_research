library(tidyverse)
library(pkgsearch)
library(plotROC)
library(ggplot2)
library(cowplot)
library(plyr)
library(RColorBrewer)
library(readr)
#######################################
#####---Pathogen distributions---######
#######################################
#-plot 1-#
 ggplot(fig_data, aes(x=Pathogens, y=values, fill=class)) +
  geom_col(position="stack") +
  theme_classic() +
  scale_fill_manual(values=brewer.pal(2, "Set2")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(face='bold', size=20)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1)) +
  theme(text=element_text(family="serif")) +
  labs(x="Pathogens_class", y="Numbers", fill="Methods") +
  coord_flip()

#-plot 2-#
ggplot(long_p2_data,aes(d=D,m=M,color=name))+geom_roc(n.cuts = 0)+
  ggsci::scale_color_lancet()+style_roc(xlab = "1 - Specificity",ylab = "sensitivity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  fig2A+annotate("text",x=.75,y=.35,label=paste("M1(AUC) =",round(calc_auc(fig2A)$AUC[1],2)))+
  annotate("text",x=.75,y=.25,label=paste("M2(AUC) =",round(calc_auc(fig2A)$AUC[2],2)))+
  annotate("text",x=.75,y=.15,label=paste("M3(AUC) =",round(calc_auc(fig2A)$AUC[3],2)))

#-plot 3-#
ggplot(long_p3_data,aes(d=D,m=M,color=name))+geom_roc(n.cuts=0)+
  ggsci::scale_fill_lancet()+style_roc(xlab = "1 - Specificity",ylab = "sensitivity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  fig2B+annotate("text",x=.75,y=.25,label=paste("AUC of M3 =",round(calc_auc(fig2B)$AUC[1],2)))+
  annotate("text",x=.75,y=.15,label=paste("AUC of M4 =",round(calc_auc(fig2B)$AUC[2],2)))

#-plot 4-#
ggplot(p5_data,aes(EV,depth))+geom_area(fill="darkgreen",alpha = 0.5,color=1,lwd = 1,linetype = 1)+
  theme(axis.title.x=element_text(vjust=1, size=10,face = "bold"))+
  theme(axis.title.y=element_text(vjust=2, size=10,face = "bold"))+ylim(0,150)+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#-plot 5-#
ggplot(p6_data,aes(E30,depth))+geom_area(fill="darkgreen",alpha = 0.5,color = 1,lwd = 1,linetype = 1)+
  theme(axis.title.x=element_text(vjust=1, size=10,face = "bold"))+
  theme(axis.title.y=element_text(vjust=2, size=10,face = "bold"))+ylim(0,150)+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#-plog 6-#
ggplot(p6_data, aes(x=Pathogens, y=values, fill=class)) +
  geom_col(position="stack") +
  theme_classic() +
  scale_fill_manual(values=brewer.pal(2, "Set2")) +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(face='bold', size=20)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1)) +
  theme(text=element_text(family="serif")) +
  labs(x="Pathogens", y="Numbers", fill="Methods") +
  coord_flip() 