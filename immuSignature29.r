setwd("D:\\项目\\TCGA_PAAD\\immuncellAI\\immuSignature29")
library(GSVA)
load("D:/项目/TCGA_PAAD/GEO/dat.1.Rdata")
head(dat.1)
#EMT_score EMTgroup OS OS.time   riskscore Riskgroup
#TCGA-HV-A7OP  0.60605467      low  0     978 0.350578789       low
#TCGA-F2-A8YN  0.74562514     high  0     517 0.628891477      high
dim(dat.1)
#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
df<-read.delim("D:\\项目\\TCGA_PAAD\\mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
exp<- as.matrix(df)    #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错

list.files()
#emt.list<-read.delim("C:\\项目\\TCGA_PAAD\\GEO\\EMT_210keygene.txt",header=T,row.names = NULL)
#emt.list<-read.delim("C:\\项目\\TCGA_PAAD\\GEO\\EMT_15keygene.txt",header=T,row.names = NULL)
sig.list<-read.delim("29-1.txt",header=T,row.names = NULL,fill=TRUE,strip.white = TRUE)
sig.list<-as.list(sig.list)

#进行gsva分析
re <- gsva(exp, sig.list, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)        #注意表达谱exp载入后需转化为matrix，前面已转换


list<-substring(colnames(re),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
colnames(re)<-list1[,1]

save(re,file="Sig29-1.score.ssgsea.RData")
sig.1<-t(re)
head(sig.1);dim(sig.1)
#
sig.list<-read.delim("29-2.txt",header=T,row.names = NULL,fill=TRUE,strip.white = TRUE)
sig.list<-as.list(sig.list)

#进行gsva分析
re2 <- gsva(exp, sig.list, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)        #注意表达谱exp载入后需转化为matrix，前面已转换


list<-substring(colnames(re2),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
colnames(re2)<-list1[,1]

save(re2,file="Sig29-2.score.ssgsea.RData")
sig.2<-t(re2)
head(sig.2);dim(sig.2)
sig.final<-cbind(sig.1,sig.2)
#####add riskscore
load("D:/项目/TCGA_PAAD/GEO/dat.1.Rdata")
head(dat.1)
sig.1<-sig.1[rownames(dat.1),]
sig.2<-sig.2[rownames(dat.1),]
sig.1.final<-cbind(sig.1,dat.1)
sig.2.final<-cbind(sig.2,dat.1)
save(sig.1.final,file="sig.1.final.RData")
save(sig.2.final,file="sig.2.final.RData")

#####plot1
library("reshape")
dat_melt<-melt(sig.1.final,id=c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup"))
colnames(dat_melt)<-c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup","celltype","Score")
head(dat_melt)

library("ggpubr")
library(RColorBrewer)
display.brewer.all()
mycol <- colorRampPalette(brewer.pal(10,'Paired'))(10)
mycol <- colorRampPalette(brewer.pal(3,'Set2'))(3)

dat_melt$Score<-as.numeric(dat_melt$Score)
dat_melt$Riskgroup<-factor(dat_melt$Riskgroup,levels=c("high","low"))
p<-ggboxplot(dat_melt, x = "celltype", y = "Score",color = "Riskgroup",
             palette = c("IndianRed1","MediumTurquoise"),
             show.legend = FALSE, short.panel.labs = TRUE)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.8,vjust =1, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Propotion",title="")

my_comparisons <- list( c("high", "low"))

#p2<-p+facet_wrap(~ celltype, scales="free",nrow = 6) + 
# stat_compare_means(comparisons = my_comparisons,label = "p.format",
#                    size=4,method="wilcox.test")+
#  theme(strip.text.x = element_text(size=11, color="black",
#                                   face="bold"))+
#  theme(strip.background = element_blank(),strip.placement = "outside")
#p2

p2<-p +stat_compare_means(aes(group = Riskgroup),label = "p.signif",
                          size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_sig29-1.score.riskgroup.4gene.pdf",plot=p2,
       device='pdf',path=".",width=7,height=3.5)
##########################################
#####plot2
library("reshape")
dat_melt<-melt(sig.2.final,id=c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup"))
colnames(dat_melt)<-c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup","celltype","Score")
head(dat_melt)

library("ggpubr")
library(RColorBrewer)
display.brewer.all()
mycol <- colorRampPalette(brewer.pal(10,'Paired'))(10)
mycol <- colorRampPalette(brewer.pal(3,'Set2'))(3)

dat_melt$Score<-as.numeric(dat_melt$Score)
dat_melt$Riskgroup<-factor(dat_melt$Riskgroup,levels=c("high","low"))
p<-ggboxplot(dat_melt, x = "celltype", y = "Score",color = "Riskgroup",
             palette = c("IndianRed1","MediumTurquoise"),
             show.legend = FALSE, short.panel.labs = TRUE)+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 0.8,vjust =1, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="Propotion",title="")

my_comparisons <- list( c("high", "low"))

#p2<-p+facet_wrap(~ celltype, scales="free",nrow = 6) + 
# stat_compare_means(comparisons = my_comparisons,label = "p.format",
#                    size=4,method="wilcox.test")+
#  theme(strip.text.x = element_text(size=11, color="black",
#                                   face="bold"))+
#  theme(strip.background = element_blank(),strip.placement = "outside")
#p2

p2<-p +stat_compare_means(aes(group = Riskgroup),label = "p.signif",
                          size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_sig29-2.score.riskgroup.4gene.pdf",plot=p2,
       device='pdf',path=".",width=6,height=3.5)
