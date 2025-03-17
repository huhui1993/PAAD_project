install.packages("oncoPredict")
####教程：https://mp.weixin.qq.com/s/QRaTd-fIsqq6sPsLmOPvIw
setwd("D:\\项目\\TCGA_PAAD\\药物敏感性分析")
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

##因为我们这个是教程，所以我就不读取自己的表达量矩阵了，直截了当的从Genomics of Drug Sensitivity in Cancer (GDSC) 
#的v2里面随机挑选10个细胞系作为要预测的矩阵。
testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
testExpr[1:4,1:4]  
colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)  
#了训练集的表达量矩阵和药物处理信息，然后也有了待预测的表达量矩阵，接下来就是一个函数的事情啦！
#这个函数calcPhenotype就是R包 oncoPredict的核心
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

####用TCGA数据进行预测
df<-read.delim("D:\\项目\\TCGA_PAAD\\mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
list<-substring(colnames(df),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
colnames(df)<-list1[,1]
df[1:3,1:4]
df<-as.matrix(df)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = df,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

result<-read.delim("D:\\项目\\TCGA_PAAD\\药物敏感性分析\\calcPhenotype_Output\\DrugPredictions.txt",
                   header=T,row.names = 1)

result.1<-t(result)
dim(result.1);result.1[1:3,1:3]
result.2<-na.omit(result.1)
dim(result.2)
#载入风险值和分组
load("D:\\项目\\TCGA_PAAD\\GEO\\dat.1.Rdata")
head(dat.1);dim(dat.1)
#EMT_score EMTgroup OS OS.time   riskscore Riskgroup
#TCGA-IB-AAUU 0.8129067     high  0     245 -0.05490632      high
#TCGA-3A-A9IL 0.1543919      low  0    2741 -0.50445004       low
#TCGA-3A-A9IX 0.6741788      low  0    1037 -0.14180242       low

comm<-intersect(rownames(dat.1),colnames(result.2))
result.3<-result.2[,comm]
dat<-dat.1[comm,]
result.final<-cbind(dat,t(result.3))


result.final<-result.final[order(result.final$Riskgroup,result.final$riskscore),]
dim(result.final);
head(result.final)
####heatmap
library("ComplexHeatmap")
library(circlize)
library(methods)
library(magrittr)

set.seed(12345)

col_fun <-  c("#ff7f00", "#1f78b4")

column_ha = HeatmapAnnotation(riskgroup = result.final$Riskgroup,
                              col = list(
                                riskgroup = c("high" = "IndianRed1", 
                                              "low" = "MediumTurquoise")
                              ))
mycol <- colorRamp2(c( -1, 0, 1), c("#ab7bf9","white","Yellow"))
data1<-t(result.final[,c(7:160)])

write.table(data1,file="GDSC2.drug.sensitivity.order.txt",
            sep="\t",row.names=T,quote=F)
#
colnum<-ncol(data1)
rownum<-nrow(data1)
pvalue<-0
for(n in 1:rownum){
  data_case<-as.numeric(data1[n,1:88])
  data_con<-as.numeric(data1[n,89:177])
  tempresult<-try(fit1<-wilcox.test(data_case,data_con, paired=FALSE),silent=TRUE)
  pvalue[n]<-tempresult$p.value
}
pvalue<-as.matrix(pvalue)
rownames(pvalue)<-rownames(data1)
write.table(pvalue,"GDSC2.drug.sensitivity.wilcox.test.pvalue.txt",
            quote = F,sep='\t')
#
tmp<-read.delim("D:\\项目\\TCGA_PAAD\\药物敏感性分析\\GDSC2.drug.sensitivity.order.filter.new.txt",
                   header=T,row.names = 1)
dim(tmp);
tmp[1:3,1:4]

data1<-data1[rownames(tmp),]

data.m<-apply(data1,1,scale)
data.m<-t(data.m)
colnames(data.m)<-colnames(data1)

write.table(data.m,file="GDSC2.drug.sensitivity.zscore.txt",
            sep="\t",row.names=T,quote=F)

mycol <- colorRamp2(c( -0.5, 0, 0.5), c("#ab7bf9","white","Yellow"))
pdf(file="heamap.GDSC2.drug.sensitivity.pdf",width=7,height=10)
Heatmap(as.matrix(data.m),cluster_columns = FALSE,cluster_rows = TRUE,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "average",
        col=mycol,top_annotation = column_ha,
        name = "Drug\nsensitivity",
        row_dend_width = unit(4, "mm"),column_names_gp=gpar(fontsize = 11),
        column_names_side = c("bottom"),row_names_side = c("left"),
        row_names_gp=gpar(fontsize = 10),
        row_title_gp = gpar(col = "black",fontsize = 12),
        show_row_names = TRUE,
        show_column_names = FALSE)
dev.off()


####boxplot
library(ggpubr) 
library("ggpubr")
library(RColorBrewer)
display.brewer.all()
mycol <- colorRampPalette(brewer.pal(10,'Paired'))(10)
mycol <- colorRampPalette(brewer.pal(3,'Set2'))(3)

drug.list<-c("Teniposide_1809","Mitoxantrone_1810","Fludarabine_1813",
             "ML323_1629","Wnt.C59_1622","Axitinib_1021","Uprosertib_2106",
             "Ibrutinib_1799","Navitoclax_1011","AZD1208_1449",
             "Linsitinib_1510","Ulixertinib_2047")

head(result.final)
data.tmp<-result.final[,-(1:5)]
head(data.tmp)
data.tmp<-data.tmp[,c("Riskgroup",drug.list)]
head(data.tmp)

library("reshape")
dat_melt<-melt(data.tmp,id=c("Riskgroup"))
head(dat_melt)
colnames(dat_melt)<-c("Riskgroup","type","value")
# palette =c(mycol),
#dat_melt.1<-dat_melt %>% subset(celltype %in% c("T cell CD4+ memory","T cell CD4+ naive"))
dat_melt$value<-as.numeric(dat_melt$value)
dat_melt$Riskgroup<-factor(dat_melt$Riskgroup,levels=c("high","low"))
p<-ggboxplot(dat_melt, x = "type", y = "value",color = "Riskgroup",
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
        axis.title.x = element_text(size=12,colour="black",face="bold"))+ylim(0,300)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="type",x="",y="drug sensitivity score",title="")

my_comparisons <- list( c("high", "low"))


p2<-p +stat_compare_means(aes(group = Riskgroup),label = "p.signif",
                          size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_risk.score.drug.sensitivity.score.pdf",plot=p2,
       device='pdf',path=".",width=7,height=4.5)
