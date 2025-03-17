#############################################################################

#gsva测试
library(GSVA)
p <- 20000    ## number of genes
n <- 30       ## number of samples
nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 100 ## maximum gene set size
X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
dim(X)
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets
es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
pheatmap::pheatmap(es.max)
pheatmap::pheatmap(es.dif)

######ssgsea
setwd("D:\\项目\\TCGA_PAAD\\GEO\\GSE62452_array缺基因")
list.files()
#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
df<-read.delim("GSE62452.exp.matrix.xls",header=T,row.names = 1)
sur.data<-read.delim("GSE62452.survival.Info.txt",header=T,row.names = 1)
df<-df[,rownames(sur.data)]
exp<- as.matrix(df)    #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错

#emt.list<-read.delim("C:\\项目\\TCGA_PAAD\\GEO\\EMT_210keygene.txt",header=T,row.names = NULL)
#emt.list<-read.delim("C:\\项目\\TCGA_PAAD\\GEO\\EMT_15keygene.txt",header=T,row.names = NULL)
emt.list<-read.delim("C:\\项目\\TCGA_PAAD\\GEO\\132keygene.txt",header=T,row.names = NULL)
emt.list<-as.list(emt.list)

#进行gsva分析
re <- gsva(exp, emt.list, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)        #注意表达谱exp载入后需转化为matrix，前面已转换

save(re,file="emt.score.ssgsea.RData")

####plot
emt.group <-  as.vector(ifelse(re > median(re), "high", "low"))#定义低、高风险

dat <- cbind(t(re), EMTgroup=as.vector(emt.group)) 
head(dat)
colnames(dat)<-c("EMT_score","EMTgroup")
dat<-as.data.frame(dat)
dat$EMT_score<-as.numeric(dat$EMT_score)
#write.table(dat, file = "TCGA_PAAD_EMT.score.txt",sep = "\t",
 #           row.names = T,col.names = NA,quote = F)

#dat=read.delim("TCGA_PAAD_EMT.score.txt",header=T, row.names=1)

head(dat)

dat.1<-cbind(dat,sur.data[,c("survival_status","survival_months","tissue.ch1")])
colnames(dat.1)<-c("EMT_score","EMTgroup","OS","OS.time","sampletype")
head(dat.1)
df_t<-as.data.frame(t(df))
df_t<-df_t[rownames(dat.1),]
df_test<-cbind(dat.1[,c("OS","OS.time")],df_t)
df_test[1:4,1:6]
write.table(df_test, file = "df_test.sur.gse62452.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)

df_module<-cbind(dat.1,df_t[,c("ANLN","CD36","CELSR3","IL20RB","VSIG1")])
df_module<-cbind(dat.1,df_t[,c("ECT2","FBXL16","KCNB1","MET")])
df_module<-cbind(dat.1,df_t[,c("AGBL4","FGF9","GPR6","NPAS4","ASIC4")])
df_module<-cbind(dat.1,df_t[,c("ASCL1","ATP1A3","DSCAM","KCNB1","NPAS4","SLC25A47")])
df_module<-cbind(dat.1,df_t[,c("BIRC5","DEPDC1","E2F1","FERMT1","RAD54L","SMOC1")])

df_module<-cbind(dat.1,df_t[,c("DEPDC1","E2F1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")])
df_module<-cbind(dat.1,df_t[,c("DEPDC1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")])
df_module<-cbind(dat.1,df_t[,c("DEPDC1B","HES2", "AGER","EPHA2")])

df_module<-cbind(dat.1,df_t[,c("DEPDC1B","HES2", "AGER","MUC16")])

save(df_module,file="df_module.GSE62452.RData")
#######预后风险评分验证
library("dplyr")
library(survival)
library(ROCR)   
library(glmnet)
library(caret)

#tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+VSIG1, data = df_module)####final

#森林图
library(survminer)
ggforest(tdmultiCox,data = df_module) #画个森林图看看
p<-ggforest(tdmultiCox,data = td.test)
ggsave(filename = "5gene.ggforest.pdf", 
       plot =p,width = 13, height =6, units = 'cm')
#计算模型风险评分

library(dplyr)
#genes = names(tdmultiCox$coefficients);length(genes)
genes<-c("NT5E","CCDC74B","KCNB1","FBXL16","ECT2","MET")
genes<-c("ECT2","FBXL16","KCNB1","MET")
genes<-c("AGBL4","FGF9","GPR6","NPAS4","ASIC4")
genes<-c("ASCL1","ATP1A3","DSCAM","KCNB1","NPAS4","SLC25A47")
genes<-c("BIRC5","DEPDC1","E2F1","FERMT1","RAD54L","SMOC1")

genes<-c("DEPDC1","E2F1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")
genes<-c("DEPDC1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")

genes<-c("DEPDC1B","HES2", "AGER","EPHA2")

genes<-c("DEPDC1B","HES2", "AGER","MUC16")

load("C:\\项目\\TCGA_PAAD\\4gene.final.tdmultiCox.RData")
fp = apply(df_module[,genes], 1,function(k)sum(tdmultiCox$coefficients * k)) #lasso回归模型的预测值就是线性加乘
#fp=apply(df_module[,genes], 1,function(k)sum(actCoef * k))
df_module$riskscore = fp
write.table(df_module,file="df_module.EMTgroup.addriskscore.txt",sep="\t",row.names=T,quote=F)

#colnames(dat.1)<-c("EMT_score","EMTgroup","OS","OS.time","ANLN","CD36","CELSR3","IL20RB","VSIG1","riskscore")
riskscore.group <-  as.vector(ifelse(df_module$riskscore > median(df_module$riskscore), "high", "low"))
df_module <- cbind(df_module, Riskgroup=as.vector(riskscore.group))
head(df_module)
save(df_module,file="df_module.GSE62452.Rdata")
write.table(df_module,file="df_module.EMTgroup.addriskscore.txt",sep="\t",row.names=T,quote=F)
####输出选择的变量以及风险评分
#基于上述风险评估得分riskScore（连续变量）和癌症患者的结局OS（二分类），可以绘制箱线图和ROC曲线，来评估Lasso选择出来的变量进行风险预测的效果如何。
  library(ggpubr) 
  ##箱线图
  #绘制x为甲状腺癌结局，y为风险评估得分的箱图，并进行统计学检验
df_module$EMTgroup<-factor(df_module$EMTgroup,levels=c("low","high"))
  m <- ggboxplot(df_module, x = "EMTgroup", y = "EMT_score",
                 color = "EMTgroup", palette = "jco",
                 add = "jitter")

  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "EMTgroup.emtscore.boxplot.pdf", width = 2.8, height = 3.2) 
  
  m <- ggboxplot(df_module, x = "OS", y = "EMT_score",
                 color = "OS", palette = "jco",
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "EMTgroup.OS.boxplot.pdf", width = 2.8, height = 3.2) 
  df_module$Riskgroup<-factor(df_module$Riskgroup,levels=c("low","high"))
  m <- ggboxplot(df_module, x = "Riskgroup", y = "EMT_score",
                 color = "Riskgroup", palette = c("RoyalBlue","Red"),
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "Riskgroup.EMT_score.boxplot.pdf", width = 2.8, height = 3.2) 
#############################################################################################
  #生存分析
  #data.final=read.table("PAAD.data.for.lasso.txt",header=T, row.names=1)
  #head(data.final)
  #dim(data.final)
  #dat.sur<-data.final
  #df_module$OS.time.year=df_module$OS.time/365
  ###五年生存状态
  #df_module1$OS.5<-df_module$OS
  #df_module$OS.5[(df_module$OS.time.year > 5)] <-  0
  load("df_module.GSE62452.Rdata")
  library(survminer)
  library(survival)
  df_module$OS.time.year<-df_module$OS.time/12
  
    #dat.sur$gene[1:6]
  df_module$EMTgroup<-factor(df_module$EMTgroup,levels=c("high","low"))
    fit2 <- survfit(Surv(OS.time, OS) ~ EMTgroup, data = df_module)
    df_module$OS.time<-as.numeric(df_module$OS.time) 
    # palette = "lancet"
    p<-ggsurvplot(
      fit2,
      data = df_module,
      censor.shape="|", censor.size = 4,
      conf.int = TRUE,
      conf.int.style = "ribbon",
      conf.int.alpha = 0.2,
      pval = TRUE,
      #surv.median.line = "hv",
      ggtheme =theme_classic(),
      legend = "top",
      xlab = "OS_time(months)",
      ylab = "Survival probability",
      title = "Survival curves",
      break.x.by = 10,
      break.y.by = 0.2,
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.height = 0.2,
      risk.table.y.text = FALSE
    )
    p
    pdf(file="Survival.curves.EMTgroup.pdf",4.5,5.8)
    p
    dev.off()
###riskscore
head(df_module)

library(survminer)
library(survival)

df_module$Riskgroup<-factor(df_module$Riskgroup,levels=c("high","low"))
fit2 <- survfit(Surv(OS.time.year, OS) ~ Riskgroup, data = df_module)
df_module$OS.time.year<-as.numeric(df_module$OS.time.year) 
# palette = "lancet"
p<-ggsurvplot(
  fit2,
  data = df_module,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  #surv.median.line = "hv",
  ggtheme =theme_classic(),
  legend = "top",
  xlab = "OS_time(years)",
  ylab = "Survival probability",
  title = "Survival curves",
  break.x.by = 1,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)
p
pdf(file="Survival.curves.Riskgroup.years.pdf",4.5,5.8)
p
dev.off()

#联合
head(df_module)
union.group <-  
  as.vector(ifelse(df_module$EMTgroup=="high" ,
                   ifelse(df_module$Riskgroup=="high",'group1','group2'),
                   ifelse(df_module$Riskgroup=="high",'group3','group4')))

df_module <- cbind(df_module, Uniongroup=as.vector(union.group))
head(df_module)
write.table(df_module,file="df_module.EMTgroup.addriskscore.txt",sep="\t",row.names=F,quote=F)

library(survminer)
library(survival)
head(df_module)
df_module$Uniongroup<-factor(df_module$Uniongroup,
                             levels=c("group1","group2","group3","group4"))
df_module.1<-df_module %>% subset(Uniongroup=="group1" | Uniongroup=="group4")
df_module.1$Uniongroup<-factor(df_module.1$Uniongroup,
                             levels=c("group4","group1"))
fit2 <- survfit(Surv(OS.time, OS) ~ Uniongroup, data = df_module.1)
df_module.1$OS.time<-as.numeric(df_module.1$OS.time) 
# palette = "lancet"
p<-ggsurvplot(
  fit2,
  data = df_module.1,
  palette="lancet",
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  #surv.median.line = "hv",
  ggtheme =theme_classic(),
  legend = "top",
  xlab = "OS_time(months)",
  ylab = "Survival probablity",
  title = "Survival curves",
  break.x.by = 10,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)
p
pdf(file="Survival.curves.Uniongroup.pdf",4.5,5.8)
p
dev.off()

################################################################

tf<-read.delim("Homo_sapiens_TF.txt",header=T,sep='\t')
intersect(comm_deg,unique(tf$Symbol))
