#############################################################################
######lasso
setwd("D:\\项目\\TCGA_PAAD")
df<-read.delim("mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
list<-substring(colnames(df),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
write.table(list1, file = "PAAD.haveRNA-Seq.sample.list.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)

exp.tmp<-t(df)
rownames(exp.tmp)<-list1[,1]

key358gene<-read.delim("C:\\项目\\TCGA_PAAD\\WGCNA\\immuneAndEMTgene_coexp_gene.230gene.txt",header=T,row.names = NULL)
#key23gene<-read.delim("C:\\项目\\TCGA_PAAD\\23gene.list.txt",header=T,row.names = NULL)
deg.exp<-exp.tmp[,key358gene$Symbol]

sur<-read.delim("TCGA-PAAD_clinical.survival.haveRNA-Seq.txt",header=T,row.names = 1)
sur.1<-sur[,c("submitter_id","OS","OS.time")]
rownames(sur.1)<-sur.1$submitter_id
sur.1[1:4,1:3]
deg.exp<-deg.exp[rownames(sur.1),]
deg.exp[1:4,1:3]
merge.data<-cbind(sur.1,deg.exp)
library("dplyr")
data.final<-dplyr::select(merge.data,-submitter_id)
data.final[1:4,1:3];dim(data.final)
data.final.358<-data.final
write.table(data.final.358, file = "PAAD.data.for.lasso.358gene.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
head(data.final)
dim(data.final)

#data.final=read.table("PAAD.data.for.lasso.358gene.txt",header=T, row.names=1)
head(data.final)
dim(data.final)
data.final[1:4,1:5]

tmp<-read.delim("survival.p.137gene.res.cat.txt",header=T, row.names=1)
dim(tmp)
head(tmp)
data.137<-data.final[,c("OS","OS.time",rownames(tmp))]
tmp2<-gsub("-",".",rownames(tmp))
data.137<-data.final[,c("OS","OS.time",tmp2)]
data.23<-data.final[,c("OS","OS.time",key23gene$gene23)]
#############################
##Lasso回归,教程https://mp.weixin.qq.com/s/EpVcu7x-EI2bzGjH5lG19Q
{
  library("glmnet")
  library("survival")
  #读取数据
  data_lasso <- data.23 
  #生存时间从day转化为year，除365
  data_lasso$OS.time=data_lasso$OS.time/365   
  #设置随机种子
  set.seed(10)   
  #设置lasso自变量和因变量
  x=as.matrix(data_lasso[,c(3:ncol(data_lasso))])
  y=data.matrix(Surv(data_lasso$OS.time,data_lasso$OS))
  #跑lasso，并绘图
  model <- glmnet(x, y, family = "cox",nlambda = 100,alpha = 1, maxit = 10000)
  plot(model)
  plot(model,xvar="lambda",label=T)
  plot(model, xvar = "lambda", label = TRUE)
  cvmodel <-  cv.glmnet(x, y, family="cox", nlambda = 100,maxit = 10000)#默认nfolds = 10
  plot(cvmodel)
}
####经过10折交叉验证。选择了两个λ（对应于图上的两条虚线）。通常左侧虚线对应的选择的变量最优。即选择了24个变量。
####################################
#最终模型系数估计
# 最终模型的系数估计
#find coefficients of best model
best_model <- glmnet(x, y, alpha = 1, family = "cox",lambda = cvmodel$lambda.min)
coef(best_model)

# 输出相关系数
coef <- coef(model, s = cvmodel$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]

geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
geneCoef     #查看相关系数最终使用的
#Gene     Coef                  
#[1,] "ECT2"   "0.00154974662234233" 
#[2,] "FBXL16" "-0.00142511256816181"
#[3,] "KCNB1"  "-0.00585712367493594"
#[4,] "MET"    "0.00299166561721053" 
#geneCoef<-as.data.frame(geneCoef)
#geneCoef<-subset(geneCoef,Gene %in% c(gene6))
save(geneCoef,file="geneCoef.3gene.final.RData")
save(actCoef,file="actCoef.3gene.final.RData")

head(data_lasso)
test_data<-data.final[,c("OS","OS.time",lassoGene)]
##逐步回归法构建最优模型
tpm<-read.delim("tdmultiCox.67gene.txt",header=T, row.names=1)
test_data<-data.final[,c("OS","OS.time",rownames(tpm))]
install.packages("My.stepwise")
library(My.stepwise)
v1<-colnames(test_data)[c(3:ncol(test_data))]
formula=Surv(OS.time, OS)
test_data<-na.omit(test_data)
My.stepwise.coxph(Time = "OS.time",
                  Status = "OS",
                  variable.list = v1,
                  data=test_data)
# 基于上述相关系数计算风险评分，并将评分中的高于中位数的定义为高风险，低于中位数的定义为低风险
#filter$gene.107
FinalGeneExp <- data_lasso[,lassoGene]
#FinalGeneExp <- data_lasso[,gene6]
#myFun <-  function(x){crossprod(as.numeric(x),as.numeric(geneCoef$Coef))}
myFun <-  function(x){crossprod(as.numeric(x),actCoef)}
riskscore <-  apply(FinalGeneExp,1,myFun)
outCol <-  c("OS.time", "OS", lassoGene)
risk <-  as.vector(ifelse(riskscore > median(riskscore), "high", "low"))#定义低、高风险
dat <- cbind(data_lasso[,outCol], riskscore=as.vector(riskscore), risk) 
dat.4gene.final<-dat
save(dat.4gene.final,file="dat.4gene.final.RData")
#data.final[1:4,285:288]
#dat<-data.final
####输出选择的变量以及风险评分
#基于上述风险评估得分riskScore（连续变量）和癌症患者的结局OS（二分类），可以绘制箱线图和ROC曲线，来评估Lasso选择出来的变量进行风险预测的效果如何。
{
  library(ggpubr) 
  ##箱线图
  #绘制x为甲状腺癌结局，y为风险评估得分的箱图，并进行统计学检验
  m <- ggboxplot(dat, x = "OS", y = "riskscore",
                 color = "OS", palette = "jco",
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "riskscore.OS.boxplot.pdf", width = 2.8, height = 3.2)
  library(ROCR)   
  library(glmnet)
  library(caret)
  ##ROC曲线
  pred <- prediction(dat$riskscore, dat$OS)
  ROC <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values #查看AUC
  #[1] 0.9269531
  plot(ROC,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  #y轴设置为灵敏度，x轴设置为1-特异度
    lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  #绘制对角线
    legend(0.4,0.8,'AUC=0.729') #添加AUC值
}

###########################################################
######单因素cox比例风险模型生存分析，教程https://mp.weixin.qq.com/s?__biz=MzA5Nzg5MzI1Ng==&mid=2649204631&idx=1&sn=45606f39e7a497601abe50c430791939&chksm=888ae970bffd60667d10e3993259c842bceeee9b8f9f63213fdbb9887f3051a9d453283b3af5#rd
data.final.358[1:4,1:5];dim(data.final.358)
dim(data.final.358)
data.final.132[1:3,1:5]
dim(data.final.132)
#OS OS.time  ABHD17C CGB3     SI
#TCGA-HV-A7OP  0     978  39.3106    0 0.0000
#TCGA-F2-A8YN  0     517 106.0922    0 0.0781
#TCGA-3A-A9IJ  0    1854  16.1221    0 0.0087
#TCGA-HZ-7918  0     969  40.3007    0 0.2050
#1.1 分析每个基因表达水平与样品所属个体生存预后的相关性。

pFilter=0.05
outResult=data.frame()
sigGenes=c("OS","OS.time")
#LCP1+NYAP1+ANLN+CD36+CELSR3+CHP2+ECT2+HOXA5+IL20RB+MUC16+VSIG1+PTGES
#gene.tmp<-c("LCP1","NYAP1","ANLN","CD36","CELSR3",
 #           "CHP2","ECT2","HOXA5","IL20RB","MUC16","VSIG1","PTGES")
#data.tmp<-data.final[,c("OS","OS.time",gene.tmp)]
for(i in colnames(data.final.358[,3:ncol(data.final.358)])){
  tdcox <- coxph(Surv(OS.time, OS) ~ data.final.358[,i], data = data.final.358)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
#1.2 得到outResult数据
head(outResult);dim(outResult)
write.table(outResult,file="UniCoxSurvival.105.txt",sep="\t",row.names=F,quote=F)
#1.3 提取有统计学意义的基因的表达数据，并输出以备后续使用
td<-data.final
UniCoxSurSigGeneExp=td[,sigGenes]
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)
#后面进行多因素分析用或者lasso
#lasso
data.final.358
data.final.107<-data.final.358[,c("OS","OS.time",outResult$id)]
##Lasso回归,教程https://mp.weixin.qq.com/s/EpVcu7x-EI2bzGjH5lG19Q
{
  library("glmnet")
  library("survival")
  #读取数据
  data_lasso <- data.final.107 
  #生存时间从day转化为year，除365
  data_lasso$OS.time=data_lasso$OS.time/365   
  #设置随机种子
  set.seed(10)   
  #设置lasso自变量和因变量
  x=as.matrix(data_lasso[,c(3:ncol(data_lasso))])
  y=data.matrix(Surv(data_lasso$OS.time,data_lasso$OS))
  #跑lasso，并绘图
  model <- glmnet(x, y, family = "cox", nlambda = 1000,maxit = 10000)
  plot(model)
  plot(model,xvar="lambda",label=T)
  plot(model, xvar = "lambda", label = TRUE)
  cvmodel <-  cv.glmnet(x, y, family="cox", nlambda = 1000,maxit = 10000)#默认nfolds = 10
  plot(cvmodel)
}
####经过10折交叉验证。选择了两个λ（对应于图上的两条虚线）。通常左侧虚线对应的选择的变量最优。即选择了24个变量。
####################################

# 输出相关系数
coef <- coef(model, s = cvmodel$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]

geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
geneCoef     #查看相关系数
#Gene      Coef                  
#[1,] "NT5E"    "0.000165375463850836"
#[2,] "CCDC74B" "-0.00880595075644125"
#[3,] "KCNB1"   "-0.01363625465643"   
#[4,] "FBXL16"  "-0.00254405013294404"
#[5,] "ECT2"    "0.0048462765089848"  
#[6,] "MET"     "0.00348054572524209" 
#geneCoef<-as.data.frame(geneCoef)
#geneCoef<-subset(geneCoef,Gene %in% c(gene6))

# 基于上述相关系数计算风险评分，并将评分中的高于中位数的定义为高风险，低于中位数的定义为低风险
#filter$gene.107
FinalGeneExp <- data_lasso[,lassoGene]
#FinalGeneExp <- data_lasso[,gene6]
#myFun <-  function(x){crossprod(as.numeric(x),as.numeric(geneCoef$Coef))}
myFun <-  function(x){crossprod(as.numeric(x),actCoef)}
save(actCoef,file="actCoef.3gene.final.RData")
riskScore <-  apply(FinalGeneExp,1,myFun)
outCol <-  c("OS.time", "OS", lassoGene)
risk <-  as.vector(ifelse(riskScore > median(riskScore), "high", "low"))#定义低、高风险
dat <- cbind(data_lasso[,outCol], riskScore=as.vector(riskScore), risk) 
dat.3gene.final<-dat
save(dat.3gene.final,file="dat.3gene.final.Rdata")
#data.final[1:4,285:288]
#dat<-data.final
####输出选择的变量以及风险评分
#基于上述风险评估得分riskScore（连续变量）和癌症患者的结局OS（二分类），可以绘制箱线图和ROC曲线，来评估Lasso选择出来的变量进行风险预测的效果如何。
{
  library(ggpubr) 
  ##箱线图
  #绘制x为甲状腺癌结局，y为风险评估得分的箱图，并进行统计学检验
  m <- ggboxplot(dat, x = "OS", y = "riskScore",
                 color = "OS", palette = "jco",
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "riskscore.OS.boxplot.pdf", width = 2.8, height = 3.2)
  library(ROCR)   
  library(glmnet)
  library(caret)
  ##ROC曲线
  pred <- prediction(dat$riskScore, dat$OS)
  ROC <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values #查看AUC
  #[1] 0.9269531
  plot(ROC,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  #y轴设置为灵敏度，x轴设置为1-特异度
    lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  #绘制对角线
    legend(0.4,0.8,'AUC=0.726') #添加AUC值
}



#################################################
##1 多因素cox分析

library(survival)
td=read.table("UniCoxSurSigGeneExp.txt",header=T, row.names=1)
#1.1 多因素不需要for循环，直接开始分析：

td.test<-data.final[,c("OS","OS.time",c(outResult$id))]
tdmultiCox=coxph(Surv(OS.time, OS) ~ ., data = data.final.107) #这里有个“.”，代表分析td数据中所有变量（列名）
write.table(summary(tdmultiCox)$coefficients,file="tdmultiCox.67gene.txt",sep="\t",row.names=T,quote=F)

tdmultiCox=coxph(Surv(OS.time, OS) ~ NYAP1+CD36+IL20RB, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ PPP1R14D+ZNF488+TMEM213+ANLN+TTC6+MFSD2A+ANKRD22+GPR78, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ PPP1R14D+ZNF488+TMEM213+ANLN+TTC6+MFSD2A+ANKRD22+GPR78, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ PIK3R5+CHI3L2+CLEC12B+CXCL6+CELSR3+FOXD1+NR1D1+PTGES+RAB26+VSIG1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CELSR3+IL20RB+PTGES+VSIG1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ KCNB1+ECT2+MET, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+PRKAR2B+VSIG1, data = data.final)


tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+VSIG1+GPR78+TTC6, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+VSIG1, data = data.final)####final
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CELSR3+IL20RB+VSIG1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+VSIG1+PTGES, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ LCP1+NYAP1+ANLN+CD36+CELSR3+CHP2+ECT2+HOXA5+IL20RB+MUC16+VSIG1+PTGES, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CELSR3+IL20RB+VSIG1+PTGES, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ANLN+CD36+CELSR3+IL20RB+VSIG1, data = data.final)####final
tdmultiCox=coxph(Surv(OS.time, OS) ~ ECT2+MET, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ NT5E+NIM1K+FBXL16+ECT2+MET, data = data.final)#NIM1Kbubiaoda
tdmultiCox=coxph(Surv(OS.time, OS) ~ NT5E+NIM1K+FBXL16+AGBL4+ECT2+MET, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ECT2+FBXL16+KCNB1+MET, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ASCL1+ATP1A3+DSCAM+GAD1+GSDMC+LDLRAD1+LRRC53+NPAS4+PNMA3+SLC25A47+AGER+ZNF114, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ASCL1+ATP1A3+DSCAM+GAD1+LRRC53+NPAS4+SLC25A47, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ATP1A3+LRRC53+NPAS4+SLC25A47, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ KCNB1+ECT2+MET, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ AGBL4+FBXL16+KCNB1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ CHRNB2+DEPDC1+GBX1+MPPED1+SLC25A47+STXBP5L+UNC80+LRRC53+AMER2+BRINP1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ GBX1+MPPED1+LRRC53+AMER2+BRINP1+STXBP5L+MYT1, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ AGBL4+FGF9+GPR6+NPAS4+ASIC4, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ AGBL4+MPPED1+SLC25A47+FGF9+DEPDC1+NPAS4+RAD54L+SOAT2+TRIM46+DSCAM, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ AGBL4+FBXL16+KCNB1+ECT2+MET, data = data.final)



tdmultiCox=coxph(Surv(OS.time, OS) ~ AGBL4+ECT2+FBXL16+KCNB1+MET+NT5E, data = data.final)#试一试
tdmultiCox=coxph(Surv(OS.time, OS) ~ GSDMC+LDLRAD1+GAD1+RAD54L, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ GSDMC+GAD1+RAD54L, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ ASCL1+ATP1A3+DSCAM+KCNB1+NPAS4+SLC25A47, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ASCL1+ATP1A3+DSCAM+KCNB1+NPAS4, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ECT2+FBXL16+KCNB1+MET, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ BIRC5+DEPDC1+E2F1+FERMT1+RAD54L+SMOC1, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ BIRC5+DEPDC1+E2F1+FERMT1+RAD54L+SMOC1, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+FERMT1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+FERMT1+KIF2C+S100A2+TGM2, data = data.final)#buxing
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+FERMT1+KIF2C+LAMC2+TGM2, data = data.final)#和7基因一样
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+FERMT1+KIF2C+LAMC2+S100A2, data = data.final)#和7个一样
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+FERMT1+LAMC2+S100A2+TGM2, data = data.final)#和7个一样
tdmultiCox=coxph(Surv(OS.time, OS) ~ E2F1+FERMT1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)#不行
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+FERMT1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)#可行final
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+E2F1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)#不行

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+FERMT1+KIF2C+AGBL4+ECT2+KCNB1+MET, data = data.final)#无法验证

tdmultiCox=coxph(Surv(OS.time, OS) ~ E2F1+FERMT1+LAMC2+S100A2+TGM2, data = data.final)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ASCL1+F2+FGF9, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ PCDH15+F2+PNMA3+DSCAM+SOAT2+GPR6, data = data.final)#无法验证
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+FERMT1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)#可行final

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+BRINP1+FERMT1+MET+FOXA2+DYNC1I1+MUC16+FAM131C+AGER, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+BRINP1+AGER, data = data.final)#AUC:0.713

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+FERMT1, data = data.final)#AUC:0.730
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER, data = data.final)#AUC:0.711

setwd("C:\\项目\\TCGA_PAAD")
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+DEPDC1, data = data.final)#AUC:0.717
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+RAD54L, data = data.final)#AUC:0.713
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+BIRC5, data = data.final)#AUC:0.707
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+MET, data = data.final)#AUC:0.750
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+DKK1, data = data.final)#AUC:0.711
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+FOXA2, data = data.final)#AUC:0.711
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+EPHA2, data = data.final)#AUC:0.726
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+DYNC1I1, data = data.final)#AUC:0.709
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+FOSL1, data = data.final)#AUC:0.701
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+MUC16, data = data.final)#AUC:0.701
tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+L1CAM, data = data.final)#AUC:0.701

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+MUC16, data = data.final)#AUC:0.735 final.new

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1B+HES2+AGER+MUC16+DYNC1I1, data = data.final)

tdmultiCox=coxph(Surv(OS.time, OS) ~ DEPDC1+FERMT1+KIF2C+LAMC2+S100A2+TGM2, data = data.final)#可行final

ph_hypo_multi <- cox.zph(tdmultiCox)
ph_hypo_multi
#write.table(ph_hypo_multi$table,file="ph_hypo_multi.67gene.txt",sep="\t",row.names=T,quote=F)
tdmultiCox
#森林图
library(survminer)
ggforest(tdmultiCox,data = data.final) #画个森林图看看
p<-ggforest(tdmultiCox,data = data.final)
ggsave(filename = "4gene.ggforest.pdf", 
       plot =p,width = 14, height =6.5, units = 'cm')
#计算模型风险评分

library(dplyr)
save(tdmultiCox,file="4gene.final.tdmultiCox.RData")

genes = names(tdmultiCox$coefficients);length(genes)
fp = apply(data.final[,genes], 1,function(k)sum(tdmultiCox$coefficients * k)) #lasso回归模型的预测值就是线性加乘
#fp = data.final[,c("MET")]*0.007541
#fp=apply(data.final[,genes], 1,function(k)sum(actCoef * k))
data.final$riskscore = fp
write.table(data.final,file="data.final.addriskscore.txt",sep="\t",row.names=F,quote=F)

#data.final[1:4,285:288]
dat<-data.final
dat.7gene.final<-dat[,c("OS","OS.time","riskscore")]
####输出选择的变量以及风险评分
#基于上述风险评估得分riskScore（连续变量）和癌症患者的结局OS（二分类），可以绘制箱线图和ROC曲线，来评估Lasso选择出来的变量进行风险预测的效果如何。

  library(ggpubr) 
  ##箱线图
  #绘制x为甲状腺癌结局，y为风险评估得分的箱图，并进行统计学检验
  m <- ggboxplot(dat, x = "OS", y = "riskscore",
                 color = "OS", palette = "jco",
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m 
  ggsave(filename = "boxplot.riskscore.pdf", 
         plot =m,width = 5.2, height =6.5, units = 'cm')
  
  library(ROCR)   
  library(glmnet)
  library(caret)
  ##ROC曲线
  pred <- prediction(dat$riskscore, dat$OS)
  ROC <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values #查看AUC
  #[1] 0.9269531
  plot(ROC,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  #y轴设置为灵敏度，x轴设置为1-特异度
    lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  #绘制对角线
    legend(0.4,0.8,'AUC=0.735') #添加AUC值



tdmultiCoxSum=summary(tdmultiCox)
outResult.m=data.frame()
outResult.m=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
head(outResult.m)
dim(outResult.m)
#1.2 多因素分析显示只有gene39和gene87的p值小于0.05，输出结果：
write.table(outResult,file="multiCox.67gene.txt",sep="\t",row.names=F,quote=F)

p.data<-read.delim("27gene.UniCoxSurvival.txt",header=T,row.names = 1)
p.data
filter<-p.data[which(p.data$pvalue <0.05),]#26gene
dim(filter)

gene<-rownames(filter)
gene#9gene
#[1] "PPP1R14D" "ZNF488"   "VSTM2L"   "TMEM213"  "ANLN"     "TTC6"     "MFSD2A"   "ANKRD22" 
#[9] "GPR78"

finalgene<-intersect(lassoGene,gene)#9 gene

#######33
###2 选取某几个变量进行多因素cox分析
##2.1 选择所有变量用.表示，若选取其中两个变量，如下命令：
tdmultiCox=coxph(Surv(surtime, surstat) ~ (gene39+gene87), data = td) # 价格（）里面写上列名，可以继续+



######################################################################################################################################
####################################################################################################
#####免疫基因和差异表达基因求取交集得到107个基因
gene.list<-read.delim("C:\\项目\\TCGA_PAAD\\免疫抑制和促进因子\\union.immport.InnateDB.immune.gene.txt",header=T,row.names = 1)
deg.exp[1:3,1:4]
deg.exp.list<-colnames(deg.exp)

length(unique(intersect(gene.list$Symbol,deg.exp.list)))
gene.107<-unique(intersect(gene.list$Symbol,deg.exp.list))
write.table(gene.107,file="C:/项目/TCGA_PAAD/immune.deg.commongene.107.txt",sep="\t",row.names=F,quote=F)
#[1] 107
#######将107个基因与与OS密切相关的基因再一次求取交集
####survival_analysis
#生存分析
data.final=read.table("PAAD.data.for.lasso.txt",header=T, row.names=1)
head(data.final)
dim(data.final)
dat.sur<-data.final
dat.sur$OS.time.year=dat.sur$OS.time/365
###五年生存状态
dat.sur$OS.5<-dat.sur$OS
dat.sur$OS.5[(dat.sur$OS.time.year > 5)] <-  0

library(survminer)
library(survival)

gene.107
gene6
gene3
gene <- c("ITGB4")#挑选进行生存分析的目标基因
#在临床矩阵中新增gene列，根据中位数将表达量（连续型变量）分为高表达和低表达两组（分类型变量）：
pdf(file=paste(gene6[r],"Survival.curves.pdf",sep='.'),6,6)

pdf(file="Survival.curves.3gene.pdf",4.5,5.8)
for(r in 1:length(gene3)){
  gene<-gene3[r]
  dat.sur$gene <- ifelse(dat.sur[,gene] > median(dat.sur[,gene]),"High","Low")
 
#dat.sur$gene[1:6]
  dat.sur$gene<-factor(dat.sur$gene,levels=c("High","Low"))
  fit2 <- survfit(Surv(OS.time.year, OS) ~ gene, data = dat.sur)
  
 # palette = "lancet"
p<-ggsurvplot(
  fit2,
  data = dat.sur,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  #surv.median.line = "hv",
  ggtheme =theme_classic(),
  legend = "top",
  xlab = "OS_time(years)",
  ylab = "Survival probablity",
  title = paste(gene3[r],"Survival curves",sep=' '),
  break.x.by = 1,
  break.y.by = 0.2,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.2,
  risk.table.y.text = FALSE
)
print(p)
}
dev.off()

p_signif<-function(p){
  if(p>0.05) label=c("")
  if(p>0.05&&p<=0.1) label=c(".")
  if(p<=0.05&&p>0.01) label=c("*")
  if(p<=0.01&&p>0.001) label=c("**")
  if(p<=0.001&&p>0.0001) label=c("***")
  if(p<=0.0001&&p>0.00001) label=c("****")
  if(p<=0.00001) label=c("*****")
  return(label)
}
p_signif(as.numeric(Result[i, c("KMp")]))

#####批量输出生存分析p值，筛选显著性基因
gene.107
gene.137<-d
length(gene.137)
gene <- c("ABO")#挑选进行生存分析的目标基因
#在临床矩阵中新增gene列，根据中位数将表达量（连续型变量）分为高表达和低表达两组（分类型变量）：
p<-c()
for(r in 1:length(gene.137)){
  gene<-gene.137[r]
  dat.sur$gene <- ifelse(dat.sur[,gene] > median(dat.sur[,gene]),"High","Low")
  model1 <- survdiff(Surv(OS.time.year, OS) ~ gene,data=dat.sur,na.action = na.exclude)
  p[r]<-model1$pvalue}
p<-as.matrix(p)
rownames(p)<-gene.137
colnames(p)<-c("pvalue")
p.final<-cbind(gene.137,p)
write.table(p.final,file="survival.p.137gene.txt",quote = F,sep="\t",row.names=T)
p.data<-read.delim("survival.p.137gene.txt",header=T,row.names = 1)
p.data
filter<-p.data[which(p.data$pvalue <0.05),]#26gene
dim(filter)

p<-c()
for(r in 1:length(gene.137)){
  gene<-gene.137[r]
  dat.sur$gene <- ifelse(dat.sur[,gene] > median(dat.sur[,gene]),"High","Low")
  model1 <- survdiff(Surv(OS.time.year, OS.5) ~ gene,data=dat.sur,na.action = na.exclude)
  p[r]<-model1$pvalue}
p<-as.matrix(p)
rownames(p)<-gene.137
colnames(p)<-c("pvalue")
p.final<-cbind(gene.137,p)
write.table(p.final,file="survival.p.137gene.5year.txt",quote = F,sep="\t",row.names=T)
p.data<-read.delim("survival.p.137gene.5year.txt",header=T,row.names = 1)
p.data
filter<-p.data[which(p.data$pvalue <0.05),]#26gene
dim(filter)

######单因素cox比例风险模型生存分析，教程https://mp.weixin.qq.com/s?__biz=MzA5Nzg5MzI1Ng==&mid=2649204631&idx=1&sn=45606f39e7a497601abe50c430791939&chksm=888ae970bffd60667d10e3993259c842bceeee9b8f9f63213fdbb9887f3051a9d453283b3af5#rd
data.final[1:4,1:5]
dim(data.final)
#OS OS.time  ABHD17C CGB3     SI
#TCGA-HV-A7OP  0     978  39.3106    0 0.0000
#TCGA-F2-A8YN  0     517 106.0922    0 0.0781
#TCGA-3A-A9IJ  0    1854  16.1221    0 0.0087
#TCGA-HZ-7918  0     969  40.3007    0 0.2050
#1.1 分析每个基因表达水平与样品所属个体生存预后的相关性。
td=read.table("PAAD.data.for.lasso.txt",header=T, row.names=1)
col.new<-gsub("\\.","-",colnames(td))
col.new<-gsub("OS-time","OS.time",col.new)
colnames(td)<-col.new
dim(data.final)

dim(td);td[1:4,1:5]
filter$gene.137
td.1<-td[,c("OS","OS.time",filter$gene.137)]

pFilter=0.05
outResult=data.frame()
sigGenes=c("OS","OS.time")
for(i in colnames(td.1[,3:ncol(td.1)])){
  tdcox <- coxph(Surv(OS.time, OS) ~ td.1[,i], data = td.1)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
#1.2 得到outResult数据
head(outResult);dim(outResult)
write.table(outResult,file="27gene.UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
#1.3 提取有统计学意义的基因的表达数据，并输出以备后续使用
td<-data.final
UniCoxSurSigGeneExp=td[,sigGenes]
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)
#后面进行多因素分析用
#########################多因素cox
##1 多因素cox分析

library(survival)
td=read.table("UniCoxSurSigGeneExp.txt",header=T, row.names=1)
#1.1 多因素不需要for循环，直接开始分析：

td=read.table("PAAD.data.for.lasso.txt",header=T, row.names=1)
dim(td);td[1:4,1:5]
filter$gene.137
td.1<-td[,c("OS","OS.time",filter$gene.137)]

tdmultiCox=coxph(Surv(OS.time, OS) ~ ., data = td.1) #这里有个“.”，代表分析td数据中所有变量（列名）
tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
head(outResult);dim(outResult)
#1.2 多因素分析显示只有gene39和gene87的p值小于0.05，输出结果：
write.table(outResult,file="multiCox.26gene.txt",sep="\t",row.names=F,quote=F)
tmp<-read.delim("multiCox.26gene.txt",header=T,row.names = 1)
multiCox.gene<-tmp[which(tmp$pvalue <0.05),]
multiCox.gene
gene10<-rownames(multiCox.gene)
gene10 #3gene
#[1] "NYAP1"  "CD36"   "IL20RB"
#[1] "PIK3R5"  "CHI3L2"  "CLEC12B" "CXCL6"   "CELSR3"  "FOXD1"   "NR1D1"   "PTGES"   "RAB26"   "VSIG1"
td.2<-td[,c("OS","OS.time",gene10)]

multiCox.gene<-read.delim("multiCox.p0.05.txt",header=T,row.names = 1,fileEncoding="UTF16LE")
gene<-rownames(multiCox.gene)
gene6<-intersect(gene,filter$gene.137)
#[1] "ITGB4"   "CD36"    "S100A2"  "ANLN"    "PTGES"   "PRKAR2B"

###2 选取某几个变量进行多因素cox分析
tdmultiCox=coxph(Surv(OS.time, OS) ~ (ITGB4+CD36+S100A2+ANLN+PTGES+PRKAR2B), data = td)
tdmultiCox=coxph(Surv(OS.time, OS) ~ (NYAP1+CD36+IL20RB), data = td.1)
tdmultiCox=coxph(Surv(OS.time, OS) ~ (PIK3R5+CHI3L2+CLEC12B+CXCL6+CELSR3+FOXD1+NR1D1+PTGES+RAB26+VSIG1), data = td.1)
tdmultiCox=coxph(Surv(OS.time, OS) ~ (CHI3L2+CLEC12B+CELSR3), data = td.1)
tdmultiCox=coxph(Surv(OS.time, OS) ~ ., data = td.2)
tdmultiCox=coxph(Surv(OS.time, OS) ~ (CHI3L2+CELSR3+PTGES+VSIG1), data = td.1)
tdmultiCox=coxph(Surv(OS.time, OS) ~ (CHI3L2+CXCL6+VSIG1), data = td.1)
tdmultiCox
#coef  exp(coef)   se(coef)      z       p
#ITGB4   -0.0004139  0.9995862  0.0005409 -0.765 0.44414
#CD36    -0.0178038  0.9823538  0.0144897 -1.229 0.21918
#S100A2   0.0002405  1.0002405  0.0002303  1.044 0.29647
#ANLN     0.0214040  1.0216347  0.0069040  3.100 0.00193
#PTGES    0.0012399  1.0012407  0.0009322  1.330 0.18348
#PRKAR2B -0.0009722  0.9990283  0.0167211 -0.058 0.95364

tdmultiCoxSum=summary(tdmultiCox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
#id        HR                  L95CI               H95CIH             pvalue               
#ITGB4   "ITGB4"   "0.99958616608517"  "0.998526997249891" "1.00064645841398" "0.444137837995851"  
#CD36    "CD36"    "0.982353781444547" "0.954848019973551" "1.01065188567405" "0.219177728646061"  
#S100A2  "S100A2"  "1.00024048964967"  "0.999789065973673" "1.000692117152"   "0.296469993398253"  
#ANLN    "ANLN"    "1.02163472432883"  "1.00790350888391"  "1.03555300755944" "0.00193365469433555"
#PTGES   "PTGES"   "1.00124066850734"  "0.999413047347545" "1.00307163182792" "0.183479354200319"  
#PRKAR2B "PRKAR2B" "0.999028260587396" "0.966817951761716" "1.03231168146355" "0.953635005332226"

tf<-read.delim("Homo_sapiens_TF.txt",header=T,sep='\t')
intersect(comm_deg,unique(tf$Symbol))
