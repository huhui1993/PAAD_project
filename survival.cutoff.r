
setwd("D:/项目/TCGA_PAAD/GEO")
####R做生存分析如何取最佳cutoff(截断)
library(survival)
library(survminer)
####计算最佳截点

dat.sur=read.table("dat.sur.357.txt",header=T, row.names=1)
head(dat.sur)


dat.sur[1:3,1:5]
#OS OS.time  ABHD17C CGB3     SI
#TCGA-HV-A7OP  0     978  39.3106    0 0.0000
#TCGA-F2-A8YN  0     517 106.0922    0 0.0781
#TCGA-3A-A9IJ  0    1854  16.1221    0 0.0087
# 1. Determine the optimal cutpoint of variables
gene7<-c("DEPDC1","E2F1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")
gene4<-c("DEPDC1B","HES2","AGER","MUC16")

gene137<-read.delim("immuneAndimmune_coexp_gene.137gene.txt",header=T,row.names = 1)
gene137$x

res.cut <- surv_cutpoint(dat.sur, #数据集
                         time = "OS.time.year", #生存时间
                         event = "OS", #生存状态
                         variables = c(gene4), #需要计算的数据列名
                         minprop=0.45
)

summary(res.cut) #查看数据最佳截断点及统计量

#cutpoint statistic
#NPL         3.9057 3.3270646
#WDFY4       0.5993 1.8141447
##展示数据分布
# 2. Plot cutpoint for NPL
# 以NPL为例
plot(res.cut, "DEPDC1", palette = "npg")

# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
write.table(res.cat,file="gene4.res.cat.txt",sep="\t",row.names=T,quote=F)

fit <- survfit(Surv(OS.time.year, OS) ~DEPDC1B, data = res.cat)#拟合生存分析
#绘制生存曲线并显示P值
ggsurvplot(fit,
           data = res.cat,
           risk.table = TRUE,
           pval = T,
           conf.int = TRUE,
           conf.int.style = "ribbon",
           conf.int.alpha = 0.2)
###批量绘制生存分析图
gene7<-c("DEPDC1","E2F1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")
gene4<-c("DEPDC1B","HES2","AGER","MUC16")
pdf(file="Survival.curves.4gene.new.pdf",4.8,5.8)
for(r in 1:length(gene4)){
  gene<-gene4[r]
  #dat.sur$gene<-factor(dat.sur$gene,levels=c("High","Low"))
  #fit <- survfit(Surv(OS.time.year, OS) ~res.cat[,gene], data = res.cat)
  fit <- survfit(Surv(OS.time.year, OS) ~get(gene), data = res.cat)
  p<-ggsurvplot(fit,
             data = res.cat,
             risk.table = TRUE,
             pval = T,
             conf.int = TRUE,
             conf.int.style = "ribbon",
             conf.int.alpha = 0.2,
             xlab = "OS_time(years)",
             ylab = "Survival probablity",
             title = paste(gene4[r],"Survival curves",sep=' '))
  print(p)
}
dev.off()

identical(rownames(res.cat),rownames(data.final))
group <-  as.vector(ifelse(data.final$riskscore > median(data.final$riskscore), "high", "low"))#定义低、高风险
dat.tmp<-cbind(res.cat[,1:2],group)

pdf(file="Survival.curves.riskscore.pdf",4.8,5.6)
fit <- survfit(Surv(OS.time.year, OS) ~group, data = dat.tmp)
p<-ggsurvplot(fit,
              data = dat.tmp,
              risk.table = TRUE,
              pval = T,
              conf.int = TRUE,
              conf.int.style = "ribbon",
              conf.int.alpha = 0.2,
              xlab = "OS_time(years)",
              ylab = "Survival probablity",
              title = paste("riskscore","Survival curves",sep=' '))
print(p)
dev.off()

################p
gene.137<-gene137$x
p<-c()
for(r in 1:length(gene.137)){
  #gene<-gene.137[r]
  group<-res.cat[,gene.137[r]]
  model1 <- survdiff(Surv(OS.time.year, OS) ~ group,data=res.cat,na.action = na.exclude)
  p[r]<-model1$pvalue}
p<-as.matrix(p)
rownames(p)<-gene.137
colnames(p)<-c("pvalue")
p.final<-cbind(gene.137,p)
write.table(p.final,file="survival.p.137gene.res.cat.txt",quote = F,sep="\t",row.names=T)
p.data<-read.delim("survival.p.137gene.res.cat.txt",header=T,row.names = 1)
p.data
filter<-p.data[which(p.data$pvalue <0.05),]
dim(filter)#62gene

#生存分析
data.final=read.table("PAAD.data.for.lasso.txt",header=T, row.names=1)
head(data.final)
col.new<-gsub("\\.","-",colnames(data.final))
col.new<-gsub("OS-time","OS.time",col.new)
colnames(data.final)<-col.new
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
gene5<-c("ANLN","CD36","CELSR3","IL20RB","VSIG1")
pdf(file="Survival.curves.5gene.pdf",4.5,5.8)
for(r in 1:length(gene5)){
  gene<-gene5[r]
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
    title = paste(gene5[r],"Survival curves",sep=' '),
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
