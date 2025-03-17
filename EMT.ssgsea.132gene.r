#############################################################################
######ssgsea
#setwd("C:\\项目\\TCGA_PAAD")
setwd("C:\\项目\\TCGA_PAAD\\GEO")

#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
df<-read.delim("C:\\项目\\TCGA_PAAD\\mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
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
write.table(dat, file = "TCGA_PAAD_EMT.score.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
######################################################################
#####################################################################
setwd("D:\\项目\\TCGA_PAAD\\GEO")
dat=read.delim("TCGA_PAAD_EMT.score.txt",header=T, row.names=1)

list<-substring(rownames(dat),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
rownames(dat)<-list1[,1]
head(dat)

#load("dat.5gene.final.Rdata")
head(dat.7gene.final)

dat<-dat[rownames(dat.7gene.final),]
dat.1<-cbind(dat,dat.7gene.final)

riskscore.group <-  as.vector(ifelse(dat.1$riskscore > median(dat.1$riskscore), "high", "low"))
dat.1 <- cbind(dat.1, Riskgroup=as.vector(riskscore.group))
head(dat.1)
save(dat.1,file="dat.1.Rdata")

####输出选择的变量以及风险评分
setwd("D:\\项目\\TCGA_PAAD\\GEO")
load("dat.1.Rdata")
head(dat.1)
sur<-read.delim("D:\\项目\\TCGA_PAAD\\TCGA-PAAD_185sample_clinical.survival.haveRNA-Seq.txt",
                header=T,row.names = 1)
head(sur)
rownames(sur)<-sur$submitter_id
sur<-sur[rownames(dat.1),]
dat.3<-cbind(dat.1,sur)
head(dat.3)
save(dat.3,file="dat.3.RData")
library(ggpubr) 
##箱线图
dat.3.1<-dat.3 %>% subset(dat.3$ajcc_pathologic_stage!="Stage I")
dat.3.1$ajcc_pathologic_stage<-factor(dat.3.1$ajcc_pathologic_stage,
                                    levels=c("Stage IA","Stage IB","Stage IIA",
                                             "Stage IIB","Stage III","Stage IV"))
m <- ggboxplot(dat.3.1, x = "ajcc_pathologic_stage", y = "riskscore",
               color = "ajcc_pathologic_stage", palette = "jco",
               add = "jitter")
m<-m+stat_compare_means(method = "anova", label.y = 3)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Stage IA")                    # Pairwise comparison against reference
m
ggsave(m, 
       filename = "boxplot.riskscore.stage.pdf", 
       width = 5, height = 3.2) 
####
dat.3.1<-dat.3 %>% subset(dat.3$stage!="NA")
dat.3.1$stage<-factor(dat.3.1$stage, 
                    levels=c("Stage I","Stage II",
                             "Stage III","Stage IV"))
m <- ggboxplot(dat.3.1, x = "stage", y = "riskscore",
               color = "stage", palette = "jco",
               add = "jitter")
m<-m+stat_compare_means(method = "anova", label.y = 3)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Stage I")                    # Pairwise comparison against reference
m
ggsave(m, 
       filename = "boxplot.riskscore.stage.short.pdf", 
       width = 4, height = 3.2) 

####
dat.3.1<-dat.3 %>% subset(dat.3$t_type!="NA")
dat.3.1$t_type<-factor(dat.3.1$t_type, 
                      levels=c("T1-3","T4"))
m <- ggboxplot(dat.3.1, x = "t_type", y = "riskscore",
               color = "t_type", palette = "jco",
               add = "jitter")
m <- m + stat_compare_means(method = "t.test")   #  Add p-value
m  
ggsave(m, 
       filename = "boxplot.riskscore.t_type.pdf", 
       width = 3, height = 3.2) 
####
dat.3.1<-dat.3 %>% subset(dat.3$n_type!="NA")
dat.3.1$n_type<-factor(dat.3.1$n_type, 
                       levels=c("N0","N1"))
m <- ggboxplot(dat.3.1, x = "n_type", y = "riskscore",
               color = "n_type", palette = "jco",
               add = "jitter")
m <- m + stat_compare_means(method = "t.test")   #  Add p-value
m  
ggsave(m, 
       filename = "boxplot.riskscore.n_type.pdf", 
       width = 3, height = 3.2) 
####m_type
dat.3.1<-dat.3 %>% subset(dat.3$m_type!="NA")
dat.3.1$m_type<-factor(dat.3.1$m_type, 
                       levels=c("M0","M1"))
m <- ggboxplot(dat.3.1, x = "m_type", y = "riskscore",
               color = "m_type", palette = "jco",
               add = "jitter")
m <- m + stat_compare_means(method = "t.test")   #  Add p-value
m  
ggsave(m, 
       filename = "boxplot.riskscore.m_type.pdf", 
       width = 3, height = 3.2) 

####gender
dat.3.1<-dat.3 %>% subset(dat.3$gender!="NA")
dat.3.1$gender<-factor(dat.3.1$gender, 
                       levels=c("female","male"))
m <- ggboxplot(dat.3.1, x = "gender", y = "riskscore",
               color = "gender", palette = "jco",
               add = "jitter")
m <- m + stat_compare_means(method = "t.test")   #  Add p-value
m  
ggsave(m, 
       filename = "boxplot.riskscore.gender.pdf", 
       width = 3, height = 3.2) 

dim(dat.3[which(dat.3$age_at_index<=60 & dat.3$Riskgroup=="low"),])
dim(dat.3[which(dat.3$age_at_index<=60 & dat.3$Riskgroup=="high"),])
dim(dat.3[which(dat.3$age_at_index>60 & dat.3$Riskgroup=="low"),])
dim(dat.3[which(dat.3$age_at_index>60 & dat.3$Riskgroup=="high"),])

write.table(dat.3,file="TCGA.dat.3.add.riskscore.txt",sep="\t",row.names=T,quote=F)
sur.test<-read.delim("C:\\项目\\TCGA_PAAD\\TCGA-PAAD_185sample_clinical.txt",
                header=T,row.names = 1)
sur.test<-sur.test[,c("submitter_id","Smoking","alcohol_history",
                      "treatments_pharmaceutical_treatment_or_therapy",
                      "treatments_radiation_treatment_or_therapy")]
rownames(sur.test)<-sur.test$submitter_id
sur.test<-sur.test[rownames(dat.3),-1]
dat.3<-cbind(dat.3,sur.test)
write.table(dat.3,file="TCGA.dat.3.add.riskscore.txt",sep="\t",row.names=T,quote=F)
save(dat.3,file="dat.3.RData")

####fisher test
x<-c(26,63,33,55); dim(x)<-c(2,2)
fisher.test(x)

x<-c(52,37,45,43); dim(x)<-c(2,2)
fisher.test(x)

x<-c(84,3,83,4); dim(x)<-c(2,2)
fisher.test(x)

x<-c(29,60,27,61); dim(x)<-c(2,2)
fisher.test(x)

x<-c(48,36,53,28); dim(x)<-c(2,2)
fisher.test(x)

x<-c(66,18,63,21); dim(x)<-c(2,2)
fisher.test(x)

x<-c(33,48,15,67); dim(x)<-c(2,2)
fisher.test(x)
#基于上述风险评估得分riskScore（连续变量）和癌症患者的结局OS（二分类），可以绘制箱线图和ROC曲线，来评估Lasso选择出来的变量进行风险预测的效果如何。
  library(ggpubr) 
  ##箱线图
  #绘制x为甲状腺癌结局，y为风险评估得分的箱图，并进行统计学检验
  m <- ggboxplot(dat, x = "EMTgroup", y = "EMT_score",
                 color = "EMTgroup", palette = "jco",
                 add = "jitter")

  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "EMTgroup.emtscore.boxplot.pdf", width = 2.8, height = 3.2) 
  
  m <- ggboxplot(dat.1, x = "OS", y = "EMT_score",
                 color = "OS", palette = "jco",
                 add = "jitter")
  m <- m + stat_compare_means()   #  Add p-value
  m  
  ggsave(m, filename = "EMTgroup.OS.boxplot.pdf", width = 2.8, height = 3.2) 

  m <- ggboxplot(dat.1, x = "Riskgroup", y = "EMT_score",
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
  dat.1$OS.time.year=dat.1$OS.time/365
  ###五年生存状态
  dat.1$OS.5<-dat.1$OS
  dat.1$OS.5[(dat.1$OS.time.year > 5)] <-  0
  
  library(survminer)
  library(survival)
  
    #dat.sur$gene[1:6]
    dat.1$EMTgroup<-factor(dat.1$EMTgroup,levels=c("high","low"))
    fit2 <- survfit(Surv(OS.time.year, OS) ~ EMTgroup, data = dat.1)
    
    # palette = "lancet"
    p<-ggsurvplot(
      fit2,
      data = dat.1,
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
      title = "Survival curves",
      break.x.by = 1,
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
  

#
    library(survminer)
    library(survival)
    
    #dat.sur$gene[1:6]
    dat.1$Riskgroup<-factor(dat.1$Riskgroup,levels=c("high","low"))
    fit2 <- survfit(Surv(OS.time.year, OS) ~ Riskgroup, data = dat.1)
    
    # palette = "lancet"
    p<-ggsurvplot(
      fit2,
      data = dat.1,
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
      title = "Survival curves",
      break.x.by = 1,
      break.y.by = 0.2,
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.height = 0.2,
      risk.table.y.text = FALSE
    )
    p
    pdf(file="Survival.curves.Riskgroup.pdf",4.5,5.8)
    p
    dev.off()
#################################################################################
#############gene exp VS riskscore
    setwd("D:\\项目\\TCGA_PAAD\\GEO")
    load("dat.1.Rdata")
    df<-read.delim("D:\\项目\\TCGA_PAAD\\mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
    list<-substring(colnames(df),1,12)
    list1<-as.data.frame(gsub("\\.", "-", list))
    
    exp.tmp<-t(df)
    rownames(exp.tmp)<-list1[,1]
    exp.tmp<-exp.tmp[rownames(dat.1),]
    dat.2<-cbind(dat.1,exp.tmp[,c("DEPDC1B","HES2","AGER","MUC16",
                                  "FOXA2","DKK1","EPHA2","FOSL1","SMOC1","S100A2")])
    
    dat.2<-cbind(dat.1,exp.tmp[,c("DEPDC1B","HES2","AGER","MUC16")])
    
    head(dat.2)
    
    genes<-c("DEPDC1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")
    genes<-c("DEPDC1B","HES2","AGER","FOXA2")
    genes<-c("DEPDC1B","HES2","AGER","MUC16")
    library(ggpubr) 
    ##箱线图
    #绘制x为风险分组，y为基因表达量的箱图，并进行统计学检验
    for(r in 1:length(genes)){
    m <- ggboxplot(dat.2, x = "Riskgroup", y = genes[r],
                   color = "Riskgroup", palette = "jco",
                   add = "jitter")
    
    m <- m + stat_compare_means()   #  Add p-value
    m  
    ggsave(m, filename = paste("boxplot.",genes[r],"exp.in.riskgroup.pdf"), width = 2.8, height = 3.2) 
    }
    #
    library("reshape")
    dat_melt2<-melt(dat.2,id=c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup"))
    colnames(dat_melt2)<-c("EMT_score","EMTgroup","OS","OS.time","riskscore","Riskgroup","type","value")
    dat_melt2$Riskgroup<-factor(dat_melt2$Riskgroup,levels=c("low","high"))
    
    dat_melt2$value<-as.numeric(dat_melt2$value)
    p<-ggboxplot(dat_melt2, x = "Riskgroup", y = "value",color = "Riskgroup",
                 palette = "jco",add = "jitter",
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
      labs(fill="type",x="",y="Gene expression",title="")
    
    my_comparisons <- list( c("high", "low"))
    
    p2<-p+facet_wrap(~ type, scales="free",nrow = 1) + 
      stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                         size=4,method="wilcox.test")+
      theme(strip.text.x = element_text(size=11, color="black",
                                        face="bold"))+
      theme(strip.background = element_blank(),strip.placement = "outside")
    p2
    ggsave(filename="boxplot.4gene.exp.in.riskgroup.pdf",plot=p2,
           device='pdf',path=".",width=7.5,height=3.2)
    
  #
    head(dat_melt2)
    dat_melt2$EMTgroup<-factor(dat_melt2$EMTgroup,levels=c("low","high"))
    
    dat_melt2$value<-as.numeric(dat_melt2$value)
    p<-ggboxplot(dat_melt2, x = "EMTgroup", y = "value",color = "EMTgroup",
                 palette = c("DodgerBlue","OrangeRed"),add = "jitter",
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
      labs(fill="type",x="",y="Gene expression",title="")
    
    my_comparisons <- list( c("high", "low"))
    
    p2<-p+facet_wrap(~ type, scales="free",nrow = 1) + 
      stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                         size=4,method="wilcox.test")+
      theme(strip.text.x = element_text(size=11, color="black",
                                        face="bold"))+
      theme(strip.background = element_blank(),strip.placement = "outside")
    p2
    ggsave(filename="boxplot.4gene.exp.in.EMTgroup.new.pdf",plot=p2,
           device='pdf',path=".",width=6,height=3.5)
    #绘制x为EMT分组，y为基因表达量的箱图，并进行统计学检验
    for(r in 1:length(genes)){
      m <- ggboxplot(dat.2, x = "EMTgroup", y = genes[r],
                     color = "EMTgroup", palette = "jco",
                     add = "jitter")
      
      m <- m + stat_compare_means()   #  Add p-value
      m  
      ggsave(m, filename = paste("boxplot.",genes[r],"exp.in.EMTgroup.pdf"), width = 2.8, height = 3.2) 
    }
    ##相关性点图
    #绘制x为EMT分数，y为基因表达量的点图，并进行统计学检验
    head(dat.2)
    
######################################################################################
#######################################################################################
setwd("C:\\项目\\TCGA_PAAD\\GEO")
load("dat.1.Rdata")
head(dat.1)

####survival_analysis
#生存分析
data.final=read.table("C:/项目/TCGA_PAAD/PAAD.data.for.lasso.358gene.txt",header=T, row.names=1)
head(data.final)
dim(data.final)
dat.sur<-data.final
dat.sur$OS.time.year=dat.sur$OS.time/365
write.table(dat.sur,file="dat.sur.357.txt",sep="\t",row.names=T,quote=F)
###五年生存状态
dat.sur$OS.5<-dat.sur$OS
dat.sur$OS.5[(dat.sur$OS.time.year > 5)] <-  0

library(survminer)
library(survival)

gene <- c("ITGB4")#挑选进行生存分析的目标基因

#在临床矩阵中新增gene列，根据中位数将表达量（连续型变量）分为高表达和低表达两组（分类型变量）：
#pdf(file=paste(gene6[r],"Survival.curves.pdf",sep='.'),6,6)
gene7<-c("DEPDC1","E2F1","FERMT1","KIF2C","LAMC2","S100A2","TGM2")
pdf(file="Survival.curves.7gene.pdf",4.5,5.8)
for(r in 1:length(gene7)){
  gene<-gene7[r]
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
    title = paste(gene7[r],"Survival curves",sep=' '),
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

#########################################################
####批量Cox单因素分析
setwd("C:\\项目\\TCGA_PAAD\\GEO")

dat.3<-read.delim("TCGA.dat.3.add.riskscore.txt",
                     header=T,row.names = 1)
sur.test<-sur.test[,c("submitter_id","Smoking","alcohol_history",
                      "treatments_pharmaceutical_treatment_or_therapy",
                      "treatments_radiation_treatment_or_therapy")]


head(dat.3)
install.packages("ezcox")
library(ezcox)
res <-  ezcox(dat.3, 
              covariates = c("riskscore","ajcc_pathologic_stage","stage", "gender",
                             "age_at_index","Smoking","alcohol_history",
                             "treatments_pharmaceutical_treatment_or_therapy",
                             "treatments_radiation_treatment_or_therapy"), ## 变量
              time = "OS.time",
              status = "OS",
              return_models = TRUE
)
library(forestmodel)
mds <- get_models(res)  ## 提取模型
show_models(mds)  ## 森林图
show_models(mds, merge_models = TRUE, drop_controls = TRUE) ## 森林图

## 直接画森林图

show_forest(dat.3, 
            covariates = c("riskscore","ajcc_pathologic_stage","stage", "gender",
                           "age_at_index","Smoking","alcohol_history",
                           "treatments_pharmaceutical_treatment_or_therapy",
                           "treatments_radiation_treatment_or_therapy"),
            time = "OS.time",
            status = "OS")
dat.4<-dat.3 %>% subset(treatments_pharmaceutical_treatment_or_therapy!="not reported")
dat.4<-dat.4 %>% subset(treatments_radiation_treatment_or_therapy!="not reported")
dim(dat.4)
tdmultiCox=coxph(Surv(OS.time, OS) ~ riskscore+ajcc_pathologic_stage+
                   gender+age_at_index+Smoking+alcohol_history+
                   treatments_pharmaceutical_treatment_or_therapy+
                   treatments_radiation_treatment_or_therapy, data = dat.4)

tdmultiCox=coxph(Surv(OS.time, OS) ~ riskscore+
                   gender+age_at_index+
                   treatments_pharmaceutical_treatment_or_therapy+
                   treatments_radiation_treatment_or_therapy, data = dat.4)

ph_hypo_multi <- cox.zph(tdmultiCox)
ph_hypo_multi
#write.table(ph_hypo_multi$table,file="ph_hypo_multi.67gene.txt",sep="\t",row.names=T,quote=F)
tdmultiCox
#森林图
library(survminer)
ggforest(tdmultiCox,data = dat.4) #画个森林图看看
p<-ggforest(tdmultiCox,data = data.final)
ggsave(filename = "4gene.ggforest.pdf", 
       plot =p,width = 14, height =6.5, units = 'cm')
