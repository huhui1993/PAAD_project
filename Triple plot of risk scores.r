setwd("./GEO")
#相关R包载入：
library(ggrisk)
library(survival)
library(survminer)

#先使用ggrisk包内置数据集进行测试：
dt <- LIRI #LIRI来自日本的肝癌ICGC数据库，包含时间、事件和四个基因
head(dt)

load("4gene.final.tdmultiCox.RData")
tdmultiCox
#coef exp(coef)  se(coef)      z        p
#DEPDC1B  0.049422  1.050663  0.013955  3.542 0.000398
#HES2     0.029827  1.030277  0.012904  2.312 0.020804
#AGER    -0.054931  0.946551  0.025128 -2.186 0.028815
#MUC16    0.014346  1.014449  0.006485  2.212 0.026951

dat.sur=read.table("dat.sur.357.txt",header=T, row.names=1)
head(dat.sur)
dat.test<-dat.sur[,c("OS","OS.time.year","DEPDC1B","HES2","AGER","MUC16")]
head(dat.test)
#使用四个基因构建多因素cox回归模型：
fit <- coxph(Surv(OS.time.year, OS) ~ DEPDC1B + HES2 + AGER + MUC16, dat.test)
fit
#风险因子联动图绘制：
ggrisk(fit)

#颜色相关参数调整：
pdf(file="Triple_plot_of_risk_scores.pdf",4.5,5)
ggrisk(
  fit,
  cutoff.value = "median",
  cutoff.x = NULL,
  cutoff.y = 1,
  cutoff.label = NULL,
  code.highrisk = "High Risk",
  code.lowrisk = "Low Risk",
  title.A.ylab = "Risk Score",
  title.B.ylab = "OS_Time",
  title.A.legend = "Risk Group",
  title.B.legend = "Status",
  title.C.legend = "Expression",
  color.A = c(low = "MediumTurquoise", high = "IndianRed1"), #图A颜色调整
  color.B = c(code.0 = "MediumTurquoise", code.1 = "IndianRed1"), #图B颜色调整
  color.C = c(low = "MediumTurquoise", median = "white", high = "IndianRed1"), #图C颜色调整
)
dev.off()
