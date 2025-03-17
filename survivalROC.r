###教程https://mp.weixin.qq.com/s/h6317qHErWualRI_rM6mLQ
setwd("C:\\项目\\TCGA_PAAD")
write.table(data.final,file="data.final.addriskscore.txt",sep="\t",row.names=F,quote=F)
head(data.final)
#"dat.5gene.final.RData"
dat.sur<-dat.1
dat.sur$OS.time.year=dat.sur$OS.time/365
rt<-dat.sur[,c("OS","OS.time.year","riskscore")]

head(rt)
#futime fustat  riskScore risk        order num
#1 7.1561644      0 -0.7054783  low TCGA-BH-A0BS   1
#2 5.2575342      0 -0.6712008  low TCGA-BH-A0AZ   2


library(tidyverse)
library(survivalROC)
#本次使用“survivalROC”这个包来计算AUC 首先定义一个函数，这个函数可以计算任意时间的ROC
survivalROC_helper <- function(t) {
  
  survivalROC(Stime=rt$OS.time.year, status=rt$OS, marker = rt$riskscore, 
              
              predict.time =t, method="KM")
  
}
###接着计算不同时间的ROC,这里选择的是1，3，5年
library(tidyverse)
library(survivalROC)

survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  
  mutate(survivalROC = map(t, survivalROC_helper),
         
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
###如果需要把三个时间点绘制在一张图上，可以这样做
survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.3f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))

pdf(file="Survival.ROC.135year.new.pdf",4,4)
survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color= AUC))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +
  theme(legend.position = c(0.8,0.2))
dev.off()
#############另外一种方法可以通过ggplot2的分面，把画出三张图