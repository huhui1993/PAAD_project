setwd("K:\\??????Ŀ\\AnimalTFDB??hTFtargetϵͳ????\\TCGA_tumor_sample\\??????֯??TCGA??֢??Ӧ֮??TF????��????")
getwd()

list.files(pattern="*.txt")
###all_gene
data<-read.delim("1556.commonTF.tissue_and_cancer.TPM.sort.final.addGroup.txt",header=T,sep='\t')
head(data)
dim(data)

library("ggplot2")
library("scales")
library("ggpubr")
library("reshape")
#colnames(data)<-c("Num","entrez_id","symbol","Sample","RSEM","cancer_type","type")
head(data)
data$log2_TPM<-log2(data$TPM+1)
head(data)
ggboxplot(data, x="Tissue", y="log2_TPM", fill="Type",
          palette = c("blue","LightSkyBlue"),bxp.errorbar = FALSE, 
          size = 0.5,notch =TRUE,width = 0.7,outlier.colour = NA) +
  labs(y="Expression level (log2TPM)",x="",title="")+theme_bw()+
  theme(axis.text.x = element_text(size = 13, angle =90, hjust = 0,vjust =0,colour = "black",lineheight=100),
        axis.text.y = element_text(size=13,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=13,colour="black"),
        axis.title.y = element_text(size=13,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())



####calculate p value

data.tmp<-data[which(data$Group=="Group4"),]
p<-ggplot(data=data.tmp, aes(x=Tissue,y=log2_TPM))+geom_boxplot(aes(fill=Type),
                                                          outlier.colour = NA,
                                                          show.legend = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  scale_fill_manual(values=co)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p+stat_compare_means(aes(label = paste0("p = ", ..p.format..)),size=4,method="anova",vjust = 0.5,label.x = 1.2)

head(data)
group_list<-unique(data$Group)
list1<-group_list[-c(4,8,10,12)]
list2<-group_list[c(4,8,10,12)]

cancer_list<-list1
rownum<- length(cancer_list)
pvalue<-0
TF_mean<-0
non_TF_mean<-0
for(n in 1:rownum){
  HNSC<-data[which(data$Group==cancer_list[n]),]
  HNSC_non_TF<-HNSC[which(HNSC$Type=="Cancer"),]
  HNSC_TF<-HNSC[which(HNSC$Type=="Tissue"),]
  tmp<-t.test(HNSC_non_TF$TPM,HNSC_TF$TPM, paired=FALSE)
  pvalue[n]<-tmp$p.value
  HNSC_non_TF_mean<-median(HNSC_non_TF$TPM)
  HNSC_TF_mean<-median(HNSC_TF$TPM)
  TF_mean[n]<-HNSC_TF_mean
  non_TF_mean[n]<-HNSC_non_TF_mean
}
pvalue<-as.matrix(pvalue)
TF_mean<-as.matrix(TF_mean)
non_TF_mean<-as.matrix(non_TF_mean)
rownames(pvalue)<-cancer_list
rownames(TF_mean)<-cancer_list
rownames(non_TF_mean)<-cancer_list
pvalue.final<-cbind(TF_mean,non_TF_mean,pvalue)
colnames(pvalue.final)<-c("Tissue","Cancer","pvalue")
write.table(pvalue.final,"t.test.pvalue.normalTssue_map2_cancer.median.txt",quote = F,sep='\t')


cancer_list<-list1
rownum<- length(cancer_list)
pvalue<-0
TF_mean<-0
non_TF_mean<-0
for(n in 1:rownum){
  HNSC<-data[which(data$Group==cancer_list[n]),]
  HNSC_non_TF<-HNSC[which(HNSC$Type=="Cancer"),]
  HNSC_TF<-HNSC[which(HNSC$Type=="Tissue"),]
  tmp<-t.test(HNSC_non_TF$TPM,HNSC_TF$TPM, paired=FALSE)
  pvalue[n]<-tmp$p.value
  HNSC_non_TF_mean<-mean(HNSC_non_TF$TPM)
  HNSC_TF_mean<-mean(HNSC_TF$TPM)
  TF_mean[n]<-HNSC_TF_mean
  non_TF_mean[n]<-HNSC_non_TF_mean
}
pvalue<-as.matrix(pvalue)
TF_mean<-as.matrix(TF_mean)
non_TF_mean<-as.matrix(non_TF_mean)
rownames(pvalue)<-cancer_list
rownames(TF_mean)<-cancer_list
rownames(non_TF_mean)<-cancer_list
pvalue.final<-cbind(TF_mean,non_TF_mean,pvalue)
colnames(pvalue.final)<-c("Tissue","Cancer","pvalue")
write.table(pvalue.final,"t.test.pvalue.normalTssue_map2_cancer.mean.txt",quote = F,sep='\t')

######???????𵥶��???
cancer_list<-list2
rownum<- length(cancer_list)
pvalue<-0
TF_mean<-0
non_TF_mean<-0

  HNSC<-data[which(data$Group==cancer_list[4]),]
  unique(HNSC$Tissue)
  
  HNSC_non_TF<-HNSC[which(HNSC$Tissue=="LUAD"),]
  HNSC_TF<-HNSC[which(HNSC$Tissue=="Lung"),]
  tmp<-t.test(HNSC_non_TF$TPM,HNSC_TF$TPM, paired=FALSE)
  pvalue<-tmp$p.value
  HNSC_non_TF_mean<-median(HNSC_non_TF$TPM)
  HNSC_TF_mean<-median(HNSC_TF$TPM)
  TF_mean<-HNSC_TF_mean
  non_TF_mean<-HNSC_non_TF_mean

pvalue<-as.matrix(pvalue)
TF_mean<-as.matrix(TF_mean)
non_TF_mean<-as.matrix(non_TF_mean)
rownames(pvalue)<-unique(HNSC_non_TF$Tissue)
rownames(TF_mean)<-unique(HNSC_non_TF$Tissue)
rownames(non_TF_mean)<-unique(HNSC_non_TF$Tissue)
pvalue.final<-cbind(TF_mean,non_TF_mean,pvalue)
pvalue.final
colnames(pvalue.final)<-c("Tissue","Cancer","pvalue")
pvalue.final


#################################
setwd("C:\\项目\\TCGA_PAAD")
library("reshape")
library("ggplot2")
#tissue_list<-unique(data$Sample)

co<-c("DeepSkyBlue","Tomato")
co<-c("blue","LightSkyBlue")
#my_comparisons <- list( c("TFs", "non-TFs"))
#p<-ggplot(data=data, aes(x=Sample,y=Value))+geom_boxplot(aes(fill=flag))+theme_bw()
normal<-read.delim("mrna_expr_tpm.normal.txt",header=T,sep='\t',row.names = 1)
tumor<-read.delim("mrna_expr_tpm.tumor.txt",header=T,sep='\t',row.names = 1)


Normal_group <- normal[,substring(colnames(normal),1,12) %in% substring(colnames(tumor),1,12)]
tumor_group<-tumor[,substring(colnames(tumor),1,12) %in% substring(colnames(normal),1,12)]

tumor_group.1 <- tumor_group[,match(substring(colnames(Normal_group),1,12),substring(colnames(tumor_group),1,12))]
colnames(Normal_group)
colnames(tumor_group.1)

expr_data<-cbind(Normal_group,tumor_group.1)

colnames(expr_data)
exp2<-expr_data
library(stringr)

rownum<- dim(exp2)[1]
colnum<-dim(exp2)[2]
pvalue<-0
T_mean<-0
N_mean<-0

for(n in 1:rownum){
  
  tmp<-t.test(as.numeric(exp2[n,1:(colnum/2)]),as.numeric(exp2[n,(colnum/2+1):colnum]), paired = TRUE)
  pvalue[n]<-tmp$p.value
  Tumor_mean<-mean(as.numeric(exp2[n,1:(colnum/2)]))
  Normal_mean<-mean(as.numeric(exp2[n,(colnum/2+1):colnum]))
  
  T_mean[n]<-Tumor_mean
  N_mean[n]<-Normal_mean
}
pvalue<-as.matrix(pvalue)
T_mean<-as.matrix(T_mean)
N_mean<-as.matrix(N_mean)
rownames(pvalue)<-rownames(exp2)
rownames(T_mean)<-rownames(exp2)
rownames(N_mean)<-rownames(exp2)
pvalue.final<-cbind(T_mean,N_mean,pvalue)
head(pvalue.final)
colnames(pvalue.final)<-c("Tumor_mean","Adjacent_mean","pvalue")
write.table(pvalue.final,"t.test.pvalue.PAAD_tumor2Adjacent.txt",quote = F,sep='\t')
