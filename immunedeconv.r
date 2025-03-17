
library(devtools)
Sys.setlocale(category ="LC_ALL",locale = "us")  #调整内码格式
install_github("icbi-lab/immunedeconv")



setwd("D:\\项目\\TCGA_PAAD\\immuncellAI")
df<-read.delim("D:\\项目\\TCGA_PAAD\\mrna_expr_tpm.tumor.txt",header=T,row.names = 1)
dim(df)
df[1:4,1:5]

#教程https://mp.weixin.qq.com/s/RRUKsHEbg9YiMmZoqTfoIQ
rm(list=ls())
#load("PAAD_tpm.Rdata")
library(immunedeconv)
method=c('quantiseq','xcell','estimate','epic','mcp_counter','timer','cibersort','cibersort_abs',)
method=c('quantiseq','xcell','estimate','epic',
         'mcp_counter','timer','cibersort','cibersort_abs')
res1=deconvolute(df,'quantiseq')

res2=deconvolute(df,'xcell')
write.table(res2, file = "PAAD.tumor.sample.immune.xcell.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
####由于临床信息行名中样本为12位数，需要将免疫细胞矩阵调整
cell_ratio<-res2
newcolname<-gsub("\\.","-",colnames(cell_ratio))
colnames(cell_ratio)<-newcolname
rownames(cell_ratio)<-cell_ratio$cell_type

#
load("D:/项目/TCGA_PAAD/GEO/dat.1.Rdata")
head(dat.1)
dim(dat.1)
#
immucell=t(cell_ratio)%>%
  as.data.frame() %>%
  tibble::add_column(ID=stringr::str_sub(rownames(.),1,12))%>%
  dplyr::filter(!duplicated(ID)) %>% # remove duplicated samples randomly
  tibble::remove_rownames(.) %>%
  tibble::column_to_rownames("ID")%>%
  dplyr::filter(rownames(.) %in% rownames(dat.1))

identical(rownames(immucell),rownames(dat.1))
#[1] TRUE
immucell<-immucell[rownames(dat.1),]
identical(rownames(immucell),rownames(dat.1))
write.table(immucell, file = "PAAD.tumor.sample.immune.xcell.final.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)

######################################
setwd("C:\\项目\\TCGA_PAAD\\immuncellAI")
immu=read.table("ImmuCellAI_abundance_result.onlyTCGAtumorsample.txt",header=T, row.names=1)
immucell=read.delim("PAAD.tumor.sample.immune.xcell.final.txt",header=T, row.names=1)
head(immucell)
list<-substring(rownames(immu),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
rownames(immu)<-list1[,1]
head(immu)

load("C:/项目/TCGA_PAAD/GEO/dat.1.Rdata")
head(dat.1)
dim(dat.1)

immucell<-immucell[rownames(dat.1),]
head(immucell);dim(immucell)

#
group<-dat.1$Riskgroup
dat<-as.data.frame(cbind(immucell[,1:36],group))
#
####相关性分析

library(psych)
#immucell<-apply(immucell,1,function(x){as.numeric(x)})
immucell.1<-immucell[,1:36]
t1<-as.matrix(immucell.1)
#test<-apply(immucell,1,function(x){as.numeric(x)})

test <- matrix(
  as.numeric(t1),ncol=36)
rownames(test)<-rownames(t1)
colnames(test)<-colnames(t1)

psych::corr.test(dat.1$riskscore, test, method = 'pearson',adjust="none")
#cor$p
correlation<-c()
p<-c()

cor <-psych::corr.test(dat.1$riskscore, test, method = 'pearson',adjust="none")
cmt <-t(as.data.frame(cor$r))%>%
  as.data.frame()%>%
  tibble::add_column(rownames(.))%>%
  tibble::remove_rownames(.)

cmt2 <-t(as.data.frame(cor$p))%>%
  as.data.frame()%>%
  tibble::add_column(rownames(.))%>%
  tibble::remove_rownames(.)

colnames(cmt)=c("Correlation_Coefficient","Immune_cell")
colnames(cmt2)=c("Pvalue","Immune_cell")
immune.final<-merge(cmt,cmt2,by="Immune_cell")
write.table(immune.final, file = "immune.xcell.and.riskscore.correlation.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
#绘图气泡图（没有基因）

#画图
library(ggplot2)
y_cols <- rep(c("#f27767","#bc9826","#53b449","#23b892","#1cb6ec","#d269a4"),each=5)

ggplot()+
  geom_point(data=df,aes(`Correlation Coefficient`,`Immune cell`,color=Software),
             shape=16,size=5.7)+
  scale_color_manual(values = c('QUANTISEQ' = "#f27767",'XCELL' = "#bc9826",
                                'EPIC'="#53b449",'MCPCOUNTER'="#23b892",
                                'TIMER'="#1cb6ec",'CIBERSORT'="#d269a4"))+
  theme(axis.text.y=element_text(color=y_cols,size=10),
        axis.text.x=element_text(size=10,color = 'black'),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

#
head(immune.final)
data<-immune.final
data<-immune.final %>% 
  subset(Immune_cell!="immune score") %>% 
  subset(Immune_cell!="stroma score") %>% 
  subset(Immune_cell!="microenvironment score")

library(ggplot2)
head(data)
data$LogP<--log10(data$Pvalue)
pp = ggplot(data,aes(Correlation_Coefficient,Immune_cell,size=Correlation_Coefficient))
p.final<-pp + geom_point(aes(color=(LogP)))+theme_bw()+
  labs(x="Correlation_Coefficient",y="",title="",
       color=expression(-log[10](Pvalue)))+
  theme(axis.text.x = element_text(size = 12, 
                                   angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=11,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        plot.title = element_text(size=13,colour="black"),
        legend.key.size = unit(5,'mm'))+
  scale_y_discrete(limit = unique(data$Description))+
  scale_colour_gradient(low="PaleTurquoise",high="Turquoise4",
                        name="-log P")
#low="LightSkyBlue",high="MediumBlue"
p.final
ggsave(filename = "immune.xcell.and.riskscore.correlation.4gene.pdf", 
       plot =p.final,width = 15, height =15, units = 'cm')
# 绘图热图（基因和免疫成分相关性）

#添加一列,来判断pvalue值范围
final<-immune.final
final$sign<-case_when(final$Pvalue<0.05 &final$Pvalue>0.01 ~"*",
                      final$Pvalue<0.01 &final$Pvalue>0.001 ~"**",
                      final$Pvalue<0.001 ~"***",
                      final$Pvalue>0.05 ~"")
ggplot(data=final,aes(x=Gene,y=celltype))+
  geom_tile(aes(fill=correlation),colour="white",size=1)+
  scale_fill_gradient2(low="#2b8cbe",mid="white",high="#e41a1c")+
  geom_text(aes(label=sign),colour="black")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  guides(fill=guide_legend(title="* p<0.05\n\n** p<0.01\n\n*** p<0.001\n\ncorrelation"))
ggsave("correlation.pdf",width=12,height=13)

#####分组
#group <-  as.vector(ifelse(data.final$riskscore > median(data.final$riskscore), "high", "low"))#定义低、高风险
#dat <- cbind(data.final, group) 
#data.final[1:3,1:4];dim(data.final)
#dim(test)
#test[test<0.001]<-0
#dat<-as.data.frame(cbind(test[,1:36],group))

group<-dat.1$Riskgroup
dat<-as.data.frame(cbind(immucell[,1:36],group))
dat<-as.data.frame(cbind(immucell[,c("Myeloid.dendritic.cell.activated","T.cell.CD8.",
                                     "T.cell.CD8..central.memory",
                                     "T.cell.CD8..effector.memory",
                                     "Class.switched.memory.B.cell",
                                     "Common.lymphoid.progenitor",
                                     "Common.myeloid.progenitor",
                                     "Endothelial.cell",
                                     "Cancer.associated.fibroblast",
                                     "Hematopoietic.stem.cell",
                                     "Macrophage.M2",
                                     "B.cell.memory",
                                     "T.cell.NK",
                                     "T.cell.CD4..Th2")],group))

dat<-as.data.frame(cbind(immucell[,c(1,4,6,9,10,12,13,14,16,18,20,23,25,29,30,35)],group))
#dat[dat<0.001]<-0
library("reshape")
dat_melt<-melt(dat,id=c("group"))
colnames(dat_melt)<-c("group","celltype","Propotion")
##boxplot
dat_melt$celltype<-factor(dat_melt$celltype,
                            levels=c("B.cell.memory",
                                     "Cancer.associated.fibroblast",
                                     "Class.switched.memory.B.cell",
                                     "Common.lymphoid.progenitor",
                                     "Common.myeloid.progenitor",
                                     "Endothelial.cell",
                                     "Hematopoietic.stem.cell",
                                     "Macrophage.M2",
                                     "Myeloid.dendritic.cell.activated",
                                     "NK.cell",
                                     "T.cell.CD4..central.memory",
                                     "T.cell.CD4..naive",
                                     "T.cell.CD4..Th2",
                                     "T.cell.CD8.",
                                     "T.cell.CD8..central.memory",
                                     "T.cell.NK"))

dat_melt$Propotion<-as.numeric(dat_melt$Propotion)
dat_melt$group<-factor(dat_melt$group,levels=c("high","low"))
p<-ggplot(data=dat_melt, aes(x=group,y=Propotion))+
  geom_boxplot(aes(fill=group),outlier.colour = NA,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p
##
library("ggpubr")
library(RColorBrewer)
display.brewer.all()
mycol <- colorRampPalette(brewer.pal(10,'Paired'))(10)
mycol <- colorRampPalette(brewer.pal(3,'Set2'))(3)

#p <- ggboxplot(dat_melt, x = "celltype", y = "Propotion",
 #              color = "group", palette = "jco",
  #             add = "jitter")
#p + stat_compare_means(aes(group = group))

# palette =c(mycol),
#dat_melt.1<-dat_melt %>% subset(celltype %in% c("T cell CD4+ memory","T cell CD4+ naive"))
dat_melt$Propotion<-as.numeric(dat_melt$Propotion)
dat_melt$group<-factor(dat_melt$group,levels=c("high","low"))
p<-ggboxplot(dat_melt, x = "celltype", y = "Propotion",color = "group",
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
        axis.title.x = element_text(size=12,colour="black",face="bold"))+ylim(0,0.5)+
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

p2<-p +stat_compare_means(aes(group = group),label = "p.signif",
                     size=4,method="wilcox.test")+
 theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_risk.score.xcell.4gene.pdf",plot=p2,
       device='pdf',path=".",width=15,height=4.5)

ggsave(filename="boxplot_risk.score.xcell.partcelltype.4gene.pdf",plot=p2,
       device='pdf',path=".",width=7.3,height=4.5)
###
p<-ggplot(dat_melt, aes(x=celltype, y=Propotion,fill=group)) +
  geom_violin(trim=FALSE,color="white",position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线???
  scale_fill_manual(values = c("IndianRed1","MediumTurquoise"))+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1,colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+ylim(0,1.6)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="",x="",y="Propotion",title="")
p
my_comparisons <- list(c("high", "low"))

#stat_compare_means(comparisons = my_comparisons
#orig.ident 
#label = "p.signif" "p.format"
p2<-p +
  stat_compare_means(aes(group = group),label = "p.signif",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="violinplot_risk.score.xcell.celltype.4gene.pdf",plot=p2,
       device='pdf',path=".",width=15,height=6)
##################################################
#######调用estimate计算免疫细胞比例

res3=deconvolute_estimate(df)
save(res3,file="res3.RData")

load("res3.RData")
cell_ratio<-res3
newcolname<-gsub("\\.","-",colnames(cell_ratio))
colnames(cell_ratio)<-newcolname
#rownames(cell_ratio)<-cell_ratio$cell_type

estimate.result=t(cell_ratio)%>%
  as.data.frame() %>%
  tibble::add_column(ID=stringr::str_sub(rownames(.),1,12))%>%
  dplyr::filter(!duplicated(ID)) %>% # remove duplicated samples randomly
  tibble::remove_rownames(.) %>%
  tibble::column_to_rownames("ID")%>%
  dplyr::filter(rownames(.) %in% rownames(dat.1))

identical(rownames(estimate.result),rownames(dat.1))
#[1] TRUE
estimate.result<-estimate.result[rownames(dat.1),]
identical(rownames(estimate.result),rownames(dat.1))
write.table(estimate.result, file = "PAAD.tumor.sample.estimate.result.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)

#estimate.result<-read.delim("PAAD.tumor.sample.estimate.result.txt",header=T,row.names = 1)

group<-dat.1$Riskgroup
#dat<-as.data.frame(cbind(immucell[,1:36],group))
#group <-  as.vector(ifelse(data.final$riskscore > median(data.final$riskscore), "high", "low"))#定义低、高风险
dat2 <- cbind(estimate.result, group) 
dat2[1:3,1:5];dim(dat2)

library("reshape")
dat_melt2<-melt(dat2,id=c("group"))
colnames(dat_melt2)<-c("group","type","value")
dat_melt2$group<-factor(dat_melt2$group,levels=c("high","low"))
p<-ggplot(data=dat_melt2, aes(x=group,y=value))+geom_boxplot(aes(fill=group),outlier.colour = NA,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p

dat_melt2$value<-as.numeric(dat_melt2$value)
p<-ggboxplot(dat_melt2, x = "group", y = "value",color = "group",
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
  labs(fill="type",x="",y="score",title="")

my_comparisons <- list( c("high", "low"))

p2<-p+facet_wrap(~ type, scales="free",nrow = 1) + 
stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                   size=4,method="wilcox.test")+
theme(strip.text.x = element_text(size=11, color="black",
                                  face="bold"))+
theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_risk.score.estimate.4gene.pdf",plot=p2,
       device='pdf',path=".",width=6.5,height=4)

###相关性点图
estimate.result<-read.delim("PAAD.tumor.sample.estimate.result.txt",header=T,row.names = 1)

#group<-dat.1$Riskgroup
#dat<-as.data.frame(cbind(immucell[,1:36],group))
#group <-  as.vector(ifelse(data.final$riskscore > median(data.final$riskscore), "high", "low"))#定义低、高风险
dat2 <- cbind(estimate.result, dat.1) 
dat2[1:3,1:5];dim(dat2);head(dat2)

library(ggplot2)
library(ggpubr)
library(ggExtra)
x<-as.numeric(dat2$riskscore)
y<-as.numeric(dat2$StromalScore)
y<-as.numeric(dat2$ImmuneScore)
y<-as.numeric(dat2$ESTIMATEScore)
y<-as.numeric(dat2$TumorPurity)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
  xlab("riskscore")+ylab("TumorPurity")+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", 
              xparams = list(fill = "Burlywood1"),
              yparams = list(fill = "SkyBlue1"))
p2
pdf(file="cor.density.TumorPurity.pdf",width=4,height=4)
print(p2)
dev.off()
#######################################################################
###############################################
setwd("C:\\项目\\TCGA_PAAD\\immuncellAI")
immu=read.table("ImmuCellAI_abundance_result.onlyTCGAtumorsample.txt",header=T, row.names=1)
list<-substring(rownames(immu),1,12)
list1<-as.data.frame(gsub("\\.", "-", list))
rownames(immu)<-list1[,1]
head(immu)

load("C:/项目/TCGA_PAAD/GEO/dat.1.Rdata")
head(dat.1)
dim(dat.1)

immu<-immu[rownames(dat.1),]
head(immu);dim(immu)

#
group<-dat.1$Riskgroup
dat<-as.data.frame(cbind(immu[,1:24],group))
#dat[dat<0.001]<-0
library("reshape")
dat_melt<-melt(dat,id=c("group"))
colnames(dat_melt)<-c("group","celltype","Propotion")
##boxplot
p<-ggplot(data=dat_melt, aes(x=group,y=Propotion))+
  geom_boxplot(aes(fill=group),outlier.colour = NA,show.legend = FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+
  labs(fill="",y="",title="")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p
##
library("ggpubr")
library(RColorBrewer)
display.brewer.all()
mycol <- colorRampPalette(brewer.pal(10,'Paired'))(10)
mycol <- colorRampPalette(brewer.pal(3,'Set2'))(3)


dat_melt.1<-dat_melt %>% 
  subset(celltype %in% c("Bcell","CD4_T",	"CD8_T","DC",
                         "Gamma_delta", 
                         "Monocyte",	"Macrophage",	
                         "NK",	"Neutrophil",	"NKT"
  ))
dat_melt.1$celltype<-factor(dat_melt.1$celltype,
                            levels=c("Bcell","CD4_T",	"CD8_T","DC",
                                     "Gamma_delta", 
                                     "Monocyte",	"Macrophage",	
                                     "NK",	"Neutrophil",	"NKT"))

dat_melt.1<-dat_melt %>% 
  subset(celltype %in% c("CD4_naive","CD8_naive",	"Central_memory",
                         "Cytotoxic","Effector_memory","Exhausted",
                         "iTreg","MAIT",	"nTreg",
                         "Tr1",	"Th1",	"Th2",	"Th17",
                         "Tfh"))
dat_melt.1$celltype<-factor(dat_melt.1$celltype,
                            levels=c("CD4_naive","CD8_naive",	"Central_memory",
                                     "Cytotoxic","Effector_memory","Exhausted",
                                     "iTreg","MAIT",	"nTreg",
                                     "Tr1",	"Th1",	"Th2",	"Th17","Tfh"))
head(dat_melt.1)
dat_melt.1$Propotion<-as.numeric(dat_melt.1$Propotion)
dat_melt.1$group<-factor(dat_melt.1$group,levels=c("high","low"))
p<-ggboxplot(dat_melt.1, x = "celltype", y = "Propotion",color = "group",
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

p2<-p +stat_compare_means(aes(group = group),label = "p.signif",
                          size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="boxplot_risk.score.immucellAI.celltype.4gene.1.pdf",plot=p2,
       device='pdf',path=".",width=5.5,height=5)
###
p<-ggplot(dat_melt, aes(x=celltype, y=Propotion,fill=group)) +
  geom_violin(trim=FALSE,color="white",position = "dodge",scale="width") +
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线???
  scale_fill_manual(values = c("IndianRed1","MediumTurquoise"))+
  theme_classic()+
  theme(legend.text=element_text(size=10,colour="black"),
        plot.title = element_text(size=15,colour="black"),
        legend.key.size = unit(4,'mm'))+
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1,colour = "black",face="bold"),
        axis.text.y = element_text(size=10,colour="black",face="bold"),
        legend.title=element_text(size=10))+
  theme(axis.title.y = element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12,colour="black",face="bold"))+ylim(0,1.6)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="",x="",y="Propotion",title="")
p
my_comparisons <- list(c("high", "low"))

#stat_compare_means(comparisons = my_comparisons
#orig.ident 
#label = "p.signif" "p.format"
p2<-p +
  stat_compare_means(aes(group = group),label = "p.signif",
                     size=4,method="wilcox.test")+
  theme(strip.text.x = element_text(size=11, color="black",
                                    face="bold"))+
  theme(strip.background = element_blank(),strip.placement = "outside")
p2

ggsave(filename="violinplot_risk.score.xcell.celltype.pdf",plot=p2,
       device='pdf',path=".",width=15,height=6)
##################################################
####相关性分析
library(psych)
#immucell<-apply(immucell,1,function(x){as.numeric(x)})
immu.1<-immu[,1:24]
t1<-as.matrix(immu.1)
#test<-apply(immucell,1,function(x){as.numeric(x)})

test <- matrix(
  as.numeric(t1),ncol=24)
rownames(test)<-rownames(t1)
colnames(test)<-colnames(t1)

identical(rownames(test),rownames(dat.1))
psych::corr.test(dat.1$riskscore, test, method = 'pearson',adjust="none")
#cor$p
correlation<-c()
p<-c()

cor <-psych::corr.test(dat.1$riskscore, test, method = 'pearson',adjust="none")
cmt <-t(as.data.frame(cor$r))%>%
  as.data.frame()%>%
  tibble::add_column(rownames(.))%>%
  tibble::remove_rownames(.)

cmt2 <-t(as.data.frame(cor$p))%>%
  as.data.frame()%>%
  tibble::add_column(rownames(.))%>%
  tibble::remove_rownames(.)

colnames(cmt)=c("Correlation_Coefficient","Immune_cell")
colnames(cmt2)=c("Pvalue","Immune_cell")
immune.final<-merge(cmt,cmt2,by="Immune_cell")
write.table(immune.final, file = "immune.immucellAI.and.riskscore.correlation.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
#绘图气泡图（没有基因）

#
head(immune.final)
data<-immune.final
#data<-immune.final %>% 
 # subset(Immune_cell!="immune score") %>% 
  #subset(Immune_cell!="stroma score") %>% 
  #subset(Immune_cell!="microenvironment score")

library(ggplot2)
head(data)
data$LogP<--log10(data$Pvalue)
pp = ggplot(data,aes(Correlation_Coefficient,Immune_cell,size=Correlation_Coefficient))
p.final<-pp + geom_point(aes(color=(LogP)))+theme_bw()+
  labs(x="Correlation_Coefficient",y="",title="",
       color=expression(-log[10](Pvalue)))+
  theme(axis.text.x = element_text(size = 12, 
                                   angle = 45, hjust = 1,
                                   vjust =1,colour = "black",lineheight=100),
        axis.text.y = element_text(size=11,colour="black"),
        legend.title=element_text(size=12,colour="black"),
        legend.text=element_text(size=11,colour="black"),
        axis.title.x = element_text(size=12,colour="black"),
        axis.title.y = element_text(size=12,colour="black"),
        plot.title = element_text(size=13,colour="black"),
        legend.key.size = unit(5,'mm'))+
  scale_y_discrete(limit = unique(data$Description))+
  scale_colour_gradient(low="PaleTurquoise",high="Turquoise4",
                        name="-log P")
#low="LightSkyBlue",high="MediumBlue"
p.final
ggsave(filename = "immune.immucellAI.and.riskscore.correlation.new.pdf", 
       plot =p.final,width = 13, height =13, units = 'cm')

