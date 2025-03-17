setwd("D:\\项目\\TCGA_PAAD")
exp<-read.delim("D:\\项目\\TCGA_PAAD\\TCGA_Gtex\\Download_from_Gtex_database\\PAAD_pancreas_aferRMBatch_effect.count.txt",header=T,sep='\t')
input<-read.delim("D:\\项目\\TCGA_PAAD\\GEO\\GSE149103_RNA_seq\\hg19.symbol.length.final.txt",header=T,sep='\t')

library(dplyr)
merge<-left_join(exp,input,by="Symbol")#根据基因那列进行合并
mycounts <- na.omit(merge)#删除错误值行
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
mycounts[1:4,1:5]
#TPM计算
kb <- mycounts$length / 1000
kb
countdata <- mycounts[,-dim(mycounts)[2]]
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.table(tpm,file="PAAD_pancreas_aferRMBatch_effect.tpm.txt",quote = F,sep='\t',row.names=T)

########################################################3
setwd("D:\\项目\\TCGA_PAAD\\WGCNA")
library("WGCNA")

options(stringsAsFactors = FALSE)
list.files()
#femData = read.table("mrna_expr_tpm.tumor.txt",header=TRUE) 

tpm<-read.delim("D:\\项目\\TCGA_PAAD\\PAAD_pancreas_aferRMBatch_effect.tpm.txt",header=T,row.names = 1)
group<-c(rep("Cancer",179),rep("Normal",332))
group<-as.data.frame(as.matrix(group))
rownames(group)<-colnames(tpm)

tpm<-read.delim("D:\\项目\\TCGA_PAAD\\PAAD_pancreas_aferRMBatch_effect.tpm.txt",header=T,row.names = 1)
femData<-tpm
femData[1:3,1:4]
datExpr = as.data.frame(t(femData))
datExpr[1:3,1:4]
#colnames(datExpr0) = rownames(femData) ###此处GENE_ID需要注意
#rownames(datExpr0) = names(femData)
#####对数据进行过滤
nSamples = nrow(datExpr)  # 统计样品数目
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))  #按列（基因）取方差
no.missingdatExpr=as.vector(apply(is.na(as.matrix(datExpr)),2, sum) )#按列（基因）统计缺失数目
KeepGenes= variancedatExpr>0 & no.missingdatExpr<0.1*nSamples  # 保留方差不等于0，且缺失低于10%的基因
table(KeepGenes)# 过滤统计

datExpr=datExpr[, KeepGenes]

#datExpr0[datExpr0==0]<-0.01
#datExpr0_1<-apply(datExpr0,2,function(x){as.numeric(x)})

#datExpr0_2<-datExpr0_1
###datExpr0_2<-log10(datExpr0_1)


datExpr0<-datExpr
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###样本聚类检测离群样本
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 250000, minSize = 0)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)

###网络构建及模块检测
key3245gene<-read.delim("D:\\项目\\TCGA_PAAD\\TCGA_Gtex\\Download_from_Gtex_database\\3245DEG.txt",header=T,row.names = NULL)
intersect(key3245gene$Symbol,colnames(datExpr))
datExpr<-datExpr[,c(intersect(key3245gene$Symbol,colnames(datExpr)))]
dim(datExpr)
# Choose a set of soft-thresholding powers筛选阈值
powers = c(c(1:20), seq(from = 22, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "./soft_cut.pdf", width = 12, height = 9);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft[[1]]
##Constructing the gene network and identifying modules is now a simple function call
net = blockwiseModules(datExpr, power = sft[[1]],
                       TOMType = "unsigned", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
table(net$colors)



# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "./cluster_dendrogram.pdf", width = 12, height = 9);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "2-networkConstruction-auto.RData")

#注意：输出模块分类结果：
module_class<-as.matrix(net$colors)
#module_class<-as.matrix(moduleLabels)
module_color<-as.matrix(moduleColors)
rownames(module_class)<-colnames(datExpr)
rownames(module_color)<-colnames(datExpr)
write.table(module_class,file="./all_module_class.xls",sep='\t',quote = F)
write.table(module_color,file="./all_module_color.xls",sep='\t',quote = F)


####moule 热图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft[[1]]);
TOM=1-dissTOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^(sft[[1]]+1);
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
#sizeGrWindow(9,9)
pdf(file = "./heatmap.pdf", width = 12, height = 9);
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()
save(TOM,file = "TOM.RData")

dissTOM = 1-TOM 
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, 
     xlab="", 
     sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, 
     hang = 0.04);
dev.off()
#######################
###计算Module Eigengenes（相关概念在上一专题进行具体阐述）并合并相似模块
#1、计算Module Eigengenes

MEList = moduleEigengenes (datExpr, colors = moduleColors) #按照模块计算每个module的ME（也就是该模块的第一主成分）
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs); # 计算根据模块特征向量基因计算模块相异度
METree = hclust(as.dist(MEDiss), method = "average"); #对不同的模块进行聚类
pdf("模块相关系数.pdf")
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                                            marHeatmap = c(3,4,2,2),
                                            plotDendrograms = FALSE,
                                            xLabelsAngle = 90)

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=0.99, col="red")
dev.off()
#2、对相似度高的模块进行合并并绘制热图。此次没有进行任何合并
merge_modules = mergeCloseModules(datExpr, moduleColors, cutHeight = 0.1, verbose = 3)
mergedColors = merge_modules$colors; # 合并后的颜色：
mergedMEs = merge_modules$newMEs; # 新模块的特征向量MEs：

pdf("module_cluster_merge.pdf")
plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors),
                    
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    
                    dendroLabels = FALSE, hang = 0.03,
                    
                    addGuide = TRUE, guideHang = 0.05)

dev.off()
######表型和模块的相关性
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs = net$MEs;
group<-as.data.frame(group[rownames(MEs),])
rownames(group)<-rownames(MEs)
colnames(group)<-c("group")
table(group$group)

datTraits = data.frame(samples=rownames(group),subtype=group$group)
design=model.matrix(~0 + datTraits$subtype)
colnames(design)=levels(factor(datTraits$subtype))

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes #计算合并后的模块的第一主成分ME。

MEs = orderMEs(MEs0); #不同颜色的模块的ME值矩阵。

moduleTraitCor = cor(MEs, design , use = "p") #计算不同模块和样本性状的相关性。

nSamples = nrow(datExpr)

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) #计算不同模块和样本性状相关性的p值。

#可视化module与性状的相关性和它们的p值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf('group_Modules_heatmap.pdf',width = 4.5, height = 6)
labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = colnames(design),
                              xLabelsAngle = 0,
                              yLabels = row.names(moduleTraitCor),
                              ySymbols = row.names(moduleTraitCor),
                              colorLabels = TRUE,
                              colors = blueWhiteRed(50),
                              textMatrix = textMatrix,
                              setStdMargins = FALSE,
                              cex.text = 1,
                              zlim = c(-1,1),
                              main = paste("Module-trait relationships"))
dev.off()
#####module网络
#########3Exporting to Cytoscape
setwd("D:\\项目\\TCGA_PAAD\\WGCNA")
# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr, power =sft[[1]]);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
uniq_module_color<-unique(moduleColors)
num<-length(uniq_module_color)
#femData$Symbol<-rownames(femData)
for(i in 1:num){
  modules =uniq_module_color[i] ;  
  ###或者 modules = "red";
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modGenes = colnames(datExpr)[match(modProbes, colnames(datExpr))];
  #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)]
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.001,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])};
#save.image("final_result.Rdata")
#CytoscapeInput-edges-blue
#CytoscapeInput-edges-turquoise
#gene.107

edgefiles<-list.files(pattern = "CytoscapeInput-edges*")

alledge=data.frame()
for (i in 1:length(edgefiles)) {
  
  list1<-read.delim(edgefiles[i],header=T,sep='\t')
  edgefilter<-list1[which(list1$weight >0.25),]
  alledge=rbind(alledge,edgefilter)
}
head(alledge);dim(alledge)

#list1<-read.delim("CytoscapeInput-edges-blue.txt",header=T,sep='\t')
#list1.1<-list1[which(list1$weight >0.2),]
alledge.1<-subset(alledge,fromNode %in% c(key132gene$key132gene))
a<-unique(alledge.1$toNode)

d<-union(a,key132gene$key132gene)
length(d)
write.table(d,file="immuneAndEMTgene_coexp_gene.358gene.txt",quote = F,sep="\t",row.names=T)

table(alledge.1$toNode)
test<-as.data.frame(as.matrix(table(alledge.1$toNode)))
test$symbol<-rownames(test)
b<-rownames(test[which(test$V1>=2),])
d<-union(b,key132gene$key132gene)
write.table(d,file="immuneAndEMTgene_coexp_gene.233gene.txt",quote = F,sep="\t",row.names=T)

gene233<-read.delim("immuneAndEMTgene_coexp_gene.233gene.txt",header=T,sep='\t',row.names = 1)
head(alledge);dim(alledge)

alledge.3<-subset(alledge,fromNode %in% c(gene233$x))
alledge.4<-subset(alledge.3,toNode %in% c(gene233$x))

write.table(alledge.4,file="immuneAndEMTgene_coexp_gene.233gene.edge.txt",quote = F,sep="\t",row.names=T)
write.table(alledge.4,file="immuneAndEMTgene_coexp_gene.233gene.edge.xls",quote = F,sep="\t",row.names=T)
