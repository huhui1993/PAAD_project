setwd("C:\\项目\\TCGA_PAAD")
normal<-read.delim("mrna_expr_counts.normal.txt",header=T,sep='\t',row.names = 1)
tumor<-read.delim("mrna_expr_counts.tumor.txt",header=T,sep='\t',row.names = 1)

###提取癌和癌旁成对的配对样本
Normal_group <- normal[,substring(colnames(normal),1,12) %in% substring(colnames(tumor),1,12)]
tumor_group<-tumor[,substring(colnames(tumor),1,12) %in% substring(colnames(normal),1,12)]

tumor_group.1 <- tumor_group[,match(substring(colnames(Normal_group),1,12),substring(colnames(tumor_group),1,12))]
colnames(Normal_group)
colnames(tumor_group.1)

test_data<-cbind(Normal_group,tumor_group.1)
write.table(test_data,file="PAAD.noraml.and.tumor.counts.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )
test_data1<-log2(test_data+1)
write.table(test_data1,file="PAAD.noraml.and.tumor.counts.log2.xls",
    quote = F,sep='\t',row.names=T,col.names =NA )

test2<-read.delim("C:\\项目\\TCGA_PAAD\\immuncellAI\\PAAD.noraml.and.tumor.counts.log2.addgroup.new.txt",
                  header=T,sep='\t',row.names = 1)

library(DESeq2)
condition <- factor(c(rep("con",dim(Normal_group)[2]),rep("case",dim(tumor_group.1)[2])), levels = c("con","case"))
colData <- data.frame(row.names=colnames(test_data), condition)
dds <- DESeqDataSetFromMatrix(test_data, colData, design= ~ condition)
dds <- DESeq(dds)
# 查看一下dds的内容
dds
#接下来，我们要查看case versus control的总体结果，并根据padj进行重新排序。利用summary命令统计显示一共多少个genes上调和下调（FDR0.1）
res = results(dds, contrast=c("condition", "con", "case"))
#或下面命令
res= results(dds)
res = res[order(res$padj),]
head(res)
summary(res)
#所有结果先进行输出,未筛选
write.csv(res,file="All_results.csv")

####将表达量矩阵按照差异结果进行排序
exp_data<-test_data[rownames(res),]

####计算正常和疾病样本表达均值(这里是成对的正常样本和癌旁)
colnum<-ncol(exp_data)
test_con<-exp_data[1:(colnum/2)]
test_case<-exp_data[(colnum/2+1):colnum]
test_con$mean_con<-apply(test_con,1,mean)
test_case$mean_case<-apply(test_case,1,mean)
####计算正常和疾病样本表达均值(这里是不成对时的正常样本和癌旁),前面表达量矩阵test_data<-cbind(Normal_group,tumor_group.1)
test_con<-exp_data[,colnames(Normal_group)]
test_case<-exp_data[,colnames(tumor_group.1)]
test_con$mean_con<-apply(test_con,1,mean)
test_case$mean_case<-apply(test_case,1,mean)

###合并计算出来的均值和差异结果
res_addMean<-data.frame(test_con$mean_con,test_case$mean_case,res)
#####输出未筛选的差异结果
write.table(res_addMean,file="nofilter.xls",quote = F,sep='\t',row.names=T,col.names =NA )
#####提取筛选后的差异表达基因
##1.去掉低表达
deseq_gene<-res_addMean[apply(res_addMean[,1:2],1,function(x){XXXX<-FALSE;if(any(x>5)){XXXX<-TRUE};return(XXXX)}),]
dim(deseq_gene)
##2.根据padj和FC进一步筛选
diff_gene_deseq2 <-subset(deseq_gene, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
##输出筛选后的差异结果
write.table(diff_gene_deseq2,file="diff_gene_deseq2.FDR0.5.FC2.Higher5.xls",
            quote = F,sep='\t',row.names=T,col.names =NA )
