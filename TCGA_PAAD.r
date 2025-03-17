setwd("C:\\??Ŀ\\TCGA_PAAD")
setwd("C:\\项目\\TCGA_PAAD")
setwd("C:\\项目\\TCGA_PAAD")
BiocManager::install("TCGAbiolinks")
library("TCGAbiolinks")
library(dplyr)
library(DT)
TCGAbiolinks:::getGDCprojects()$project_id
TCGAbiolinks:::getProjectSummary("TCGA-PAAD")

#clinical????
query <- GDCquery(
  project = "TCGA-PAAD", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

patient_info <- clinical.BCRtab.all$clinical_patient_paad
dim(patient_info)
#[1] 1099  112
patient_info[1:6, 1:6]
######################333
clinical_data <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")
write.table(clinical_data,file="TCGA-PAAD_185sample_clinical.txt",quote = F,sep="\t",row.names=T)


query <- GDCquery(
  project = "TCGA-PAAD", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab",
  file.type = "radiation"
)
GDCdownload(query)
clinical.BCRtab.radiation <- GDCprepare(query)
clinical.BCRtab.all$clinical_drug_acc  %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

#Clinical indexed data
#In this example we will fetch clinical indexed data (same as showed in the data portal).

clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
clinical %>%
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)
#????????ǰ??query??
query <- GDCquery(project = "TCGA-LUAD",
                  
                  data.category = "Transcriptome Profiling",
                  
                  data.type = "Gene Expression Quantification",
                  
                  workflow.type = "STAR - Counts")
data <- getResults(query)
View(data)
#???????????У??????? harmonized ???ݿ????????????? GBM ?????ļ׻???????

query <- GDCquery(
  project = c("TCGA-LUAD"),
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  sample.type = "Recurrent Tumor"
)

data <- getResults(query)
View(data) 

#?? harmonized ???ݿ??зֱ???ȡ???а????׻??????ݣ?450k???ͱ??????ݵ????ٰ???????Ϣ

query.met <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450")
)
query.exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
# ??ȡ?ȼ????׻????ּ?????????????
common.patients <- intersect(
  substr(getResults(query.met, cols = "cases"), 1, 12),
  substr(getResults(query.exp, cols = "cases"), 1, 12)
)

length(common.patients)
#[1] 459

#???? 459 ???????????????? barcode ??????ֻѡ?񲿷???????????
save(query.exp,file="query.exp.RData")
query.exp$results
write.table(unclass(query.exp$results),file="TCGA-LUAD_sample_barcode_info.txt",quote = F,sep="\t",row.names=T)
#???غ???ȡ????????
GDCdownload(query.exp)
exps<-GDCprepare(query = query.exp)
exp.matrix<-SummarizedExperiment::assay(exps)#assay????��Դ??SummarizedExperiment??
expdat <- GDCprepare(
  query = query.exp,
  save = TRUE, 
  save.filename = "exp.rda"
)

###########################################
##miRNA Expression Quantification
library(TCGAbiolinks)
query.mirna <- GDCquery(
  project = "TARGET-LUAD", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
  data.type = "miRNA Expression Quantification"
)
GDCdownload(query.mirna)
mirna <- GDCprepare(
  query = query.mirna,
  save = TRUE, 
  save.filename = "mirna.rda"
)
###############################################

query.met <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  barcode = common.patients[1:5]
)

#############################################################
list.files()
###????id mapping
#??ȡjson?ļ?????ȡ????fromJSON��??R?Դ?R??jsonlite????
meta_dt <- jsonlite::fromJSON("TCGA-LUAD.metadata.cart.json")
#?鿴????list?е??????еĵ?һ??list(???????ݿ?):
meta_dt[[3]][[1]]
id <- meta_dt$associated_entities
id[[1]]
id[[1]][,1] #??????????Ҫ??????
head(meta_dt[[4]])#??????????Ҫ??��ƥ?????ļ???
ids <- sapply(id,function(x){x[,1]}) #??ȡ????????id??????һ??????????id??????��????ǰ??lapply????????lapply???ص???list??
ids2 <- data.frame(file_name = meta_dt$file_name,
                   id = ids) #????һ????????????id???ļ?????��?ж?Ӧ??ϵ?????ݿ????ļ?????Ϊ???????????ļ???
head(ids2)
#################################################################
#??��??ȡ????tsv??ʽ??׺???ļ?????
exp_all <- dir("GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification",
               pattern = "*.rna_seq.augmented_star_gene_counts.tsv$",#*??ʾ????ǰ׺??$??ʾ?̶???׺
               recursive = T)
head(exp_all)
#?????Զ??庯??,??????��????????tcv?ļ???unstranded?У?counts??????fpkm_unstranded(fpkm)??
exp_dt <- function(x){
  result <- data.table::fread(file.path("GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/",x),
                              select = c("fpkm_unstranded"))
  return(result)
}
exp <- lapply(exp_all,exp_dt)#??ȡ????.tsv??unstranded?У???exp_allת??Ϊlist??????list?е?ÿһ??Ԫ?ض?Ӧ?ú???exp_dt??
exp <- do.call(cbind,exp)#??exp???кϲ???????listת??Ϊdata.table
exp[1:6,1:6]

##ע?⣬data.table::fread??ʽ?????ļ???????ָ???????????Բ?ʹ?????ַ???
x1 <- data.table::fread("GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/0b73f593-b6e8-4897-9b9b-f09e65e436b7/30a6ea89-8eda-47cf-b336-e070fcd08740.rna_seq.augmented_star_gene_counts.tsv",
                        select = c("gene_id"))
head(x1)
#????????(gene_id),???б?????????gene_id????ͬ?ģ?
exp <- as.data.frame(exp)#data.table??ʽ?????Զ?????????????????ת??Ϊ???ݿ?
rownames(exp) <- x1$gene_id
exp[1:8,1:4]
#?á?/????һ???ַ??????֣??ѱ??????????ļ???(exp_all)?ĺ????κ?file_name??ƥ?䣺
library(stringr)
exp_all2 <- str_split(exp_all,"/",simplify = T)#??/?????֣?simplify = T?ѽ??????سɾ???
head(exp_all2)
#ȡ?????ڶ??У?
exp_all2 <- exp_all2[,2]
#?ðٷְ?in?????Ƿ????ֺ????ļ???????ȫ??????file_name??(??ʱ˳????ͬ?????˲?????==)
exp_all2 %in% ids2$file_name#ȫ??????TRUE????????????

#??match????????file_name??˳????exp_all2һ?£???Ϊ???Ǻϲ??ı????????ǰ???exp_all2??˳?????кϲ??ģ?
ids3 <- ids2[match(exp_all2,ids2$file_name),]
#?鿴??????˳????
head(ids3$file_name)
head(exp_all2)
identical(ids3$file_name,exp_all2)#????TRUE??ȷ??????
#?޸?????Ϊ????????
colnames(exp) <- ids3$id
exp[1:8,1:2]#???????????ɺϲ???
dim(exp)

##??????һ??tsvȡgene_name??һ?У??????ǵ?ת??Ϊ???ݿ򣬲?Ȼ?????ϲ??ᶪʧexp??????????
list <- as.data.frame(data.table::fread("GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/0b73f593-b6e8-4897-9b9b-f09e65e436b7/30a6ea89-8eda-47cf-b336-e070fcd08740.rna_seq.augmented_star_gene_counts.tsv",
                                        select = c("gene_name")))
exp2 <- cbind(list,exp)
exp2[1:10,1:3]#????gene_name?еľ???over??

write.table(exp,file="TCGA-LUAD_598sample_exp_FPKM.txt",quote = F,sep="\t",row.names=T)
write.table(exp2,file="TCGA-LUAD_598sample_exp_FPKM.addGene.txt",quote = F,sep="\t",row.names=T)
#????ֻ??????50%?????ж??????Ļ?????
#gene filter
#exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ] # only keep genes express in half of samples
#dim(exp)
#[1] 31784   598
#write.table(exp,file="TCGA-LUAD_598sample_exp_count.geneFilter.txt",quote = F,sep="\t",row.names=T)
########
###???ֱ????ͷǱ???????????
#####??????(ֱ????%in%?жϲ?ɸѡ????)??
##?????????????????Ȱѱ?????????gene-name??gene_type???ߺϲ???ͨ??gene_typeɸѡ???????Ѳ???Ҫ????ȥ????
list22 <- as.data.frame(data.table::fread("GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/0b73f593-b6e8-4897-9b9b-f09e65e436b7/30a6ea89-8eda-47cf-b336-e070fcd08740.rna_seq.augmented_star_gene_counts.tsv",
                                          select = c("gene_name","gene_type")))
exp3 <- cbind(list22,exp) #??gene_type?кϲ???��
exp3[1:8,1:3]
#????ֱ??ɸѡ???ɣ?
mRNA_exp2 <- exp3[exp3$gene_type %in% c("protein_coding"),]
lncRNA_exp2 <- exp3[exp3$gene_type %in% c("lncRNA"),]

#????ȥ??????ɸѡ??gene_type??һ?У?
mRNA_exp2 <- mRNA_exp2[,-2]
lncRNA_exp2 <- lncRNA_exp2[,-2]

write.table(mRNA_exp2,file="TCGA-LUAD_598sample_mRNAexp_FPKM.addGene.txt",quote = F,sep="\t",row.names=T)

write.table(lncRNA_exp2,file="TCGA-LUAD_598sample_lncRNAexp_FPKM.addGene.txt",quote = F,sep="\t",row.names=T)

load("C:\\??Ŀ\\TCGA_PAAD\\output_expr\\TCGA-PAAD_mrna_expr_tpm.rdata")
mrna_expr_tpm[1:4,1:4];dim(mrna_expr_tpm)
write.table(mrna_expr_tpm, file = "mrna_expr_tpm.txt",sep = "\t",
            row.names = T,col.names = NA,quote = F)
#group information
library(stringr)
table(str_sub(colnames(exp),14,15))
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)
# normal  tumor 
#     59    539

#################################
##??Ϊnormal??tumor????
sample <- colnames(mRNA_exp2)

normal <- c()
tumor <- c()

for (i in 1:length(sample)){
  if((substring(colnames(mRNA_exp2)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
    normal <- append(normal,sample[i])
  } else {
    tumor <- append(tumor,sample[i])
  }
}

tumor_matrix <- mRNA_exp2[,tumor]
normal_matrix <- mRNA_exp2[,normal]

write.table(tumor_matrix,file="TCGA-LUAD_598sample_tumor_mRNAexp_FPKM.addGene.txt",
            quote = F,sep="\t",row.names=T)
write.table(normal_matrix,file="TCGA-LUAD_598sample_normal_mRNAexp_FPKM.addGene.txt",
            quote = F,sep="\t",row.names=T)
#
sample <- colnames(lncRNA_exp2)

normal <- c()
tumor <- c()

for (i in 1:length(sample)){
  if((substring(colnames(lncRNA_exp2)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
    normal <- append(normal,sample[i])
  } else {
    tumor <- append(tumor,sample[i])
  }
}

tumor_matrix <- lncRNA_exp2[,tumor]
normal_matrix <- lncRNA_exp2[,normal]

write.table(tumor_matrix,file="TCGA-LUAD_598sample_tumor_lncRNAexp_FPKM.addGene.txt",quote = F,sep="\t",row.names=T)
write.table(normal_matrix,file="TCGA-LUAD_598sample_normal_lncRNAexp_FPKM.addGene.txt",quote = F,sep="\t",row.names=T)

#????????
project="TCGA-LUAD"
if(!dir.exists("data"))dir.create("data")
save(exp,Group,project,file = paste0("data/",project,"_TCGA_LUAD_gdc.Rdata"))


# ??ȡ?ض?????
patient_id<-read.table("patient_id.txt",header=T,sep='\t',row.names=1)
sub_group <- exp[,colnames(exp) %in% rownames(patient_id)]

####??????50%?????ж??????Ļ??򣬴???ȥ??????��ȫ??????1??count???Ļ???
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ] # only keep genes express in half of samples
dim(exp)
exp_1=exp[apply(exp,1,function(x){XXXX<-FALSE;if(any(x>1)){XXXX<-TRUE};return(XXXX)}),];dim(exp_1) #ȥ??????��ȫ??????1??count???Ļ???
dim(exp_1)

getTCGAexpr <- function(project){
  
  dir.create("output_expr")
  
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts"
  )
  
  GDCprepare(query,save = T,save.filename = paste0("output_expr/",project,"_expr.rdata"))
  
  load(file = paste0("output_expr/",project,"_expr.rdata"))
  
  se <- data
  
  clin_info <- as.data.frame(colData(se))
  save(clin_info, file = paste0("output_expr/",project, "_clinical.rdata"))
  
  rowdata <- rowData(se)
  
  se_mrna <- se[rowdata$gene_type == "protein_coding",]
  se_lnc <- se[rowdata$gene_type == "lncRNA"]
  
  
  expr_counts_mrna <- assay(se_mrna,"unstranded")
  
  expr_tpm_mrna <- assay(se_mrna,"tpm_unstrand")
  
  expr_fpkm_mrna <- assay(se_mrna,"fpkm_unstrand")
  
  expr_counts_lnc <- assay(se_lnc,"unstranded")
  
  expr_tpm_lnc <- assay(se_lnc,"tpm_unstrand")
  
  expr_fpkm_lnc <- assay(se_lnc,"fpkm_unstrand")
  
  symbol_mrna <- rowData(se_mrna)$gene_name
  
  symbol_lnc <- rowData(se_lnc)$gene_name
  
  
  mrna_expr_counts <- cbind(data.frame(symbol_mrna),as.data.frame(expr_counts_mrna))
  mrna_expr_counts <- mrna_expr_counts %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_counts, file = paste0("output_expr/",project, "_mrna_expr_counts.rdata"))
  write.table(mrna_expr_counts, file = "mrna_expr_counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  
  
  mrna_expr_tpm <- cbind(data.frame(symbol_mrna),as.data.frame(expr_tpm_mrna))
  mrna_expr_tpm <- mrna_expr_tpm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_tpm, file = paste0("output_expr/",project, "_mrna_expr_tpm.rdata"))
  
  
  mrna_expr_fpkm <- cbind(data.frame(symbol_mrna),as.data.frame(expr_fpkm_mrna))
  mrna_expr_fpkm <- mrna_expr_fpkm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_mrna,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_mrna") %>% 
    as.data.frame()
  save(mrna_expr_fpkm, file = paste0("output_expr/",project, "_mrna_expr_fpkm.rdata"))
  
  
  lncrna_expr_counts <- cbind(data.frame(symbol_lnc),as.data.frame(expr_counts_lnc))
  lncrna_expr_counts <- lncrna_expr_counts %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_counts, file = paste0("output_expr/",project, "_lncrna_expr_counts.rdata"))
  
  
  lncrna_expr_tpm <- cbind(data.frame(symbol_lnc),as.data.frame(expr_tpm_lnc))
  lncrna_expr_tpm <- lncrna_expr_tpm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_tpm, file = paste0("output_expr/",project, "_lncrna_expr_tpm.rdata"))
  
  
  lncrna_expr_fpkm <- cbind(data.frame(symbol_lnc),as.data.frame(expr_fpkm_lnc))
  lncrna_expr_fpkm <- lncrna_expr_fpkm %>% 
    as_tibble() %>%
    mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
    arrange(desc(meanrow)) %>% 
    distinct(symbol_lnc,.keep_all=T) %>% 
    select(-meanrow) %>% 
    column_to_rownames(var = "symbol_lnc") %>% 
    as.data.frame()
  save(lncrna_expr_fpkm, file = paste0("output_expr/",project, "_lncrna_expr_fpkm.rdata"))
  
}

getTCGAexpr("TCGA-PAAD")

##########################################################################
#???????ֱ???��??count??TPM?????????????ͼ???????
list.files(path="./output_expr")
#[1] "TCGA-PAAD_clinical.rdata"           "TCGA-PAAD_expr.rdata"              
#[3] "TCGA-PAAD_lncrna_expr_counts.rdata" "TCGA-PAAD_lncrna_expr_fpkm.rdata"  
#[5] "TCGA-PAAD_lncrna_expr_tpm.rdata"    "TCGA-PAAD_mrna_expr_counts.rdata"  
#[7] "TCGA-PAAD_mrna_expr_fpkm.rdata"     "TCGA-PAAD_mrna_expr_tpm.rdata" 

getTCGAexprtumorcontrol <- function(project){
  
  load(paste0("./output_expr/",project, "_mrna_expr_counts.rdata"))
  write.table(mrna_expr_counts, file = "mrna_expr_counts.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(mrna_expr_counts)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(mrna_expr_counts)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- mrna_expr_counts[,tumor]
  normal_matrix <- mrna_expr_counts[,normal]
  write.table(tumor_matrix, file = "mrna_expr_counts.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "mrna_expr_counts.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  
  load(paste0("./output_expr/",project, "_mrna_expr_tpm.rdata"))
  write.table(mrna_expr_tpm, file = "mrna_expr_tpm.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(mrna_expr_tpm)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(mrna_expr_tpm)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- mrna_expr_tpm[,tumor]
  normal_matrix <- mrna_expr_tpm[,normal]
  write.table(tumor_matrix, file = "mrna_expr_tpm.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "mrna_expr_tpm.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  

  load(paste0("./output_expr/",project, "_mrna_expr_fpkm.rdata"))
  write.table(mrna_expr_fpkm, file = "mrna_expr_fpkm.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(mrna_expr_fpkm)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(mrna_expr_fpkm)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- mrna_expr_fpkm[,tumor]
  normal_matrix <- mrna_expr_fpkm[,normal]
  write.table(tumor_matrix, file = "mrna_expr_fpkm.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "mrna_expr_fpkm.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  
  load(paste0("./output_expr/",project, "_lncrna_expr_counts.rdata"))
  write.table(lncrna_expr_counts, file = "lncrna_expr_counts.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(lncrna_expr_counts)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(lncrna_expr_counts)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- lncrna_expr_counts[,tumor]
  normal_matrix <- lncrna_expr_counts[,normal]
  write.table(tumor_matrix, file = "lncrna_expr_counts.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "lncrna_expr_counts.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  

  load(paste0("./output_expr/",project, "_lncrna_expr_tpm.rdata"))
  write.table(lncrna_expr_tpm, file = "lncrna_expr_tpm.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(lncrna_expr_tpm)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(lncrna_expr_tpm)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- lncrna_expr_tpm[,tumor]
  normal_matrix <- lncrna_expr_tpm[,normal]
  write.table(tumor_matrix, file = "lncrna_expr_tpm.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "lncrna_expr_tpm.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  
  load(paste0("./output_expr/",project, "_lncrna_expr_fpkm.rdata"))
  write.table(lncrna_expr_fpkm, file = "lncrna_expr_fpkm.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
  sample <- colnames(lncrna_expr_fpkm)
  normal <- c()
  tumor <- c()
  for (i in 1:length(sample)){
    if((substring(colnames(lncrna_expr_fpkm)[i],14,15)>10)){    #14??15λ?ô???10??Ϊnormal????
      normal <- append(normal,sample[i])
    } else {
      tumor <- append(tumor,sample[i])
    }
  }
  tumor_matrix <- lncrna_expr_fpkm[,tumor]
  normal_matrix <- lncrna_expr_fpkm[,normal]
  write.table(tumor_matrix, file = "lncrna_expr_fpkm.tumor.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  write.table(normal_matrix, file = "lncrna_expr_fpkm.normal.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  
}

getTCGAexprtumorcontrol("TCGA-PAAD")

