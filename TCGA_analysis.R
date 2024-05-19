setwd('./TCGA')
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(pheatmap)
library(Matrix)
library(biomaRt)
library(limma)
library(edgeR)
library(ggpubr)
library(ggplot2)

### RETRIEVE SNP ###

data_type <- "Masked Somatic Mutation"
data_category <- "Simple Nucleotide Variation"
cancer_type <- "TCGA-LIHC"

for (l in cancer_type){
  query_maf <- GDCquery(project = l, 
                        data.category = data_category, 
                        data.type =  data_type, 
  )
  
  GDCdownload(query_maf, method = "api")
  snp <- GDCprepare(query = query_maf)
  write.csv(snp,paste0(l,".all.snp.csv"))
  snp1 <- snp[snp$Hugo_Symbol=="ARID1A",]
  write.csv(snp1,paste0(l,".arid1a.snp.csv"))
}

#Set patients without ARID1A mutation as ctrl
snp <- read.csv('./TCGA/TCGA-LIHC.all.snp.csv',row.names = 1)
res <- snp[1,]
for(i in unique(snp$X1)){
  tmp <- subset.data.frame(snp,snp$X1==i)
  if(!any(grepl('ARID1A',tmp$Hugo_Symbol))){
    res <- rbind(res,tmp[1,])
  }
}
res <- res[-1,]
write.csv(res,'ctrls_LIHC.csv')

### RETRIEVE mRNA###

mutlist <- c("TCGA-ZS-A9CD",	"TCGA-UB-A7MD", "TCGA-BW-A5NQ")
ctrlinfo <- read.csv('ctrls_LIHC.csv')
ctrlinfo$barcode <- substring(ctrlinfo$Tumor_Sample_Barcode,1,12)
patient_list <- c(mutlist,ctrlinfo$barcode)

clinical <- read.csv('TCGA-LIHC_clinical.csv',row.names = 2) #Downloaded from TCGA
clinical <- clinical[patient_list,]

clinical <- clinical[!is.na(clinical$age_at_index),]
clinical <- clinical[clinical$age_at_index>=63&clinical$age_at_index<=73,] 
clinical <- clinical[!is.na(clinical$ajcc_pathologic_stage),]
clinical <- clinical[clinical$ajcc_pathologic_stage!="Stage IIIA"&clinical$ajcc_pathologic_stage!="Stage IIIB",]
clinical <- clinical[clinical$gender=='male',]
clinical$ajcc_stage <- clinical$ajcc_pathologic_stage

patient_list <- intersect(patient_list,rownames(clinical))

data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "STAR - Counts"
cancer_type <- "TCGA-LIHC"
barcode_type <- patient_list

query_TranscriptomeCounts <- GDCquery(project = cancer_type, 
                                      data.category = data_category, 
                                      data.type =  data_type, 
                                      workflow.type = workflow_type,
                                      barcode = barcode_type,
)

GDCdownload(query_TranscriptomeCounts, method = "api")
expdat <- GDCprepare(query = query_TranscriptomeCounts)
count_matrix=assay(expdat)
counts <- as.data.frame(count_matrix)
counts$EnsemblID <- strsplit2(rownames(count_matrix),'\\.')[,1]
counts <- counts[!duplicated(counts$EnsemblID),]

#Use biomart to convert gene symbol
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",
                 host = "https://jan2020.archive.ensembl.org/")
res<-getBM(attributes = c("ensembl_gene_id","external_gene_name"),
           filters = "ensembl_gene_id",
           values = counts$EnsemblID,
           mart = human)

#Continue
counts <- na.omit(counts[match(res$ensembl_gene_id,counts$EnsemblID),])
counts$Symbol <- res$external_gene_name
counts <- counts[!duplicated(counts$Symbol),]
rownames(counts)<-counts$Symbol
counts <- counts[,1:(ncol(counts)-2)]
for(l in 1:ncol(counts)){counts[,l]<-as.numeric(counts[,l])}
counts <- as.matrix(counts)
write.csv(counts,'./TCGA-LIHC.rawcnts.csv')
cpmcounts <- cpm(counts,log = T)
write.csv(cpmcounts,'./TCGA-LIHC.cpmcnts.csv')

#Select c2 genes
genelist <- read.csv('./RNA_genelist.csv')
c2counts <- na.omit(cpmcounts[match(genelist$cluster2,rownames(cpmcounts)),])
write.csv(c2counts,'TCGA-LIHC_c2.cpmcnts.csv')

### CHECK EXP ###
c2counts <- read.csv('./TCGA-LIHC_c2.cpmcnts.csv',row.names = 1)
MeanExp <- colMeans(c2counts)
colmean <- as.data.frame(MeanExp)
colmean$Group <- c(rep("Mutant",length(mutlist)),
                   rep("Ctrl",ncol(c2counts)-length(mutlist)))
colmean$PatientID <- gsub("\\.","-",substring(rownames(colmean),1,12))
rownames(colmean)<-1:nrow(colmean)
clinical$PatientID <- rownames(clinical)
exp <- merge(colmean,clinical,by='PatientID')

q <- ggplot(exp,aes(x=Group,y=MeanExp))+
  geom_violin()+
  geom_jitter(width = 0.2,aes(col=ajcc_stage))+
  ggtitle(paste0("TCGA_LIHC"," n=",nrow(exp)))+
  ylab("MeanExp (C2)")+
  stat_compare_means(method = "t.test",
                     label.x = 1.4,
                     label.y = min(exp$MeanExp))+
  theme_bw()

