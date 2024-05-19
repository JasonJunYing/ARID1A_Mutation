library(DiffBind)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

meta<-read.csv('./ATAC/meta.csv',row.names = 1)
design<-model.matrix(~0+Group,data = meta)

#ReadIn
atac <- dba(sampleSheet='./ATAC/SampleSheet_ATAC.csv')
dba.plotHeatmap(atac) # UnNormed

#Count
atac <- dba.count(atac) #to calculate binding matrix

info <- dba.show(atac)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

atac <- dba.normalize(atac) # Normalization (default:libsize)

norm <- dba.normalize(atac, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
dba.plotHeatmap(atac) # Normed

profiles <- dba.plotProfile(atac,merge = NULL)
dba.plotProfile(profiles)
saveRDS(atac,'./ATAC/atac_Uncontrasted.rds')

#Make contrast
names(meta)[2]<-"Condition"
design <- model.matrix(~0+Condition,data = meta)
atac <- readRDS('./ATAC/atac_Uncontrasted.rds')
atac <- dba.contrast(atac, minMembers = 2,
                          reorderMeta=list(Condition="WT"))
atac <- dba.analyze(atac,method = DBA_ALL_METHODS)
edgeRlist <- dba.analyze(atac,bRetrieveAnalysis = DBA_EDGER_GLM)
dba.show(atac, bContrasts=TRUE)
plot(atac, contrast=2)

#Annotation
atac.WT.V.all <- dba.report(atac,contrast = 4,bUsePval = FALSE,th=1,method = DBA_EDGER) #All
atac.Vsh.V.all <- dba.report(atac,contrast = 6,bUsePval = FALSE,th=1,method = DBA_EDGER) #All

peakanno <- annotatePeak(atac.WT.V.all,
                         tssRegion = c(-2000, 2000), 
                         TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         level = "transcript", 
                         assignGenomicAnnotation = TRUE, 
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                         annoDb = "org.Hs.eg.db", 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000, 
                         sameStrand = FALSE, 
                         ignoreOverlap = FALSE, 
                         ignoreUpstream = FALSE, 
                         ignoreDownstream = FALSE, 
                         overlap = "TSS", 
                         verbose = TRUE)
ann.all.WT.V <- as.data.frame(peakanno)
write.csv(ann.all.WT.V,'./ATAC/Final_ann_WT_V.csv')

peakanno <- annotatePeak(atac.Vsh.V.all,
                         tssRegion = c(-2000, 2000), 
                         TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         level = "transcript", 
                         assignGenomicAnnotation = TRUE, 
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                         annoDb = "org.Hs.eg.db", 
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 5000, 
                         sameStrand = FALSE, 
                         ignoreOverlap = FALSE, 
                         ignoreUpstream = FALSE, 
                         ignoreDownstream = FALSE, 
                         overlap = "TSS", 
                         verbose = TRUE)
ann.all.V.WT <- as.data.frame(peakanno)
write.csv(ann.all.V.WT,'./ATAC/Final_ann_Vsh_V.csv')

VWT <- read.csv('./ATAC/Final_ann_WT_V.csv')
VWT$Fold = (-VWT$Fold) # trans to V.WT
VWT_Down <- as.data.frame(VWT[VWT$Fold<(-log2(1.5)),])
VWT_Down_TSS <- VWT_Down[abs(VWT_Down$distanceToTSS)<2000,]
VWT_Up <- as.data.frame(VWT[VWT$Fold>log2(1.5),])
VWT_Up_TSS <- VWT_Up[abs(VWT_Up$distanceToTSS)<2000,]
VWT_NS <- VWT[abs(VWT$Fold)<log2(1.5),]
VWT_NS_TSS <- VWT_NS[abs(VWT_NS$distanceToTSS)<2000,]

VshV <- read.csv('./ATAC/Final_ann_Vsh_V.csv')
VshV_Up <- as.data.frame(VshV[VshV$Fold>log2(1.5),])
VshV_Down <- as.data.frame(VshV[VshV$Fold<(-log2(1.5)),])
VshV_Up_TSS <- VshV_Up[abs(VshV_Up$distanceToTSS)<2000,]
VshV_Down_TSS <- VshV_Down[abs(VshV_Down$distanceToTSS)<2000,]
VshV_NS <- VshV[abs(VshV$Fold)<log2(1.5),]
VshV_NS_TSS <- VshV_NS[abs(VshV_NS$distanceToTSS)<2000,]

Res_Vsh_V_all <- VWT_Down[match(VshV_Up$start,VWT_Down$start),]
Res_Vsh_V_all <- Res_Vsh_V_all[!is.na(Res_Vsh_V_all$seqnames),]
Res_Vsh_V <- VWT_Down_TSS[match(VshV_Up_TSS$start,VWT_Down_TSS$start),]
Res_Vsh_V <- Res_Vsh_V[!is.na(Res_Vsh_V$seqnames),]
Res_Vsh_V <- Res_Vsh_V[abs(Res_Vsh_V$distanceToTSS)<2000,]
write.csv(Res_Vsh_V,'./ATAC/Res_Vsh_V.csv')
