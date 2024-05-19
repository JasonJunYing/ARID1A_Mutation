library(DiffBind)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

meta<-read.csv('./CHIP/meta.csv',row.names = 1)
design<-model.matrix(~0+Group,data = meta)

#ReadIn
samples <- read.table('./CHIP/SampleSheet_CHIP.txt',header = T)
chip2 <- dba(sampleSheet=samples)
dba.plotHeatmap(chip2) # UnNormed

#Count
chip2 <- dba.count(chip2) #to calculate binding matrix

info <- dba.show(chip2)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, 
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

chip2 <- dba.normalize(chip2) # Normalization (default:libsize)

norm <- dba.normalize(chip2, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
dba.plotHeatmap(chip2) # Normed

profiles <- dba.plotProfile(chip2,merge = NULL)
dba.plotProfile(profiles)
saveRDS(chip2,'./CHIP/CHIP_Uncontrasted.rds')

#Make contrast
names(meta)[2]<-"Condition"
design <- model.matrix(~0+Condition,data = meta)
chip2 <- readRDS('CHIP_Uncontrasted.rds')
chip2 <- dba.contrast(chip2, minMembers = 2,
                          reorderMeta=list(Condition="WT-shEV"))
chip2 <- dba.analyze(chip2,method = DBA_ALL_METHODS)
edgeRlist <- dba.analyze(chip2,bRetrieveAnalysis = DBA_EDGER_GLM)
dba.show(chip2, bContrasts=TRUE)
plot(chip2, contrast=2)

chip2.V.WT.all <- dba.report(chip2,contrast = 9,bUsePval = FALSE,th=1,method = DBA_EDGER)
chip2.Vsh.V.all <- dba.report(chip2,contrast = 2,bUsePval = FALSE,th=1,method = DBA_EDGER)

#Annotation
peakanno <- annotatePeak(chip2.V.WT.all,
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
write.csv(ann.all.V.WT,'./CHIP/Final_ann_V_WT.csv')

peakanno <- annotatePeak(chip2.Vsh.V.all,
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
ann.all.Vsh.V <- as.data.frame(peakanno)
write.csv(ann.all.Vsh.V,'./CHIP/Final_ann_Vsh_V.csv')

#Define rescued regions
VWT <- read.csv('./CHIP/Final_ann_V_WT.csv')
VWT_Down <- as.data.frame(VWT[VWT$Fold<(-log2(1.5)),])
VWT_Down_TSS <- VWT_Down[abs(VWT_Down$distanceToTSS)<2000,]
VWT_Up <- as.data.frame(VWT[VWT$Fold>log2(1.5),])
VWT_Up_TSS <- VWT_Up[abs(VWT_Up$distanceToTSS)<2000,]
VWT_NS <- VWT[abs(VWT$Fold)<log2(1.5),]
VWT_NS_TSS <- VWT_NS[abs(VWT_NS$distanceToTSS)<2000,]

VshV <- read.csv('./CHIP/Final_ann_Vsh_V.csv')
VshV_Up <- as.data.frame(VshV[VshV$Fold>log2(1.5),])
VshV_Down <- as.data.frame(VshV[VshV$Fold<(-log2(1.5)),])
VshV_Up_TSS <- VshV_Up[abs(VshV_Up$distanceToTSS)<2000,]
VshV_Down_TSS <- VshV_Down[abs(VshV_Down$distanceToTSS)<2000,]
VshV_NS <- VshV[abs(VshV$Fold)<log2(1.5),]
VshV_NS_TSS <- VshV_NS[abs(VshV_NS$distanceToTSS)<2000,]

Res_Vsh_V_all <- VWT_Down[match(VshV_Up$start,VWT_Down$start),]
Res_Vsh_V <- VWT_Down_TSS[match(VshV_Up_TSS$start,VWT_Down_TSS$start),]
Res_Vsh_V_all <- Res_Vsh_V_all[!is.na(Res_Vsh_V_all$seqnames),]
Res_Vsh_V <- Res_Vsh_V[!is.na(Res_Vsh_V$seqnames),]
Res_Vsh_V <- Res_Vsh_V[abs(Res_Vsh_V$distanceToTSS)<2000,]
write.csv(Res_Vsh_V_all,'./CHIP/Res_Vsh183_allpeaks.csv')

#Make bed files
Res_Vsh_V <- read.csv('./CHIP/Res_Vsh183_peaks.csv',row.names = 1)
ResBed <- Res_Vsh_V[,1:5]
names(ResBed)[4]<-"PeakID"
ResBed$PeakID <- rownames(ResBed)
write.table(ResBed, './CHIP/Res_Vsh183_peaks.bed', sep = '\t', quote = FALSE,row.names = F)

#Annotate Rescued regions
chip2.Vsh183.Res <- chip2.V.WT.all[match(Res_Vsh_V$Fold,chip2.V.WT.all$Fold),]
peakanno <- annotatePeak(chip2.Vsh183.Res,
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
ann.Res <- as.data.frame(peakanno)
write.csv(ann.Res,'./CHIP/Final_ann_Res.csv')