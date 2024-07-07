# ARID1A_Mutation
## Deciphering the regulatory mechanisms and biological implications of ARID1A C-terminal missense mutations in cancer

We combined ATAC-seq and ChIP-seq to investigate the influence of ARID1A mutation on the regulatory landscape in Huh7 cells using differential analysis and motif enrichment. Raw sequencing data and bed files were uploaded to GEO database (GSE254864 and GSE254865). The procedures of preprocessing and peak calling can be found in: [Preprocessing_and_peak_calling.sh](./Preprocessing_and_peak_calling.sh). 
1. ATAC-seq analysis : [Diffbind_ATAC.R](./Diffbind_ATAC.R)
(1) Calculate binding matrix

(2) Normalization (method: libsize as default)
(3) Differential enrichment analysis (Figure 6B)
(4) Peak region annotation
(5) Define rescued region
3. ChIP-seq analysis : [Diffbind_CHIP.R](./Diffbind_CHIP.R)
   (1) Calculate binding matrix
   (2) Normalization (method: libsize as default)
   (3) Differential enrichment analysis
   (4) Peak region annotation (Figure 6C)
   (5) Define rescued region (Figure 6D)
   (6) Enhancer annotation (using the published H3K27ac and H3K4me3 enhancer regions of Huh7 cells (GSE89212) [Enhancer_annotation.sh](Enhancer_annotation.sh))
   (7) Motif enrichment analysis (Figure 6K) [Motif_Analysis.sh](Motif_Analysis.sh);[MotifAnalysis_Resc.R](MotifAnalysis_Resc.R)

Furthermore, we extracted the RNA-seq data of liver cancer patients with different types of ARID1A mutations from the TCGA-LIHC database, and investigated the expression of the genes potentially regulated by ARID1A.
3. TCGA analysis: [TCGA_analysis.R](TCGA_analysis.R)
   (1) Retrieve SNP information (maf)
   (2) Control sample extraction (patients without ARID1A mutation)
   (3) Retrieve RNA-seq expression matrices and convert gene symbols
   (4) CPM normalization
   (5) Select Cluster 2 genes (as shown in Figure 6E) and calculate mean expression
   (6) Comparison between the mutated and control samples (Figure 6M)
