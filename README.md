# ARID1A_Mutation
Source codes for: Deciphering the regulatory mechanisms and biological implications of ARID1A C-terminal missense mutations in cancer

We combined ATAC-seq and ChIP-seq to investigate the influence of ARID1D mutation on the regulatory landscape in HEK293T cells using differential analysis and motif enrichment. Raw sequencing data and bed files was uploaded to GEO database. The analysis pipelines of the ATAC-seq and the ChIP-seq data can be found in: [Diffbind_ATAC.R](./Diffbind_ATAC.R), [Diffbind_CHIP.R](./Diffbind_CHIP.R), which were used to generate Figure 5B and 5C. The procedures of preprocessing, peak calling can be found in: [Linux shell codes](./Linux shell codes.sh). 

Furthermore, we extracted the RNA-seq data of liver cancer patients with different types of ARID1A mutations from the TCGA-LIHC database, and investigated the expression of the genes potentially regulated by ARID1A [TCGA_analysis.R](TCGA_analysis.R)
