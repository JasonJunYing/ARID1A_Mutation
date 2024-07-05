#Linux shell codes

### liftover (chain downloaded from https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
liftOver GSE89212_enhancers.bed hg19ToHg38.over.chain GSE89212_enhancers_hg38.bed GSE89212_enhancers_unmatched.bed

### Enhancer annotation
bedtools intersect -a Res_Vsh183_allpeaks.bed -b GSE89212_enhancers_hg38.bed -wa > Res_Vsh_V.Enhancer.bed
