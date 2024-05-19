#Linux shell codes

#Motif Analysis
findMotifsGenome.pl Res_Vsh183_peaks.bed hg38

#liftover (chain downloaded from ucsc)
liftOver H3K4me1_HUH7.enhancer.bedgraph hg19ToHg38.over.chain H3K4me1_HUH7.enhancer.hg38.bedgraph H3K4me1_HUH7.unmatched.bedgraph