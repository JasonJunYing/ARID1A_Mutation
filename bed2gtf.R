# R script to convert bed to gtf
## xw251, 10Mar2022
## input: bed file
## output: gtf file

bedfile = "all_merged.bed"
BED = read.table(bedfile, sep="\t", header=FALSE)
BED[,1]<-as.character(BED[,1])
GTF = cbind(BED[,1], "Merged", "PEAKS", BED[,2], BED[,3], ".", "+", ".", paste0('peak_id "peak_', 1:nrow(BED), '"'))
write.table(GTF, "final.gtf", quote=FALSE, row.names=F, col.names=F, sep="\t")

