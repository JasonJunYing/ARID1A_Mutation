#Linux shell codes

### Raw data processing

#Preprocess
FILE1=*_R1.fq
FILE2=*_R2.fq

fn=$(basename "$FILE1")
fn="${fn%_R1.fq}"

mkdir fastQC
echo "Running fastQC for $fn ..."
fastqc -o fastQC -f fastq $FILE1

echo "trimming"
trimmomatic PE -phred33 $FILE1 $FILE2 $fn\.R1_trimmed.fq $fn\.R1_unpairtrimmed.fq $fn\.R2_trimmed.fq $fn\.R2_unpairtrimmed.fq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

fastqc -o fastQC -f fastq $fn\.R1_trimmed.fq

echo "Running alignment with bowtie2 ... "
bowtie2 -N 1 -k 1 -p 20 -x ~/refgen/hg38_bt2/hg38 -1 $fn\.R1_trimmed.fq -2 $fn\.R2_trimmed.fq -S $fn\.sam

echo "Converting to bam..."
samtools view -@ 8 -b -T ~/refgen/hg38.fa -S $fn\.sam > $fn\.bam
samtools sort -@ 8 $fn\.bam > $fn\.sorted.bam

echo "Running Picard tools for duplicates removing..."
java -jar ~/Apps/picard.jar MarkDuplicates I=$fn\.sorted.bam O=$fn\.dupremoved.bam M=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=./tmp

echo "Remove low quality reads..."
samtools view -@ 8 -h $fn\.dupremoved.bam | awk 'substr($0,1,1) == "@" || $2 !=4 || $5>=30 {print}' | samtools view -@ 8 -bS - > $fn\.filtered.bam
samtools view -@ 8 $fn\.filtered.bam | egrep -v "chrM|chrUn" | samtools view -@ 8 -bT ~/refgen/hg38.fa - -o $fn\.processed.bam

#Macs2 call peaks
for FILE in *.processed.bam
	do
		fn=$(basename "$FILE")
		fn="${fn%.processed.bam}"

		echo "MACS2 peak calling..."
		#q=0.1
		macs2 callpeak -q 0.1 -t $FILE -c HUH7-24KO-WT-shEV-input-1.processed.bam -f BAM -g hs --outdir ./Redundant -n $fn
		awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4}' ./Redundant/$fn\_peaks\.narrowPeak > ./Redundant/$fn\_peaks_new\.bed
		bedtools intersect -v -a ./Redundant/$fn\_peaks_new\.bed -b GRCh38_unified_blacklist.bed > ./MacsPeaks/$fn\_peaks\.bed
	done

cat *.bed > all.bed
bedtools sort -i all.bed > all_sorted.bed
bedtools merge -i all_sorted.bed > all_merged.bed
(bed->gtf by Rscript)
annotatePeaks.pl all_merged.bed hg38 > annedpeak.txt
