#!/bin/bash
#Written by Jigar N. Bandaria
# Usage ./analyze_bam_files.sh chr_num.txt
# Desktop computer has low memory. Hence, running the code one chromosome at a time.
# The code reads the BAM alignment files from Bowtie. 20 nt that are on the 5' end are read. 
# Only sequences with 35 <= GC% <=75% are kept.
#fa_path - Replace this with there the genome fasta files for hg19 are stored.

clear 
while read line1; do
	echo ANALYZING THE CHROMOSOME $line1
	echo
	
	echo MAKING BAM, SLOPPING AND KEEPING UNIQUES FOR $line1
	echo
	time bedtools bamtobed -i merged.$line1.bam |bedtools slop -i stdin -g chrom_hg19.sizes -s -l 20 -r -3 | uniq  > slop_PAM.bed  
	echo
	
	echo GETTING GC DATA FOR $line1
	echo
	time bedtools nuc -fi /home/Desktop/Genome_Fasta_Files/$line1.fa -bed slop_PAM.bed > GC_quant.bed
	echo
	rm -r slop_PAM.bed
	
	export line1
	echo RUNNING PYTHON
	time python reading_file.py
	echo
	
	rm -r GC_quant.bed
	echo 
	
	echo GETTING FASTA
	time bedtools getfasta -tab -s -fi /home/Desktop/Genome_Fasta_Files/$line1.fa -bed remaining.bed -fo PAM1_$line1.fa
	echo
	
	awk '!/N/' PAM1_$line1.fa > tmp && mv tmp PAM1_$line1.fa
	rm -r remaining.bed
	echo
	echo $line1 DONE DONE DONE '!!!' 
	echo
done < "$1"

while read line1; do
	
	echo NUMBER of READS $line1 : `cat PAM1_$line1.fa | wc -l `
	
	echo
done < "$1"

