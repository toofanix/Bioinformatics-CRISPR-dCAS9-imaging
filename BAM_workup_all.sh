#!/bin/bash
# Written by Jigar N. Bandaria
# Usage-Run the command below
# ./BAM_workup_all.sh PAM_sequences.txt chr_num.txt

while read line1; do
	#Change to where Bowtie is installed, and genome fasta files are kept
	/home/Desktop/Bioinformatics/bowtie-1.1.2/bowtie-build /home/Desktop/Genome_Fasta_Files/$line1.fa $line1
	echo
	echo Index built : $line1
	while read line2; do
 		echo  chromosome number : $line1
		echo  PAM sequence : $line2
		echo 
		/home/Desktop/Bioinformatics/bowtie-1.1.2/bowtie -t -S -a -v 0 -y -c $line1 $line2 > $line1.$line2.sam
		samtools sort $line1.$line2.sam -o $line1.$line2.sorted

	done < "$1"

	samtools merge  merged.$line1.bam $line1.*.sorted
	rm -rf *.ebwt *.sam $line1.*.sorted
	echo
	echo
	echo Removed SAMFILES
	echo Remove individual BAMFILES
	echo Removed index : $line1
	
	
done < "$2"


