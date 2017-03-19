#!/bin/bash

samtools view no_rev_comp_grc37.sam -o no_rev_comp_grc37.bam

bedtools bamtobed -i no_rev_comp_grc37.bam |bedtools slop -i stdin -g chrom_hg19.sizes -s -l 3 -r 3 | uniq  > slop_PAM.bed

bedtools getfasta -tab -s -fi /home/jigar/Desktop/Bioinformatics/Genome_Fasta_Files/hg19.fa -bed slop_PAM.bed -fo slop_fasta.fa
