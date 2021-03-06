# Bioinformatics-CRISPR-dCAS9-imaging
This respository contains some of the code that I had used in our publication:

"Live cell imaging of low- and non-repetitive chromosome loci using CRISPR-Cas9",Peiwu Qin, Mahmut Parlak, Cem Kuscu, Jigar Bandaria, Mustafa Mir, Karol Szlachta,Ritambhara Singh, Xavier Darzacq, Ahmet Yildiz, & Mazhar Adli, *Nature Communications*, **2017**.

This repository contains the following files.

Bash scripts :
1. "BAM_workup_all.sh", "analyze_bam_files.sh" and "compare_ngg.sh". The first two scripts are used to generate all the potential sgRNAs in the human genome (hg19). These are further analyzed by the jupyter notebooks in the repository. "compare_ngg.sh" is used at the end of the analysis to confirm all that the sequences agree with criteria in the manuscript. In this scripts I am using Bowtie, SAMTOOLS and Bedtools.

2. Jupyter notebooks "sgRNA analysis 1.ipynb","sgRNA analysis 2.ipynb" and ""sgRNA analysis 3.ipynb" perform some preliminary analysis on the sequences that are generated by the first 2 bash scripts. "Counting sgRNA.ipynb" performs some initial statistics on the sequences.

3. "Remove_reverse_complement.ipynb" and "Rechecking_10kb_after_remove_rev_comp.ipynb" are used to remove sequences that are reverse complements so that they are not counted twice. The second file also creates 2 csv files that are the final summary of the entire analysis.

4. "Statistics_10kb.ipynb" is used to perform the final statistics on the sequences and make plots.
