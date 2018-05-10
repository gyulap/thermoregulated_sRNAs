#!/usr/bin/zsh

p=$(egrep -c '^processor' '/proc/cpuinfo')
m=$(egrep 'MemTotal' '/proc/meminfo' | awk '{if ($1 > 1048576) {printf "%iK\n", $2/8} else {print "768M"} }')
genomefile='./Auxiliary_files/TAIR10_chr_all.fa'
outdir='./sRNA-seq/ShortStack_results'
reads=(./sRNA-seq/Processed_sequences/*processed.fastq.gz)

#Prediction and characterization of sRNA loci with ShortStack

ShortStack --genomefile $genomefile --readfile $reads --outdir $outdir --bowtie_cores $p --sort_mem $m
samtools view -H "${outdir}/merged_alignments.bam" | awk -F "\t" '/^@RG/{print substr($2, 4, length($2))}' > "${outdir}/rg_list.txt"

#Creating the mapping statistics (Table S1A)

./Scripts/PCE_sRNA_mapping_statistics.sh

#Differential expression of the sRNA loci with DESeq2

Rscript './Scripts/PCE_DESeq2.R'
Rscript './Scripts/PCE_MA-plot.R'

#Extract the main miRNA and miRNA* sequences from ShortStack-predicted miRNA loci

./Scripts/PCE_extract_main_miRNAs_from_ShortStack_MIRNA_files.sh

#Creating genome browser track for the 24-nt sRNAs from the ShortStack alignment file

./Scripts/PCE_bedgraph.sh

#Creating the expression table of all the thermoregulated sRNA loci

#Creating heatmaps and expression tables of the 21-nt phasiRNA and 24-nt siRNA loci

#Creating a sequence count table from the ShortStack alignment file

./Scripts/PCE_Raw_count_table.sh

#Getting the miRNA counts

#Creating heatmap and expression table of the thermoregulated miRNAs
