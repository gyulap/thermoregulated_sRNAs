#!/usr/bin/zsh

p=$(egrep -c '^processor' '/proc/cpuinfo')
m=$(egrep 'MemTotal' '/proc/meminfo' | awk '{if ($1 > 1048576) {printf "%iK\n", $2/8} else {print "768M"} }')
url='https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas'
genomefile="./Auxiliary_files/$(basename $url)ta"
outdir='./sRNA-seq/ShortStack_results'
reads=(./sRNA-seq/Processed_sequences/*processed.fastq.gz)
mirnas='./sRNA-seq/ShortStack_results/MIRNAs/miRBase_ShortStack_main_collapsed.fasta'
degreads='./Degradome-seq/All_degradome_merged.fa'

#Downloading the TAIR10 genome file from the TAIR site

if [[ ! -f $genomefile ]]; then
  wget $url -O - | awk '{if ($1 ~ /^>/) {print ">chr"toupper(substr($1, 2, 1))} else {print $0}}' > $genomefile
fi

#Prediction and characterization of sRNA loci with ShortStack

if [[ ! -f "${outdir}/merged_alignments.bam" ]]; then
  ShortStack --genomefile $genomefile --readfile $reads --outdir $outdir --bowtie_cores $p --sort_mem $m --bowtie_m 1000 --mincov 5 --ranmax 'none'
  awk 'BEGIN{FS=OFS="\t"}{if (NR == 1 || $2 ~ /Cluster/) {print $0}}' 'Counts.txt' > 't' && mv -f 't' 'Counts.txt'
  samtools view -H "${outdir}/merged_alignments.bam" | awk -F "\t" '/^@RG/{print substr($2, 4, length($2))}' > "${outdir}/rg_list.txt"
fi

#Creating the mapping statistics (Table S1A)

if [[ ! -f "${outdir}/Table_S1A.txt" ]]; then
  ./Scripts/PCE_sRNA_mapping_statistics.sh
fi

#Extract the main miRNA and miRNA* sequences from ShortStack-predicted miRNA loci
#and merging them with the miRBase miRNAs

if [[ ! -f "${outdir}/MIRNAs/Main_miRNAs_collapsed.fasta" ]]; then
  ./Scripts/PCE_extract_main_miRNAs_from_ShortStack_MIRNA_files.sh
fi

#Creating a sequence count table from the ShortStack alignment file

if [[ ! -f "${outdir}/Raw_count_table.txt" ]]; then
  ./Scripts/PCE_Raw_count_table.sh
fi

#Getting the miRNA counts



#Creating heatmap and expression table of the thermoregulated miRNAs


#Differential expression of the sRNA loci with DESeq2

if [[ ! -d "${outdir}/DESeq2" ]]; then
  Rscript './Scripts/PCE_DESeq2.R'
  Rscript './Scripts/PCE_MA-plot.R'
fi

#Creating the expression table of all the thermoregulated sRNA loci



#Creating heatmaps and expression tables of the 21-nt phasiRNA-producing loci


#Determining the phase initiating miRNA for the phasiRNA-producing loci

#PhaseTank --genome '.Auxiliary_files/TAIR10_chr_all.fa' --lib '21nt_NR.fasta' --miR $mirnas --degradome $degreads --trigger --size 21 --dir './sRNA-seq/PhaseTank_results'


#Creating genome browser track for the 24-nt sRNAs from the ShortStack alignment file

if [[ ! -d "${outdir}/Genome_browser_tracks" ]]; then
  ./Scripts/PCE_bedgraph.sh
fi

#Creating heatmaps and expression tables of the 24-nt siRNA loci


