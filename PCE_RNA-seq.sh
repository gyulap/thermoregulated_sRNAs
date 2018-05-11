#! /usr/bin/zsh

url='https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_cdna_20101214_updated'
index="./Auxiliary_files/$(basename $url)"
fasta="${index}.fasta"
reads=(./RNA-seq/Processed_sequences/*mRNA_processed.fastq.gz)
outdir='./RNA-seq/kallisto_results'

#Downloading the TAIR10 transcriptome from the TAIR site

if [[ ! -f $fasta ]]; then
  wget -q $url -O $fasta
fi

#Making the kallisto index for pseudoalignment

if [[ ! -f $index ]]; then
  kallisto index -i $index $fasta
fi

#Quantifying transcript abundances with kallisto

p=$(egrep -c '^processor' '/proc/cpuinfo')

if [[ ! -d './RNA-seq/kallisto_results/kallisto_files' ]]; then
  mkdir -p './RNA-seq/kallisto_results/kallisto_files'
fi

for i in $reads
  do
    kallisto quant -i $index -b 100 --single -l 200 -s 20 -t $p -o "${outdir}/kallisto_files/${${i%_mRNA_processed.fastq.gz}##*/}" $i
  done

#Normalizing transcript abundances with sleuth

Rscript './Scripts/PCE_sleuth.R'

#Annotating the count table

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{a[$1]=$3; next}
     {split($1, s, "[.]");
      if (FNR==1 || s[1] in a)
        {print $0, a[s[1]]}
     }' "./Auxiliary_files/TAIR10_functions_simplified.txt" "${outdir}/kallisto_isoform_table.txt" > "${outdir}/kallisto_isoform_table_annotated.txt"
