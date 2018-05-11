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

  for i in $reads
    do
      kallisto quant -i $index -b 100 --single -l 200 -s 20 -t $p -o "${outdir}/kallisto_files/${${i%_mRNA_processed.fastq.gz}##*/}" $i
    done
fi

#Normalizing transcript abundances with sleuth

Rscript './Scripts/PCE_sleuth.R'

#Annotating the count table

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR{a[$4]=$8"\t"$9"\t"$10"\t"$11"\t"$12; next}
     {split($1, b, "."); s = b[1];
     if ($1 == "target_id") {
       print $0"\tAnnotation1\tAnnotation2\tAnnotation3\tAnnotation4\tAnnotation5"
     }
     else if (s in a) {
       print $0, a[s]
     }
     else {
       print $0"\tNA\tNA\tNA\tNA\tNA"
     }
     }' './Auxiliary_files/TAIR10_genes_annotated8.bed' "${outdir}/kallisto_isoform_table.txt" > "${outdir}/kallisto_isoform_table_annotated.txt"
