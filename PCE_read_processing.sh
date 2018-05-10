#!/usr/bin/zsh

#Creating the directory structure

for i in 'sRNA-seq' 'Degradome-seq' 'RNA-seq'
  do
    mkdir -p "${i}/Raw_sequences/FastQC_raw_sequences" "${i}/Processed_sequences/FastQC_processed_sequences"
  done

#Determining the number of cores on the computer to set the number of threads for adapter trimming.

p=$(grep -c '^processor' '/proc/cpuinfo')

while read line
  do
    Assay_type=$(echo $line | cut -f1)
    Sample_Name=$(echo $line | cut -f6)
    Run=$(echo $line | cut -f9)

    case $Assay_type in
      'miRNA-Seq')
        out='./sRNA-seq/Raw_sequences'
        m=18
        M=35
      ;;
      'OTHER')
        out='./Degradome-seq/Raw_sequences'
        m=20
        M=21
      ;;
      'RNA-Seq')
        out='./RNA-seq/Raw_sequences'
      ;;
    esac

    rawname="${out}/${Sample_Name}_raw.fastq.gz"
    procname="${out%/*}/Processed_sequences/${Sample_Name}_processed.fastq.gz"

#Downloading raw sequences from the SRA database, renaming them by the sample name and placing them into the appropriate directory.

    echo "Downloading $Run ${Sample_Name}.fastq.gz from SRA"

    fastq-dump --gzip $Run -Z > $rawname

#Trimming the Illumina TruSeq Small RNA 3' adapter (RA3) and doing some filtering steps.
#A quality check is performed before and after read processing using FastQC.

    fastqc -t $p -o "${out}/FastQC_raw_sequences" $rawname &&
    
    if [[ ($Assay_type == 'miRNA-Seq') ||  ($Assay_type == 'OTHER') ]]; then
      cutadapt -j $p -a 'TGGAATTCTCGGGTGCCAAGG' -m $m -M $M -q 20 --max-n=0 --discard-untrimmed $rawname | pigz -p $p > $procname
    elif [[ $Assay_type == 'RNA-Seq' ]]; then
      cutadapt -j $p -q 20 --max-n=0 $rawname | pigz -p $p > $procname
    fi

    fastqc -t $p -o "${out%/*}/Processed_sequences/FastQC_processed_sequences" $procname
  done < <(awk 'NR>1{print}' './Auxiliary_files/SraRunTable.txt')
