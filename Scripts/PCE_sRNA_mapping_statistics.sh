#!/usr/bin/zsh

globalqc=("Unprocessed" "Passed" "Discarded" "Mapped" "Unmapped" "Uniquely mapped" "Multi-mapped")
red=("Redundant" "Non-Redundant")
bamfile='./sRNA-seq/ShortStack_results/merged_alignments.bam'
rgfile='./sRNA-seq/ShortStack_results/rg_list.txt'
outputfile='Table_S1A.txt'

{
 printf "%s\n\n" 'Table S1A - Mapping statistics of the small RNA libraries'
 printf "\t%s\t" $globalqc
 printf "\n\t" ""
 for t in $globalqc
   do
     printf "%s\t" $red
   done
 ;} | cut -f1,2,4-6,8- > $outputfile

while read rg
  do
    {
    printf "%s" $(echo ${rg%_*_processed} | sed 's/_/ /g;s/\(\(^.*[0-9][0-9] \)\([12]$\)/\1°C\2/'); printf "%s\t" ""
    total_unprocessed_R=$(gunzip -c "./sRNA-seq/Raw_sequences/${rg%_processed}_raw.fastq.gz" | sed -n '2~4p' | wc -l)
    printf "%.0f\t" $total_unprocessed_R
    awk -v unproc="$total_unprocessed_R" '
    BEGIN {FS=OFS="\t"}
    {
     passedarray[$2]++
     if ($1 == "0" || $1 == "16") {mappedcount++; mappedarray[$2]++};
     if ($1 == "4") {unmappedcount++ ; unmappedarray[$2]++};
     if ($3 ~ /U/) {uniquecount++; uniquearray[$2]++};
     if ($3 ~ /P/) {multicount++; multiarray[$2]++};
    }
    END {
     printf "%.0f\t", NR;
     printf "%.0f\t", length(passedarray);
     printf "%.0f\t", unproc-NR;
     printf "%.0f\t", mappedcount;
     printf "%.0f\t", length(mappedarray);                    
     printf "%.0f\t", unmappedcount;
     printf "%.0f\t", length(unmappedarray);          
     printf "%.0f\t", uniquecount;
     printf "%.0f\t", length(uniquearray);
     printf "%.0f\t", multicount;
     printf "%.0f\t", length(multiarray);
    }' <(samtools view -r $rg $bamfile | cut -f2,10,16)
   ;} >> $outputfile
  done < $rgfile
