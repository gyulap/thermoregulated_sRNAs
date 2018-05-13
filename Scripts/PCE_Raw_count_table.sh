#! /usr/bin/zsh

mydir=$PWD
cd './sRNA-seq/ShortStack_results'

#Extracting all the unique sequences from the ShortStack alignment file considering only the mapped reads and sorting them by total abundance.

samtools fasta -F4 'merged_alignments.bam' |\
awk '
  NR%2==0{a[$0]++}
  END{print "Sequence";
      PROCINFO["sorted_in"] = "@val_num_desc";
      for (i in a)
        {print i}
     }' > 'NR_sequences.txt'

#Counting every unique sequences per samples.

while read rg
  do
    samtools view -bu -F4 -r $rg 'merged_alignments.bam' |\
    samtools fasta - |\
    awk '
      NR==FNR{if (NR%2==0)
               {a[$0]++}
             }
      END{for (i in a)
           {print i"\t"a[i]}
         }' > "${rg%_sRNA_processed}_NR.txt"
  done < 'rg_list.txt'

#Making a new sequence count table for every sample that matches the sequence content and order in the total unique sequence table.
#In case of a missing sequence, the count is set to 0.

for i in *NR.txt
  do
    awk -v i="${i%_NR.txt}" '
      BEGIN{FS=OFS="\t"}
           NR==FNR{a[$1]=$2; next}
           {if ($1 == "Sequence")
              {print i}
            else if ($1 in a)
              {print a[$1]}
            else
              {print 0}
           }' $i 'NR_sequences.txt' > "${i%.txt}2.txt"
  done

#Pasting all the samples together.

paste 'NR_sequences.txt' Seedling*NR2.txt Root*NR2.txt Leaf*NR2.txt Flower*NR2.txt > 'Raw_count_table.txt'

rm -f *NR*

cd $mydir
