#! /usr/bin/zsh

mydir=$PWD
cd "./sRNA-seq/ShortStack_results/MIRNAs"

for i in Cluster*.txt
  do
    awk -v name="${i%_Y.txt}" '
      BEGIN {FS=" "}
      {ind=match($1, "[ACGUacgu]")
       if (NR==3) {sl=length($1)/2}
       else if (NR>=5 && NR<=6 && match($1, "[.]")==1 && ind < sl) {gsub("[.]", "", $1); count5p[toupper($1)]=substr($NF, 3)}
       else if (NR>=5 && NR<=6 && match($1, "[.]")==1 && ind > sl) {gsub("[.]", "", $1); count3p[toupper($1)]=substr($NF, 3)}
      }
      END {
        for (a in count5p) {
         print name"-5p\t"a >> i"_5p.txt"
         }
        for (b in count3p) {
         print name"-3p\t"b >> i"_3p.txt"
         }
        }' $i
  done

rm -f 'Main_miRNAs.fasta'

for p in *_[35]p.txt
  do
    cat $p | sort -k 1,1r | awk 'BEGIN{FS=OFS="\t"}{gsub("U", "T", $2); gsub("u", "t", $2); print ">"$1"\n"$2}' >> 'Main_miRNAs.fasta'
  done

awk 'BEGIN{RS=">";FS="\n"}{a[$2]=a[$2]","$1}END{for (t in a) {gsub("^,", "", a[t]); print ">"a[t]"\n"t}}' 'Main_miRNAs.fasta' > 'Main_miRNAs_collapsed.fasta'

rm -f *_5p.txt *_3p.txt

cd $mydir
