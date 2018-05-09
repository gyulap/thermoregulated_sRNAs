#!/usr/bin/zsh

#Degradome analysis with a modified version of CleaveLand4

PCE_Cleaveland_modified.pl -e 'All_degradome_merged.fa' -u 'miRBase_ShortStack_main_collapsed.fasta' -n 'TAIR10_cdna_simplified.fasta'

#Performing degradome analysis by degradome category separately. Only category 0, 1, and 2 are considered.

for i in {0..2}
  do PCE_CleaveLand4_modified.pl -d 'All_degradome_merged.fa_dd.txt' -g 'miRBase_ShortStack_main_collapsed.fasta_GSTAr.txt' -t -o "Degradome_T-plots_c${i}" -c $i -p 0.1 >> 'Degradome_output.txt'

     #Creating a single pdf file from the png images for every category.
     #The ImageMagick tool was used that can be downloaded from http://www.imagemagick.org

     montage -mode 'concatenate' -tile '3x4' -page 'A4' ./Degradome_T-plots_c${i}/*.png "T-plots_c${i}.pdf"

  done

#Creating a single pdf file with a header containing the legend and explaining the figure elements. The T-plots are groupped by degradome category.



pdfunite 'T-plots_header.pdf' T-plots_c*.pdf 'PCE_Figure S5.pdf'
