#!/usr/bin/zsh

outdir='./Degradome-seq'
degreads="${outdir}/All_degradome_merged.fasta"
mirnas='./sRNA-seq/ShortStack_results/MIRNAs/miRBase_ShortStack_main_collapsed.fasta'
transcripts='./Auxiliary_files/TAIR10_noMIR_cdna_simplified.fasta'

#Preparing input files for Cleaveland

gunzip "${transcripts}.gz"
gunzip -c ${outdir}/Processed_sequences/*.fastq.gz | sed -n '2~4p' | awk '{count++; print ">deg_"count"\n"$0}' > $degreads

#Performing degradome analysis for all miRNAs. Only category 0, 1, and 2 are considered.

./Scripts/PCE_CleaveLand4_modified.pl -e $degreads -u $mirnas -n $transcripts -t -o "${outdir}/Degradome_T-plots" -c 2 -p 0.1 >> "${outdir}/Degradome_output.txt"

#Creating a single pdf file from the T-plot images for every category.
#The ImageMagick tool was used that can be downloaded from http://www.imagemagick.org

for i in {0..2}
  do
    montage -mode 'concatenate' -tile '3x4' -page 'A4' ${outdir}/Degradome_T-plots/*c${i}.png "${outdir}/Degradome_T-plots/T-plots_c${i}.pdf"
  done

#Creating a single pdf file with a header containing the legend and explaining the figure elements. The T-plots are groupped by degradome category.

pdfunite './Auxiliary_files/T-plots_header.pdf' ${outdir}/Degradome_T-plots/T-plots_c*.pdf "${outdir}/PCE_Figure S5.pdf"
