#!/usr/bin/zsh

annotfile='./Auxiliary_files/TAIR10_fornorm.bed'
genomefile='./Auxiliary_files/TAIR10_chr_all.txt'
outdir='./sRNA-seq/ShortStack_results'
bamfile="${outdir}/merged_alignments.bam"
rgfile="${outdir}/rg_list.txt"

if [[ ! -d "${outdir}/Genome_browser_tracks" ]]; then
  mkdir "${outdir}/Genome_browser_tracks"
fi

while read rg
  do
    mergedrg="${rg%_1_sRNA_processed}"
    normfactor=$(samtools view -F4 -L $annotfile $bamfile | egrep -c $mergedrg | awk '{printf "%.4f\n", 1000000/$0}')
    printf "%s\t%.4f\n" $mergedrg $normfactor >> "${outdir}/norm_factors_merged_RPM.txt"
    temperature="${mergedrg#*_}"

# Sets the track colors by temperature. Color codes are RGB (Red, Green, Blue).

    case $temperature in
      "15")
        color='51,153,255'
      ;;
      "21")
        color='51,255,51'
      ;;
      "27")
        color='255,153,51'
      ;;
    esac

# Creates the tracks for every readlength, positive and negative strands separately and then merges them into one track file.

    for readlength in {24..24}
      do
        trackname="$(echo ${mergedrg} | sed 's/\(^.*\)_\([0-9]\{2\}\)/\1 \2 °C/') ${readlength}nt"
        trackline="track type=bedGraph name=\"${trackname}\" visibility=full color=${color} graphType=bar viewLimits=-200.0:200.0"
        plus=$(bedtools genomecov -bg -strand + -ibam <(samtools view -bu -r ${mergedrg}* $bamfile | bamtools filter -length $readlength;) -g $genomefile -scale $normfactor;)
        minus=$(bedtools genomecov -bg -strand - -ibam <(samtools view -bu -r ${mergedrg}* $bamfile | bamtools filter -length $readlength;) -g $genomefile -scale $normfactor | awk 'BEGIN{FS=OFS="\t"}{$4=-$4; print $0}';)
        cat <(echo $plus;) <(echo $minus;) | bedtools sort | sed "1i$trackline" > "${outdir}/Genome_browser_tracks/${mergedrg}_${readlength}nt_norm.bedgraph"
      done
  done < <(egrep '_1_' $rgfile)
