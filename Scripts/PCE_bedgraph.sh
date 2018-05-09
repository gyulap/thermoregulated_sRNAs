#!/usr/bin/zsh

normfile='norm_factors_merged.txt'

for bamfile in ${bamdir}/*_[12][157]_butter_mm1.bam
  do
    bamname=${${bamfile%_butter_mm1.bam}##*/}
    normfactor=$(awk -v bamname="${bamname}" 'BEGIN{FS=OFS="\t"} {if ($1 ~ bamname) printf "%.4f", $2}' $normfile;)
    temperature=${bamname#*_}

# Sets the track colors by temperature. Color codes are RGB (Red, Green, Blue).

    case $temperature in
      "15")
        color='51,153,255'
        altColor='102,153,204'
      ;;
      "21")
        color='51,255,51'
        altColor='102,204,102'
      ;;
      "27")
        color='255,153,51'
        altColor='204,153,102'
      ;;
    esac

# Creates the tracks for every readlength, positive and negative strands separately and then merges them into one track file.

    for readlength in {24..24}
      do
        if [[ ! -f "$outdir/${bamname}_${readlength}nt_norm.bedgraph" ]]
          then
            echo "${bamname}_${readlength}nt_norm.bedgraph"
            genomfile="../TAIR10_chr_all.txt"
            trackname="$(echo ${bamname} | sed 's/\(^.*\)_\([0-9]\{2\}\)/\1 \2°C/') ${readlength}nt"
            trackline="track type=bedGraph name=\"${trackname}\" visibility=full color=${color} altColor=${altColor} graphType=bar viewLimits=-2000.0:2000.0"
            plus=$(bedtools genomecov -bg -strand + -ibam <(bamtools filter -in $bamfile -length $readlength;) -g $genomfile -scale $normfactor;)
            minus=$(bedtools genomecov -bg -strand - -ibam <(bamtools filter -in $bamfile -length $readlength;) -g $genomfile -scale $normfactor | awk 'BEGIN{FS=OFS="\t"}{$4=-$4; print $0}';)
          cat <(echo $plus;) <(echo $minus;) | bedtools sort | sed "1i$trackline" > $outdir/${bamname}_${readlength}nt_norm.bedgraph
          unset plus minus
        fi
      done
  done
