#!/bin/bash

module load meme
# this is for motif enrichment using ame from the meme suite
# it will iterate through every .fa file within a directory (including subdirs)
# requires there to be a background file named "xxx/background/xxx.fa"


indir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/motif_enrichment/"

cd $indir
for d in */;
  do 
  echo $d
  fa=${d/%\/}
  bckg="$indir/$d/background/$fa.fa"
  if [ ! -f $bckg ];
    then
    echo "no background file for $fa"
    continue
  fi
  #cd "$indir/$d"
  for file in `find "$indir/$d" -name "*.fa" -type f`;
    do
    echo "input: $file"
    # should check if file is same as background and if so skip.....
    outfile=${file/$indir/$outdir}
    mkdir -p $outfile
    if [ -f "$outfile/ame.txt" ];
      then
      echo "ame file already exists, skipping to next"
      continue
    fi
    echo "output: $outfile"
    ame --oc $outfile --method fisher --pvalue-report-threshold 1 --scoring avg --control $bckg $file $pwm
    done
done
