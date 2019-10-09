#!/bin/bash

module load meme
# this is for motif enrichment using ame from the meme suite
# it will iterate through every .fa file within a directory (including subdirs)
# background file is the same for all and needs to be specified 


indir="/bi/scratch/Stephen_Clark/gastrulation_data/H3K27ac"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/H3K27ac/motif_enrichment"
bckg="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/H3K27ac_distal_E7.5_Ect_intersect12/background/H3K27ac_distal_E7.5_Ect_intersect12.fa"

cd $indir
# for d in */;
for d in `find -name "*fasta" -type d`; 
  do 
  echo $d
  fa=${d/%\/}
  #cd "$indir/$d"
  for file in `find "$indir/$d" -name "*.fa" -type f`;
    do
    echo "input: $file"
    # skip background files....
    if [[ $file == *"background"* ]];
      then
      echo "file is background - skipping"
      continue
    fi
    
    outfile=${file/$indir/$outdir}
    mkdir -p $outfile
    # skip if output file exists already...
    if [ -f "$outfile/ame.txt" ];
      then
      echo "ame file already exists, skipping to next"
      continue
    fi
    echo "output: $outfile"
    ame --oc $outfile --method fisher --pvalue-report-threshold 1 --scoring avg --control $bckg $file $pwm
    done
done
