#!/bin/bash

module load meme
# this is for motif enrichment using ame from the meme suite
# it will iterate through every .fa file within a directory (including subdirs)
# background file is the same for all and needs to be specified 


indir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/motif_enrichment/background_all_k27ac"
bckg="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/H3K27ac_distal_E7.5_union_intersect12/background/H3K27ac_distal_E7.5_union_intersect12.fa"

cd $indir
# for d in */;
for d in `find -name "*intersect12" -type d`; # only use intersect12 files which are all k27ac distal
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
