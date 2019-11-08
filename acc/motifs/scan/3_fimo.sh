#!/bin/bash

module load meme
# this is for motif scanning using fimo
# iterates over every .fa file within the specified directory

indir="/bi/scratch/Stephen_Clark/gastrulation_data/features/fasta"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/motif_scan/"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/motif_scan/pwm_subset.txt"

mkdir -p $outdir

cd $indir

for file in *intersect12.fa
do
  out="$outdir/$file"
  echo -e $file
#  if [ -f "$out/fimo.txt" ];
  if [ -f "$out/cisml.xml" ];
      then
      echo "fimo file already exists, skipping to next"
      continue
    fi
  fimo --oc $out --no-qvalue $pwm $file
done