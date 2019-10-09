#!/bin/bash

module load meme
# this is for motif scanning using fimo


indir="/bi/scratch/Stephen_Clark/gastrulation_data/features/fasta"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/test/motif_scan/"
#pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/test/motif_scan/pwm_subset.txt"

mkdir -p $outdir

cd $indir

file="H3K27ac_distal_E7.5_End_intersect12.fa"
out="$outdir/$file"
fimo --oc $out --no-qvalue $pwm $file