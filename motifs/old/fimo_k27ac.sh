#!/bin/bash

module load meme
# this is for motif scanning using fimo
# iterates over every .fa file within the specified directory

indir="/bi/scratch/Stephen_Clark/gastrulation_data/features/fasta"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/features/motifs"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"

mkdir -p $outdir

cd $indir

for file in *.fa
do
out="$outdir/$file"
echo -e $file
fimo --oc $out $pwm $file
done