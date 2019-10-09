#!/bin/bash

module load meme
# this is for motif scanning of the whole mouse genome

indir="/bi/scratch/Genomes/Mouse/GRCm38"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/features/fimo/whole_genome/new"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"

mkdir -p $outdir

cd $indir

for file in Mus_musculus.GRCm38.68.dna.chromosome.*
do
out="$outdir/$file"
echo -e $file
fimo --text --oc  $out $pwm $file 
done