#!/bin/bash

module load meme
# this is for motif enrichment


indir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/H3K27ac_distal_E7.5_union_intersect12/up"
#background="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/H3K27ac_distal_E7.5_union_intersect12/background/H3K27ac_distal_E7.5_union_intersect12.fa"
background="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/fasta/K27ac_distal_all.fa"

outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/motif_enrichment/up/all_k27ac_back/H3K27ac_distal_E7.5_union_intersect12"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"

mkdir -p $outdir

cd $indir

for file in *.fa
do
back="$background"
out="$outdir/$file"
echo -e $file
echo -e $back
ame --oc $out --method fisher --pvalue-report-threshold 1 --scoring avg --control $back $file $pwm
done