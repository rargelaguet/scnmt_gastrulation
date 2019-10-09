#!/bin/bash

module load meme
# this is for motif enrichment using ame from the meme suite
# it will iterate through every .fa file within a directory (including subdirs)
# background file is the same for all and needs to be specified 


file="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/test/fasta/H3K27ac_distal_E7.5_Ect_intersect12/up//E7.5EpiblastEctoderm_vs_E7.5MesodermEndoderm_H3K27ac_distal_E7.5_Ect_intersect12.fa"
pwm="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level/test/motif_enrichment/background_all_k27ac"
bckg="/bi/scratch/Stephen_Clark/gastrulation_data/acc/differential/feature_level//fasta/H3K27ac_distal_E7.5_union_intersect12.fa"


ame --oc $outfile --method fisher --pvalue-report-threshold 1 --scoring avg --control $bckg $file $pwm


