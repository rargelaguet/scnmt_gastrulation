#!/bin/bash


# this is for motif discovery (not enrichment)

#indir="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/significant_sites/fasta/test"
indir="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/significant_sites/fasta"
outdir="/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/significant_sites/motifs"
mkdir -p $outdir

cd $indir

for file in *.fa
do
out="$outdir/$file"
meme $file -oc $out -dna -revcomp -p 8 -maxsize 1000000000
done