#!/bin/bash

# CpG methylation 
# indir="/hps/nobackup/stegle/users/ricard/gastrulation_NMT/metacc/original/lane5913"
# outdir="/hps/nobackup/stegle/users/ricard/gastrulation_NMT/metacc/raw/test"


indir="/hps/nobackup/stegle/users/ricard/gastrulation/met/original/E5.5_scBS/2683_2684/cov"
outdir="/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/E5.5_scBS/2683_2684"
cd $indir

## GpC methylation

# CpG_report.txt.gz
# for file in *CpG_report.txt.gz
# do
# 	stage=$( echo $file | cut -d "_" -f 3,4 ); stage=${stage/'_'/'.'}
# 	plate=$( echo $file | cut -d "_" -f 5 )
# 	sample=$( echo $file | cut -d "_" -f 8 )
# 	R=$( echo $file | cut -d "_" -f 10 )
# 	echo "Processing ${plate}_${sample}_${R}..."
# 	zmore ${indir}/${file} | awk '{ if ($4>0 || $5>0) {print $1,$2,$4,$5}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${plate}_${sample}_${R}_CpG.tsv.gz
# done

# CpG.cov.gz
# for file in *CpG.cov.gz
# do
# 	stage=$( echo $file | cut -d "_" -f 3,4 ); stage=${stage/'_'/'.'}
# 	plate=$( echo $file | cut -d "_" -f 5 )
# 	sample=$( echo $file | cut -d "_" -f 8 )
# 	R=$( echo $file | cut -d "_" -f 10 )
# 	echo "Processing ${stage}_${plate}_${sample}_${R} for CpG..."
# 	zmore ${indir}/${file} | awk '{ if ($5>0 || $6>0) {print $1,$2,$5,$6}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${stage}_${plate}_${sample}_${R}_CpG.tsv.gz
# done


for file in *cov.gz
do
	stage=$( echo $file | cut -d "_" -f 3,4 ); stage=${stage/'_'/'.'}
	plate=$( echo $file | cut -d "_" -f 5 )
	sample=$( echo $file | cut -d "_" -f 8 )
	R=$( echo $file | cut -d "_" -f 10 )
	echo "Processing ${stage}_${plate}_${sample}_${R} for CpG..."
	zmore ${indir}/${file} | awk '{ if ($5>0 || $6>0) {print $1,$2,$5,$6}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${stage}_${plate}_${sample}_${R}_CpG.tsv.gz
done






## GpC accessibility

# GpC.cov.gz
for file in *GpC.cov.gz
do
	stage=$( echo $file | cut -d "_" -f 3,4 ); stage=${stage/'_'/'.'}
	plate=$( echo $file | cut -d "_" -f 5 )
	sample=$( echo $file | cut -d "_" -f 8 )
	R=$( echo $file | cut -d "_" -f 10 )
	echo "Processing ${stage}_${plate}_${sample}_${R} for GpC..."
	zmore ${indir}/${file} | awk '{ if ($5>0 || $6>0) {print $1,$2,$5,$6}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${stage}_${plate}_${sample}_${R}_GpC.tsv.gz
done

# GpC_report.txt.gz
# for file in *GpC_report.txt.gz
# do
# 	stage=$( echo $file | cut -d "_" -f 3,4 ); stage=${stage/'_'/'.'}
# 	plate=$( echo $file | cut -d "_" -f 5 )
# 	sample=$( echo $file | cut -d "_" -f 8 )
# 	R=$( echo $file | cut -d "_" -f 10 )
# 	echo "Processing ${plate}_${sample}_${R}..."
# 	zmore ${indir}/${file} | awk '{ if ($4>0 || $5>0) {print $1,$2,$4,$5}}' | sort -k 1,1 -k2,2n | gzip > ${outdir}/${stage}_${plate}_${sample}_${R}.GpC.tsv.gz
# done