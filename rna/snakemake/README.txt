#################
## Run locally ##
#################

snakemake --use-conda --cores 1
snakemake --forceall --use-conda --cores 1
snakemake --forceall --use-conda --cores 1 --dry-run

############################
## Run on the EBI cluster ##
############################

bsub -M 50 snakemake -j 1 --cluster "bsub -M {params.memory} -n 1 -q research-rh74"
bsub -M 50 snakemake -j 99 --cluster "bsub {params.other} -M {params.memory} -R rusage[mem={params.rusage}] -n 1 -o {params.outfile} -e {params.errfile} -q research-rh74"

#################################
## Run on the Babraham cluster ##
#################################

sbatch -n 1 --mem 5G snakemake --forceall -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem 5G"
snakemake --forceall -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem 12G"
snakemake -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem {params.memory}G"