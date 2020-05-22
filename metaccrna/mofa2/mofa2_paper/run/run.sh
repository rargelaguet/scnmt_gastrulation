###################
## Load settings ##
###################

source /Users/ricard/mofa2_gastrulation/run/settings.sh

#########
## Run ##
#########

for j in $(seq 1 $ntrials); do
	outfile="$output_folder/model_$j.hdf5"
	cmd="python run.py --input_file $input_file --outfile $outfile --factors $factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq --seed $j"
	# job 13 3 research-rh7 $cmd
	# gpujob 12 1 research-rh7 $cmd
	eval $cmd
done