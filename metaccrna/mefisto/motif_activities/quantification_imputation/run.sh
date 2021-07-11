job() {
	command="bsub -M $(( $1 * 1024 )) -n $2 -q research-rh74 ${@:3}"
	echo $command
	eval $command
}

# # Stochastic settings
number_samples_to_mask=( 25 50 75 100 )
seed=( 1 2 3 )
# tau=( 0.05 0.15 0.25 0.5 0.75 1.0 )
# forgetting_rate=( 0.75 0.5 0.25 0 )
# 
for i in "${number_samples_to_mask[@]}"; do
	for j in "${seed[@]}"; do
		cmd="Rscript run.R --number_samples_to_mask $i --seed $j --test_mode"
		job 10 1 $cmd
	done
done

