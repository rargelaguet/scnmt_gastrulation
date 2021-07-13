job() {
	command="bsub -M $(( $1 * 1024 )) -n $2 -q research-rh74 ${@:3}"
	echo $command
	eval $command
}

# # Stochastic settings
number_samples_to_mask=( 100 150 200 250 )
seed=( 1 2 3 )

# number_samples_to_mask=( 25 50 )
# seed=( 1 )
 
for i in "${number_samples_to_mask[@]}"; do
	for j in "${seed[@]}"; do
		cmd="Rscript run.R --number_samples_to_mask $i --seed $j"
		job 10 1 $cmd
	done
done

