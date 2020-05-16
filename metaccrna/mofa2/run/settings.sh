job() {
	command="bsub -M $(( $1 * 1024 )) -n $2 \
	-q $3 ${@:4}"
	echo $command
	eval $command
}

gpujob() {
	command="bsub -M $(( $1 * 1024 )) -n $2 \
	-q $3 -P gpu -gpu - ${@:4}"
	echo $command
	eval $command
}

#########
## I/O ##
#########

if echo $HOSTNAME | grep -q "ricard"; then
	input_file="/Users/ricard/data/gastrulation/mofa2/data.txt.gz"
	output_folder="/Users/ricard/data/gastrulation/mofa2/hdf5/test"
elif echo $HOSTNAME | grep -q "ebi"; then
	input_file="/hps/nobackup2/research/stegle/users/ricard/gastrulation/mofa/data.txt.gz"
	output_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation/mofa/hdf5_2"
else
	echo "Computer not recognised"; exit
fi

#############
## Options ##
#############

# Number of trials
ntrials=1

# Number of factors
factors=15

# Convergence settings
convergence_mode="medium"
maxiter=500

# ELBO settings
elbofreq=5
start_elbo=150



## TEST ##
factors=10
maxiter=10
elbofreq=100
## TEST ##

