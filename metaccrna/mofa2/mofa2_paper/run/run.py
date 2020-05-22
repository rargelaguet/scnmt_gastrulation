from mofapy2.run.entry_point import entry_point
import pandas as pd
import argparse

######################
## Define arguments ##
######################

p = argparse.ArgumentParser( description='' )

# I/O options
p.add_argument( '--input_file',               type=str,                required=True,           help='Input data file (matrix format)' )
p.add_argument( '--outfile',               type=str,                required=True,           help='Output file to store the model (.hdf5)' )

# Model options
p.add_argument( '--factors',            type=int,              default=25,             help='Number of factors' )

# Training options
p.add_argument( '--seed',                  type=int,                default=0,               help='Random seed' )
p.add_argument( '--verbose',  action="store_true",                              help='Do stochastic inference?' )
p.add_argument( '--start_elbo',            type=int,              default=1,             help='Start of ELBO computation' )
p.add_argument( '--elbo_freq',            type=int,              default=1,             help='Frequency of ELBO computation' )
p.add_argument( '--iterations',       type=int,              default=100,              help='Number of iterations')
p.add_argument( '--convergence_mode',       type=str,              default="medium",              help='Convergence mode')

# GPU acceleration
p.add_argument( '--gpu_mode',  action="store_true",                              help='GPU mode?' )

# Stochastic inference options
# p.add_argument( '--stochastic_inference',  action="store_true",                              help='Do stochastic inference?' )
# p.add_argument( '--batch_size',            type=float,              default=0.5,             help='Batch size (fraction of samples)' )
# p.add_argument( '--learning_rate',            type=float,              default=0.75,             help='Batch size (fraction of samples)' )
# p.add_argument( '--forgetting_rate',       type=float,              default=0.,              help='Forgetting rate for stochastic inference')

args = p.parse_args()

###############
## Load data ##
###############

data = pd.read_csv(args.input_file, delimiter="\t", header=0)

#####################
## Train the model ##
#####################

# initialise entry point    
ent = entry_point()

# Set data
ent.set_data_df(data)
# ent.set_data_matrix(data, views_names=views, groups_names=groups)

# Set model options
ent.set_model_options(factors=args.factors, spikeslab_weights=True, ard_factors=True, ard_weights=True)

# Set training options
if args.seed == 0: args.seed = None
ent.set_train_options(iter = args.iterations, convergence_mode = args.convergence_mode, 
	startELBO = args.start_elbo, freqELBO = args.elbo_freq, 
	gpu_mode = args.gpu_mode, verbose = args.verbose, seed = args.seed, dropR2 = None)

# Build the model
ent.build()

# Train the model
ent.run()

# Save the model
ent.save(args.outfile)