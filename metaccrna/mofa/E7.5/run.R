suppressMessages(library(MOFA))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(scater))
suppressMessages(library(reticulate))
suppressMessages(library(argparse))

# print(py_config())
# print(.libPaths())

# Define arguments
p <- ArgumentParser(description='')
p$add_argument('-i','--trial', type="integer", help='Trial number')
args <- p$parse_args(commandArgs(TRUE))
if (is.null(args$trial)) args$trial <- 1

# Load settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation/metaccrna/mofa/E7.5/load_settings.R")
  use_python("/Users/ricard/anaconda3/bin/python")
} else {
  source("/homes/ricard/gastrulation/metaccrna/mofa/E7.5/load_settings.R")
  use_python("/nfs/software/stegle/users/ricard/anaconda/bin/python")
}

# Load data
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation/metaccrna/mofa/E7.5/load_data.R")
} else {
  source("/homes/ricard/gastrulation/metaccrna/mofa/E7.5/load_data.R")
}


# Create MOFAobject
MOFAobject <- createMOFAobject(all_matrix_list)

# Data processing options
DataOptions <- getDefaultDataOptions()

# Model options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 10

# Training options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 5000
TrainOptions$tolerance <- 0.01
TrainOptions$DropFactorThreshold <- 0
TrainOptions$seed <- args$trial

# Prepare MOFAobject for training
MOFAmodel <- prepareMOFA(MOFAobject, 
  DataOptions = DataOptions, 
  ModelOptions = ModelOptions, 
  TrainOptions = TrainOptions
)

# Train the model
outfile <- sprintf("%s/hdf5/model_%d.hdf5",io$outdir,args$trial)
model <- runMOFA(MOFAmodel, outfile)
