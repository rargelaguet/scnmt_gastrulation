suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--number_samples_to_mask',  type="integer",  default=0, help='Number of samples to mask')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$number_samples_to_mask <- 10
## END TEST

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  reticulate::use_python("/nfs/research1/stegle/users/ricard/conda-envs/R4/bin/python", required=TRUE)
  source("/homes/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else {
  stop()
}

###############
## Load data ##
###############

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_data.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_data.R")
} else {
  stop()
}

############################
## load pre-computed UMAP ##
############################

umap.dt <- fread(io$umap) %>%
  .[sample%in%unique(data$sample)]

sample_metadata <- sample_metadata %>% merge(umap.dt,by="sample")

data <- data[sample%in%unique(umap.dt$sample)]

################################################
## Mask epigenetic values from entire samples ##
################################################

# nsamples = sample_metadata[pass_metQC==T & pass_accQC==T,.N]
args$number_samples_to_mask <- 50
if (args$number_samples_to_mask>0) {
  print(sprintf("Masking epigenetic values from N=%s samples",length(samples_to_mask)))
  samples_to_mask <- sample(sample_metadata$sample, size=args$number_samples_to_mask)
  data[sample%in%samples_to_mask & view%in%c("motif_acc","motif_met"),value:=NA]
}

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_df(data)

# Visualise data structure
plot_data_overview(MOFAobject)

# Add covariates for MEFISTo
# covariates.dt <- umap.dt %>% .[sample%in%unlist(samples_names(MOFAobject))] %>% setkey(sample) %>% .[unlist(samples_names(MOFAobject))] %>%
#   melt(id.vars=c("sample"), variable.name="covariate")
# MOFAobject <- set_covariates(MOFAobject, covariates = covariates.dt)

####################
## Define options ##
####################

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

# Training options
train_opts <- get_default_training_options(MOFAobject)
# train_opts$maxiter <- 25

# MEFISTO options
mefisto_opts <- get_default_mefisto_options(MOFAobject)
# mefisto_opts$sparseGP <- TRUE
# mefisto_opts$frac_inducing <- 0.50


#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  model_options = model_opts,
  training_options = train_opts,
  mefisto_options = mefisto_opts
)

#####################
## Train the model ##
#####################

# io$hdf5.outfile <- paste0(io$outdir,"/test.hdf5")
# mofa <- run_mofa(MOFAobject, outfile = io$hdf5.outfile)
mofa <- run_mofa(MOFAobject)

stop()

###############################
## Add metadata to the model ##
###############################

samples_metadata(mofa) <- sample_metadata %>% 
  .[sample%in%unlist(samples_names(mofa))] %>%
  .[,group:="single_group"] %>%
  setkey(sample) %>% .[unlist(samples_names(mofa))]

##########
## Save ##
##########

io$mofa.outfile <- sprintf("%s/mefisto_imputation_N%s.rds",io$outdir,args$number_samples_to_mask)
saveRDS(list("samples_masked"=samples_to_mask, "model"=mofa), io$mofa.outfile)
