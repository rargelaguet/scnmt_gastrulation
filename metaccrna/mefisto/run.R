suppressPackageStartupMessages(library(MOFA2))

########################
## Create MOFA object ##
########################

# MOFAobject <- create_mofa_from_SingleCellExperiment(sce_filt.downsample, extract_metadata = TRUE)
MOFAobject <- create_mofa_from_df(data.downsampled)

# Visualise data structure
plot_data_overview(MOFAobject)

# Add covariates for MEFISTo
covariates.dt <- umap.dt %>% .[sample%in%unlist(samples_names(MOFAobject))] %>% setkey(sample) %>% .[unlist(samples_names(MOFAobject))] %>%
  melt(id.vars=c("sample"), variable.name="covariate")
MOFAobject <- set_covariates(MOFAobject, covariates = covariates.dt)

####################
## Define options ##
####################

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options
train_opts <- get_default_training_options(MOFAobject)
# train_opts$maxiter <- 1

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

mofa <- run_mofa(MOFAobject)

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

saveRDS(mofa, io$mofa.outfile)
