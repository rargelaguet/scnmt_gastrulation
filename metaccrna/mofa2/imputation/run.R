suppressPackageStartupMessages(library(MOFA2))

###################
## Load settings ##
###################

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/imputation/load_settings.R")
io$outdir <- paste0(io$basedir,"/metaccrna/mofa2")

###############
## Load data ##
###############

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/imputation/prepare_data.R")

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

# Data options
data_opts <- get_default_data_options(MOFAobject)

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42


# Prepare MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train the model
MOFAobject <- run_mofa(MOFAobject)

#########################
## Downstream analysis ##
#########################

# Add sample metadata 
cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
samples_metadata(MOFAobject) <- sample_metadata %>% setkey(sample) %>% .[cells]

# Save
saveRDS(MOFAobject, paste0(io$outdir,"/model.rds"))
# MOFAobject <- readRDS(paste0(io$outdir,"/model.rds"))

