suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(MOFA2))

############################
## Define I/O and options ##
############################

source("/Users/ricard/mofa2_gastrulation/load_settings.R")

if (file.exists(paste0(io$basedir,"/rds/model.rds"))) {
  
  model <- readRDS(paste0(io$basedir,"/rds/model.rds"))
  
} else {
  
  ################
  ## Load model ##
  ################
  
  file <- paste0(io$basedir,"/hdf5/test_1.hdf5")
  model <- load_model(file)
  
  ##################
  ## Rename views ##
  ##################
  
  opts$rename.views <- c(
    "met_Promoters" = "Promoter methylation",
    "met_E7.5 enhancers" = "Enhancer methylation",
    "acc_Promoters" = "Promoter accessibility",
    "acc_E7.5 enhancers" = "Enhancer accessibility",
    "RNA" = "RNA expression"
  )
  views(model) = stringr::str_replace_all(views(model), opts$rename.views)
  
  # Order views
  model <- subset_views(model, views=unname(opts$rename.views) )
  
  #########################
  ## Add sample metadata ##
  #########################
  
  cells <- as.character(unname(unlist(MOFA2::samples(model))))
  
  sample_metadata_filt <- sample_metadata %>% setkey(sample) %>% .[cells] %>%
    .[,group:=stage]
  
  stopifnot(all(cells==sample_metadata_filt$sample))
  samples_metadata(model) <- sample_metadata_filt
  
  ####################
  ## Subset factors ##
  ####################
  
  # model.orig <- model
  
  opts$min.r2 <- 0.01
  r2 <- Reduce("+", model@cache$variance_explained$r2_per_factor)
  model <- subset_factors(model, which(r2[,"RNA expression"]>opts$min.r2))
  factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")
  
  ##########
  ## Save ##
  ##########
  
  saveRDS(model, paste0(io$basedir,"/mofa/rds/model.rds"))
  
}
