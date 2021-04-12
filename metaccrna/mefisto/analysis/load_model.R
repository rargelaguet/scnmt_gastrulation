suppressMessages(library(RColorBrewer))
suppressMessages(library(MOFA2))

source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/load_settings.R")
mefisto <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/mefisto_model.rds")

##################
## Rename views ##
##################

# opts$rename.views <- c(
#   "met_Promoters" = "Promoter methylation",
#   "met_E7.5 enhancers" = "Enhancer methylation",
#   "acc_Promoters" = "Promoter accessibility",
#   "acc_E7.5 enhancers" = "Enhancer accessibility",
#   "RNA" = "RNA expression"
# )
# views(model) = stringr::str_replace_all(views(model), opts$rename.views)
# model <- subset_views(model, views=unname(opts$rename.views) )

#########################
## Add sample metadata ##
#########################

# cells <- as.character(unname(unlist(MOFA2::samples(model))))
# sample_metadata_filt <- sample_metadata %>% setkey(sample) %>% .[cells] %>%
#   .[,group:=stage]
# stopifnot(all(cells==sample_metadata_filt$sample))
# samples_metadata(model) <- sample_metadata_filt

####################
## Subset factors ##
####################

# opts$min.r2 <- 0.01
# r2 <- Reduce("+", model@cache$variance_explained$r2_per_factor)
# model <- subset_factors(model, which(r2[,"RNA expression"]>opts$min.r2))
# factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")
