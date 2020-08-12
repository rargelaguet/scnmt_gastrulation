suppressPackageStartupMessages(library(MOFA2))

###################
## Load settings ##
###################

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/test/load_settings.R")
io$outdir <- paste0(io$basedir,"/metaccrna/mofa2/test")

###############
## Load data ##
###############

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/test/prepare_data.R")

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

# Training options
train_opts <- get_default_training_options(MOFAobject)

# Prepare MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
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
# saveRDS(MOFAobject, paste0(io$outdir,"/model.rds"))
# MOFAobject <- readRDS(paste0(io$outdir,"/model.rds"))

######################
## Data exploration ##
######################

opts$rename.views <- c(
  "met_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (met)",
  "met_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (met)",
  "met_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (met)",
  "acc_H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers (acc)",
  "acc_H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers (acc)",
  "acc_H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers (acc)"
)
views_names(MOFAobject) = stringr::str_replace_all(views_names(MOFAobject), opts$rename.views)

plot_variance_explained(MOFAobject, x="view", y="factor") +
  theme(
    axis.text.x = element_text(colour="black",size=rel(1), angle=25, hjust=1),
  )

# plot_factor(MOFAobject, factors=1, dot_size = 3, color_by = "stage_lineage")

# plot_factor(MOFAobject, factors=1, dot_size = 3, color_by = "lineage10x_2", dodge=T)

plot_factors(MOFAobject, factors=c(1,3), dot_size = 3, color_by = "lineage10x_2")# +
  # scale_fill_manual(values=opts$colors_lineages)

plot_data_heatmap(MOFAobject, 
  factor = 1, 
  view = "Mesoderm enhancers (acc)", 
  imputed = F,
  denoise = T,
  annotation_samples = "stage_lineage",
  show_rownames = F, show_colnames = F,
  cluster_rows = T, cluster_cols = F
)


plot_data_scatter(MOFAobject, factor=1)
