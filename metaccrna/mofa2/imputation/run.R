suppressPackageStartupMessages(library(MOFA2))

###################
## Load settings ##
###################

source("/Users/ricard/scnmt_gastrulation/metrna/mofa2/load_settings.R")
io$outdir <- paste0(io$basedir,"/metaccrna/mofa2")

###############
## Load data ##
###############

source("/Users/ricard/scnmt_gastrulation/metrna/mofa2/prepare_data.R")

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
# saveRDS(MOFAobject, paste0(io$outdir,"/model.rds"))
MOFAobject <- readRDS(paste0(io$outdir,"/model.rds"))

################
## Imputation ##
################

met.views <- grep("met",views_names(MOFAobject), value=T)
acc.views <- grep("acc",views_names(MOFAobject), value=T)

# Predictions
pred.dt <- lapply(acc.views, function(m) 
  as.data.frame(predict(MOFAobject, views=m)[[1]][[1]]) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","feature") %>% 
    melt(id.vars="feature", variable.name="sample", variable.factor=F) %>% .[,view:=m]
) %>% rbindlist %>% .[,c("view","feature","sample","value")] %>%
  .[,value:=100*2**value/(1+2**value)]
hist(pred.dt$value)
fwrite(pred.dt, paste0(io$outdir,"/imputed_acc_data.txt.gz"))


# Imputation
MOFAobject <- impute(MOFAobject)
# imputed.dt <- get_imputed_data(MOFAobject, views=met.views, as.data.frame = T) %>% 
#   as.data.table %>% .[,group:=NULL] %>%
#   .[,value:=100*2**value/(1+2**value)]
# hist(imputed.dt$value)
# fwrite(imputed.dt, paste0(io$outdir,"/imputed_data.txt.gz"))


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

plot_factor(MOFAobject, factors=1, dot_size = 3, group_by = "lineage10x_2")

plot_factor(MOFAobject, factors=1, dot_size = 3, color_by = "lineage10x_2", dodge=T)

plot_factors(MOFAobject, factors=c(1,2), dot_size = 3, color_by = "lineage10x_2") +
  scale_fill_manual(values=opts$colors_lineages)

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
