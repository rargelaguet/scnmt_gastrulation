suppressPackageStartupMessages(library(MOFA2))

#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/imputation/load_settings.R")
io$outdir <- paste0(io$basedir,"/metaccrna/mofa2")

################
## Load model ##
################

MOFAobject <- readRDS(paste0(io$outdir,"/model.rds"))

# plot_data_overview(MOFAobject)

samples_metadata(MOFAobject)$technology <- ifelse(is.na(samples_metadata(MOFAobject)$pass_metQC),"scRNA-seq","scNMT-seq")
# table(samples_metadata(MOFAobject)$technology )

################
## Imputation ##
################

MOFAobject <- impute(MOFAobject)

######################
## Data exploration ##
######################

plot_variance_explained(MOFAobject, x="view", y="factor", max_r2 = 10) +
  theme(
    axis.text.x = element_text(colour="black",size=rel(1), angle=25, hjust=1),
  )

plot_factor(MOFAobject, factors=5, dot_size = 3, group_by = "stage")

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

MOFAobject <- run_umap(MOFAobject, factors="all")
MOFAobject <- run_tsne(MOFAobject, factors="all")
plot_dimred(MOFAobject, method="TSNE", color_by = "stage")
plot_dimred(MOFAobject, method="TSNE", color_by = "stage_lineage")


###############
## OLD STUFF ##
###############

# met.views <- grep("met",views_names(MOFAobject), value=T)
# acc.views <- grep("acc",views_names(MOFAobject), value=T)

# Predictions
# pred.dt <- lapply(acc.views, function(m) 
#   as.data.frame(predict(MOFAobject, views=m)[[1]][[1]]) %>% 
#     as.data.table(keep.rownames = T) %>% setnames("rn","feature") %>% 
#     melt(id.vars="feature", variable.name="sample", variable.factor=F) %>% .[,view:=m]
# ) %>% rbindlist %>% .[,c("view","feature","sample","value")] %>%
#   .[,value:=100*2**value/(1+2**value)]
# hist(pred.dt$value)
# fwrite(pred.dt, paste0(io$outdir,"/imputed_acc_data.txt.gz"))


# Imputation
# imputed.dt <- get_imputed_data(MOFAobject, views=met.views, as.data.frame = T) %>% 
#   as.data.table %>% .[,group:=NULL] %>%
#   .[,value:=100*2**value/(1+2**value)]
# hist(imputed.dt$value)
# fwrite(imputed.dt, paste0(io$outdir,"/imputed_data.txt.gz"))