suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(scran))

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/load_settings.R")
} else {
  stop()
}

opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast",
  # "E4.5_Primitive_endoderm",
  
  # E5.5
  "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm",
  
  # E6.5
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  # "E6.5_Visceral_endoderm",
  # "E6.5_Mesoderm",
  
  # E7.5
  # "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  # "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
  # "E7.5_Visceral_endoderm"
)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage] %>%
  # .[pass_rnaQC==T & pass_accQC==T & pass_metQC==T & stage_lineage%in%opts$stage_lineage] %>%
  # .[,c("id_rna","stage","lineage10x_2","stage_lineage")] %>%
  droplevels

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  features = NULL, 
  cells = sample_metadata$id_rna, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Feature selection
decomp <- modelGeneVar(sce)
hvgs <- rownames(decomp)[decomp$p.value<0.05 & decomp$mean > 0.01]

sce_filt <- sce[hvgs,]

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce_filt) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("id_rna") %>%
  .[colnames(sce_filt),] %>% DataFrame()

############################
## load pre-computed UMAP ##
############################


umap.dt <- fread("/Users/ricard/data/gastrulation/metaccrna/mofa/all_stages/umap_coordinates.txt") %>%
  merge(sample_metadata[,c("id_rna","sample")]) %>% .[,sample:=NULL] %>%
  .[id_rna%in%colnames(sce_filt)]

opts$cells <- intersect(umap.dt$id_rna,colnames(sce_filt))

sce_filt <- sce_filt[,opts$cells]
umap.dt <- umap.dt[id_rna%in%opts$cells] %>% setkey(id_rna) %>% .[opts$cells]
stopifnot(umap.dt$id_rna==colnames(sce_filt))

sce_filt$V1 <- umap.dt$V1
sce_filt$V2 <- umap.dt$V2

################
## Downsample ##
################

sce_filt.downsample <- sce_filt[,sample(1:ncol(sce_filt),size = 1000)]

########################
## Create MOFA object ##
########################

# sce_filt$stage_numeric <- as.numeric(as.factor(sce_filt$stage))

MOFAobject <- create_mofa_from_SingleCellExperiment(sce_filt.downsample, extract_metadata = TRUE)
# MOFAobject <- create_mofa_from_df(data)
MOFAobject

###############################
## Add metadata to the model ##
###############################

# samples_metadata(MOFAobject) <- sample_metadata %>% 
#   .[sample%in%unlist(samples_names(mofa))] %>%
#   .[,stage_numeric:=as.numeric(as.factor(stage))] %>%
#   setkey(sample) %>% .[unlist(samples_names(mofa))]

# Add covariates for MEFISTo
MOFAobject <- set_covariates(MOFAobject, covariates = c("V1","V2"))

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

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
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts,
  mefisto_options = mefisto_opts
)

#####################
## Train the model ##
#####################

mofa <- run_mofa(MOFAobject)

# Save
saveRDS(mofa, io$mofa.outfile)

###############################
## Add metadata to the model ##
###############################

samples_metadata(mofa) <- sample_metadata %>% 
  .[sample%in%unlist(samples_names(mofa))] %>%
  .[,group:="group"] %>%
  setkey(sample) %>% .[unlist(samples_names(mofa))]

#############################
## Plot variance explained ##
#############################

plot_variance_explained(mofa, x="view", y="factor")

##################
## Plot factors ##
##################

plot_factor_cor(mofa)

plot_factor(mofa, factors = 3, groups="all", color_by = "stage_lineage", group_by = "stage", add_violin = T, dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factor(mofa, factors = 3, color_by = "stage", group_by = "stage", add_violin = T, add_boxplot = T, add_dots=F, dodge=T) +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(mofa, factors = c(1,2), color_by = "stage")# +
  # scale_fill_manual(values=opts$celltype.colors)

plot_factors(mofa, factors = c(1,2), color_by = colMeans(mofa@data$RNA$single_group))


##########
## UMAP ##
##########

mofa <- run_umap(mofa, factors = "all")
# mofa <- run_umap(mofa, factors = 2:get_dimensions(mofa)[["K"]])
plot_dimred(mofa, method="UMAP", color_by = "stage", stroke=0.1)# +
  # scale_fill_manual(values=opts$celltype.colors)

##################
## Plot weights ##
##################

plot_weights(mofa, factor = 9, nfeatures = 10, scale = F)

###############
## Plot data ##
###############

plot_data_heatmap(mofa,
                  factor = 1,
                  features = 25,
                  view = 1,
                  denoise = FALSE,
                  legend = TRUE,
                  # min.value = 0, max.value = 6,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  scale="row",
                  annotation_samples = "Category",  annotation_colors = list("Category"=opts$colors), annotation_legend = F
)

plot_data_scatter(mofa, factor=1, view=1, color_by = "embryo", features = 8, dot_size = 2)

#############
## MEFISTO ##
#############

plot_smoothness(mofa)

to.plot <- mofa@samples_metadata[,c("sample","V1","V2")] %>% as.data.table %>%
  merge(get_factors(mofa)[[1]] %>% as.data.table(keep.rownames = T) %>% setnames("rn","sample"), by="sample") %>%
  melt(id.vars=c("V1","V2","sample"), variable.name="factor") %>%
  .[,value:=value/max(abs(value)),by="factor"]

ggscatter(to.plot, x="V1", y="V2", color="value") +
  facet_wrap(~factor) +
  scale_color_gradient2(low = "gray50", mid="gray90", high = "red")
