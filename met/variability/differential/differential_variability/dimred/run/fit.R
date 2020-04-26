suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(MOFA2))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--anno',        type="character", nargs='+',   help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('--hvg',         type="integer",                help='Number of highly variable features')
p$add_argument('--outprefix',   type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# args <- list()
# args$anno <- c("prom_2000_2000")
# args$outprefix <- "/Users/ricard/data/gastrulation/met/results/variability/differential/differential_variability/dimred/test"
# args$hvg <- 100

############################
## Define I/O and options ##
############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run/load_settings.R")
} else {
  stop("Computer not recognised")
}

#########################################################
## Load precomputed differential variability estimates ##
#########################################################

diff_dt <- args$anno %>%
  map(~ fread(sprintf("%s/%s_epsilon.txt.gz",io$scmet.diff,.))
) %>% rbindlist

###########################
## Load methylation data ##
###########################

met_dt <- args$anno %>%
  map(~ fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,.), showProgress=F)) %>% 
  rbindlist %>% setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>%
  .[id_met%in%sample_metadata$id_met]

#################
## Filter data ##
#################

# Features have been filtered a priori when running scMET

# Filter features by minimum number of CpG sites
# met_dt <- met_dt[Ntotal>=opts$min.CpGs]

# Filter features by coverage
# met_dt <- met_dt[,N:=.N,by="id"] %>% .[N>=opts$min.cells] %>% .[,N:=NULL]

# Filter features by differential variability
diff_dt <- diff_dt[res_disp_diff_test=="E7.5+"] %>%
  .[res_disp_diff_tail_prob>=opts$tail_prob_threshold & abs(res_disp_LOR)>=opts$LOR_threshold] %>% 
  .[,abs_epsilon_diff:=abs(res_disp_A-res_disp_B)] %>%
  setorder(-abs_epsilon_diff)
hvgs <- diff_dt$id %>% head(args$hvg)

met_dt <- met_dt[id%in%hvgs]

################
## Parse data ##
################

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# prepare data for MOFA
met_dt <- met_dt %>% 
  .[,c("id_met","id","m")] %>%
  setnames(c("sample","feature","value"))

####################
## Fit MOFA model ##
####################

object <- create_mofa(met_dt)

data_opts <- get_default_data_options(object)
model_opts <- get_default_model_options(object)
model_opts$likelihoods[1] <- "gaussian"
model_opts$num_factors <- 5
model_opts$ard_weights <- FALSE
model_opts$spikeslab_weights <- FALSE

train_opts <- get_default_training_options(object)
train_opts$convergence_mode <- "fast"

object <- prepare_mofa(
  object = object, 
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

model <- run_mofa(object, sprintf("%s.hdf5",args$outprefix), save_data = FALSE)

##########################
## Extract MOFA factors ##
##########################

# Remove factors that explain no variance
factors <- names(which(model@cache$variance_explained$r2_per_factor[[1]][,1]>0.001))
model <- subset_factors(model, factors)

# Fetch factors
Z.mofa <- get_factors(model)[[1]]

############################################
## Do non-linear dimensionality reduction ##
############################################

# Run UMAP
# umap.out <- umap(Z.mofa)
# Z.umap <- umap.out$layout %>% as.data.table %>% .[,sample:=rownames(Z.mofa)]

# Run t-SNE
# tsne <- Rtsne::Rtsne(Z, check_duplicates=FALSE, pca=FALSE, theta=0.5, dims=2)
# Z.out <- tsne$Y %>% as.data.table %>% .[,sample:=rownames(Z)]

# Save
Z.mofa <- Z.mofa %>% as.data.table %>% .[,sample:=rownames(Z.mofa)]
# fwrite(Z.mofa, sprintf("%s_mofa.txt.gz",args$outprefix))
# fwrite(Z.umap, sprintf("%s_umap.txt.gz",args$outprefix))

################################################################
## Plot dimensionality reduction coloured by cell type labels ##
################################################################

to.plot <- Z.mofa %>% 
  .[,c("sample","Factor1","Factor2")] %>% setnames(c("sample","V1","V2")) %>%
  merge(sample_metadata[,c("id_met","lineage10x_2")] %>% setnames("id_met","sample"),by="sample")

pdf(sprintf("%s_scatterplot.pdf",args$outprefix), useDingbats=F, width=5.5, height=4)
plot_dimred(to.plot, color.by="lineage10x_2") +
  scale_fill_manual(values=opts$colors_lineages)
dev.off()

################
## Clustering ##
################

source("clustering.R")
