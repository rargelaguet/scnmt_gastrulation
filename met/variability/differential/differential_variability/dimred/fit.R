suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
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
# args$model <- "gaussian"
# args$anno <- c("distal_H3K27ac_cortex")
# args$outprefix <- "/Users/ricard/data/Ecker_2017/mouse/variability/dimred_vs_hvg/test/test"
# args$hvg <- 100

############################
## Define I/O and options ##
############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/met/variability/differential/dimred/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/met/variability/differential/dimred/load_settings.R")
} else {
  stop("Computer not recognised")
}

#########################################################
## Load precomputed differential variability estimates ##
#########################################################

diff_dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s_epsilon.txt.gz",io$scmet.diff,.))
) %>% rbindlist

# id,overall_mean,overall_res_disp,res_disp_A,res_disp_B,res_disp_LOR,res_disp_OR,res_disp_diff_tail_prob,res_disp_diff_test,evidence_thresh,efdr,efnr,anno

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

# Filter features by minimum number of CpG sites
# met_dt <- met_dt[Ntotal>=opts$min.CpGs]

# Filter features by coverage
# met_dt <- met_dt[,N:=.N,by="id"] %>% .[N>=opts$min.cells] %>% .[,N:=NULL]

# Filter features by differential variability
diff_dt <- diff_dt[res_disp_diff_test=="E7.5+"] %>%
  .[res_disp_diff_tail_prob>=opts$tail_prob_threshold & abs(res_disp_LOR)>=opts$LOR_threshold] %>% 
  .[,abs_epsilon_diff:=abs(res_disp_A-res_disp_B)] %>%
  setorder(-abs_epsilon_diff) %>%
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
  setnames(c("id_met","feature","value"))

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

# to.plot <- Z.umap %>% 
#   merge(sample_metadata,by="sample")

# pdf(sprintf("%s_1.pdf",args$outprefix), useDingbats=F, width=5.5, height=4)
# plot_dimred(to.plot, color.by="`Neuron type 1`")
# dev.off()

# pdf(sprintf("%s_2.pdf",args$outprefix), useDingbats=F, width=5.5, height=4)
# plot_dimred(to.plot, color.by="`Neuron type 2`")
# dev.off()

################
## Clustering ##
################

source("clustering.R")
