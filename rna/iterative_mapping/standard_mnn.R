# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--test',            action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$atlas_stages <- c(
#   # "E6.5"
#   # "E6.75",
#   "E7.0"
#   # "E7.25",
#   # "E7.5",
#   # "E7.75",
#   # "E8.0",
#   # "E8.25",
#   # "E8.5"
#   # "mixed_gastrulation"
# )
# args$test <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/rna/iterative_mapping/utils.R")
  # io$atlas.marker_genes <- "/Users/ricard/data/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/Users/ricard/scnmt_gastrulation/rna/iterative_mapping/load_data.R"
} else {
  source("/homes/ricard/scnmt_gastrulation/rna/settings.R")
  source("/homes/ricard/scnmt_gastrulation/rna/iterative_mapping/utils.R")
  # io$atlas.marker_genes <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/homes/ricard/scnmt_gastrulation/rna/iterative_mapping/load_data.R"
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/rna/results/iterative_mapping")

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

source(io$script_load_data)

#############
## Run MNN ##
#############


# Joint normalisation
sce.all <- joint.normalisation(sce.query, sce.atlas, cosineNorm = TRUE)

# Select HVGs
sce.atlas <- logNormCounts(sce.atlas)
genes <- getHVGs(sce.atlas, p.value = 0.05)
# genes <- getHVGs(sce.atlas, block=as.factor(sce.atlas$sample), p.value = 0.10)

# MNN
mapping.dt <- mnn.fn(sce.all, sce.query, sce.atlas, genes = genes, npcs = 30, k = 15)

# table(mapping_dt$celltype_mapped)
# foo <- mapping.dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

##########
## Save ##
##########

fwrite(mapping.dt, sprintf("%s/standard_mnn.txt.gz",io$outdir), sep="\t")

