suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

here::i_am("rna/processing/3_seurat_to_SCE.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
p$add_argument('--normalise',       action="store_true",                 help='Log-Normalise?')
p$add_argument('--seurat',         type="character", help='Seurat object (input)')
p$add_argument('--metadata',         type="character", help='Metadata file')
p$add_argument('--outfile',         type="character", help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args <- list()
args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
args$seurat <- file.path(io$basedir,"processed/rna/seurat.rds")
args$normalise <- FALSE
args$outfile <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
args$test <- FALSE
## END TEST ##

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% 
	.[pass_rnaQC==TRUE] %>% .[,cell:=NULL] %>% setnames("id_rna","cell")

# Load seurat
seurat <- readRDS(args$seurat)[,sample_metadata$cell]

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- as.SingleCellExperiment(seurat)

# remove logcounts assays
sce@assays@data[["logcounts"]] <- NULL

# Add metadata
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

##########################
## Compute size factors ##
##########################

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# clusts <- as.numeric(quickCluster(sce))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

###################
## Log Normalise ##
###################

if (args$normalise) {
	sce <- logNormCounts(sce)
}

##########
## Save ##
##########

saveRDS(sce, args$outfile)
