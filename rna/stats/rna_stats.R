###########################################################
## Gastrulation scNMT-seq: script to calculate RNA stats ##
###########################################################

library(data.table)
library(purrr)
library(scater)
library(ggplot2)
library(RColorBrewer)

################
## Define I/O ##
################

io <- list()
io$rna <- "/Users/ricard/data/gastrulation/rna/SingleCellExperiment.rds"
io$metadata.file <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$outdir <- "/Users/ricard/data/gastrulation/rna/stats"

## Define options ##
opts <- list()

# Define stage and lineage
opts$stage_lineage <- "all"

# Define which cells to use
if (opts$stage_lineage[1] == "all") {
  opts$cells <- fread(io$metadata.file) %>% 
    .[pass_rnaQC==T,id_rna]
} else {
  opts$cells <- fread(io$metadata.file) %>% 
    .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
    .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage,id_rna]
}

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$metadata.file) %>% .[id_rna %in% opts$cells]

# SingleCellExperiment object
sce <- readRDS(io$rna)[,opts$cells]

##############################
## Calculate RNA statistics ##
##############################

stats <- data.table(
  id_rna = colnames(sce), 
  total_counts = sce$total_counts, 
  num_genes = sce$total_features_by_counts
)

# Save results
fwrite(stats, paste0(io$outdir,"/rna_stats.txt"), sep="\t", quote=F)
