suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))

##############
## Settings ##
##############

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$path2atlas <- "/Users/ricard/data/gastrulation10x"
  io$path2scNMT <- "/Users/ricard/data/scnmt_gastrulation"
  io$outdir <- "/Users/ricard/data/scnmt_gastrulation/rna/results/mapping"
  source("/Users/ricard/scnmt_gastrulation/rna/mapping/mapping_functions.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$path2atlas <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$path2scNMT <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation"
  io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation/rna/results/mapping"
  source("/homes/ricard/scnmt_gastrulation/rna/mapping/mapping_functions.R")
} else {
  stop("Computer not recognised")
}


opts <- list()
opts$atlas_stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)


####################
## Load 10x atlas ##
####################

# Load atlas metadata
meta_atlas <- fread(paste0(io$path2atlas,"/sample_metadata.txt.gz")) %>%
  .[stripped==F & doublet==F & stage%in%opts$atlas_stages]

# Load atlas SingleCellExperiment
sce_atlas  <- readRDS(paste0(io$path2atlas,"/processed/SingleCellExperiment.rds"))[,meta_atlas$cell]

##########################
## Load scNMT-seq query ##
##########################

# Load query metadata
meta_nmt <- fread(paste0(io$path2scNMT,"/sample_metadata.txt")) %>%
  .[pass_rnaQC==T & stage%in%c("E6.5","E7.5")]
table(meta_nmt$stage)

# Load query SingleCellExperiment
sce_nmt  <- readRDS(paste0(io$path2scNMT,"/rna/SingleCellExperiment.rds"))[,meta_nmt$id_rna]

#############
## Prepare ## 
#############

# Data structure required for the mapping...
meta_scnmt <- list()
meta_scnmt$cell <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$cells <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$stage <- meta_nmt$day[match(colnames(sce_nmt), meta_nmt$id_rna)]

# Filter out genes with little expression
nmt.genes <- names(which(rowSums(counts(sce_nmt))>25))
atlas.genes <- names(which(rowSums(counts(sce_atlas))>10))

# Subset genes that are present in both data sets
genes <- intersect(nmt.genes, atlas.genes)
sce_nmt  <- sce_nmt[genes,]
sce_atlas <- sce_atlas[genes,]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  query_sce = sce_nmt, query_meta = meta_scnmt, 
  k = 25
)

##########
## Save ##
##########

saveRDS(mapping, paste0(io$outdir,"/mapping.rds"))
fwrite(mapping$mapping.dt, paste0(io$outdir,"/mapping.txt.gz"), sep="\t")

