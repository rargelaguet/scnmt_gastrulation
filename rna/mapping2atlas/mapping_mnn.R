library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/gastrulation/rna/mapping2atlas/Mapping2gastrulationAtlas.R")

path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
path2scNMT <- "/Users/ricard/data/gastrulation"

####################
## Load 10x atlas ##
####################

sce_atlas  <- readRDS(paste0(path2atlas, "/SingleCellExperimentAtlas.rds"))
meta_atlas <- readRDS(paste0(path2atlas, "/sample_metadata_atlas.rds"))


####################
## Load scNMT-seq ##
####################

sce_nmt  <- readRDS(paste0(path2scNMT, "/rna/parsed/SingleCellExperiment.rds"))
meta_nmt <- read.table(file = paste0(path2scNMT, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)

# Filter
meta_nmt <- meta_nmt[meta_nmt$pass_rnaQC==T & meta_nmt$stage%in%c("E6.5","E7.5") ,]
sce_nmt <- sce_nmt[,meta_nmt$id_rna] 

meta_scnmt <- list()
meta_scnmt$cell <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$cells <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$stage <- meta_nmt$stage[match(colnames(sce_nmt), meta_nmt$id_rna)]


#############
## Prepare ## 
#############

sce_nmt  <- sce_nmt[!is.na(match(rownames(sce_nmt), rownames(sce_atlas))), ]
sce_atlas <- sce_atlas[match(rownames(sce_nmt), rownames(sce_atlas)), ]


#########
## Map ##
#########

mapping <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas, map2stage_x = c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75"),
  map_sce = sce_nmt, map_meta = meta_scnmt, mapstage_x = c("E6.5","E7.5"),
  k = 30
)


##########
## Save ##
##########

saveRDS(mapping, "/Users/ricard/data/gastrulation/mapping_10x/out/mapping.rds")
