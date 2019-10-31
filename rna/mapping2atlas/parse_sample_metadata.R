###############################################################
## Script to add lineage information to sample metadata file ##
###############################################################

library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
io$mapping.results <- "/Users/ricard/data/gastrulation/rna/mapping_10x/mapping.rds"
io$sample_metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"

####################
## Define options ##
####################

opts <- list()

# Level_1 to Level_2 lineage annotations
opts$aggregated_lineages <- c(
  # Mesoderm lineages
  "Pharyngeal_mesoderm" = "Mesoderm",
  "Paraxial_mesoderm" = "Mesoderm",
  "ExE_mesoderm" = "Mesoderm",
  "Mesenchyme" = "Mesoderm",
  "Mixed_mesoderm" = "Mesoderm",
  "Blood_progenitors_2" = "Mesoderm",
  "Haematoendothelial_progenitors" = "Mesoderm",
  "Mature_mesoderm" = "Mesoderm",
  "Intermediate_mesoderm" = "Mesoderm",
  "Somitic_mesoderm" = "Mesoderm",
  "Caudal_mesoderm" = "Mesoderm",
  "Nascent_mesoderm" = "Mesoderm",
  
  # Endoderm lineages
  "Gut" = "Endoderm",
  "Embryonic_endoderm" = "Endoderm",
  "Def._endoderm" = "Endoderm",
  "Notochord" = "Endoderm",
  
  # ExE endoderm lineages
  "Parietal_endoderm" = "Visceral_endoderm",
  "ExE_endoderm" = "Visceral_endoderm",
  
  # Primitive Streak lineages
  "Caudal_epiblast" = "Primitive_Streak",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  
  # Ectoderm lineages
  "Rostral_neurectoderm" = "Ectoderm",
  "Rostral_neuroectoderm" = "Ectoderm",
  "Surface_ectoderm" = "Ectoderm"
)

###########################
## Load mapping results  ##
###########################

# The output of the mapping algorithm is a data.frame with three important columns:
# - cell: cell ID (id_rna in sample_metadata.txt)
# - celltype.mapped: lineage label
# - stage.mapped: stage label
mapping <- readRDS(io$mapping.results)$mapping %>% as.data.table %>%
  .[,c("cell","celltype.mapped","stage.mapped")] %>%
  setnames("cell","id_rna")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample_metadata)

###############################################
## Merge mapping results and sample metadata ##
###############################################

# NOTE: THERE IS NO NEED TO RUN THIS AGAIN

# The mapping was done only for E6.5 and E7.5 cells, so we need to include
# the E4.5 cells that were not used for the mapping. Hence the all.y=TRUE
# foo <- merge(mapping,sample_metadata, by="id_rna", all.y=T)  %>%
#   .[is.na(celltype.mapped),c("celltype.mapped","stage.mapped"):=list(lineage10x_2,stage)] %>%
#   (...)

################################################################
## Create level 2 annotations by aggregating similar lineages ##
################################################################

sample_metadata %>%
  .[,lineage10x_3:=lineage10x] %>%
  .[,lineage10x_3:=stringr::str_replace_all(lineage10x,opts$aggregated_lineages)]

################################
## Some additional processing ##
################################

sample_metadata %>%
  # PGC are not expected
  .[lineage10x%in%c("PGC"), lineage10x_2:=NA] %>%
  # At E6.5 we do not expect any ectoderm, they are just epiblast cells
  .[stage=="E6.5" & lineage10x%in%c("Rostral_neurectoderm","Surface_ectoderm"), lineage10x_2:="Epiblast"]

##########
## Save ##
##########

fwrite(sample_metadata, io$sample_metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)
