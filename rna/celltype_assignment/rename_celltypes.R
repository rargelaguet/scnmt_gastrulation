suppressPackageStartupMessages(library(argparse))

here::i_am("rna/mapping/run/parse_sample_metadata_after_mapping.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$outfile <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

###############################################
## Rename cell types from class 1 to class 2 ##
###############################################

opts$rename.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm",
  "Parietal_endoderm" = "ExE_endoderm",
  "Visceral_endoderm" = "ExE_endoderm"
)

sample_metadata %>% .[,celltype2:=stringr::str_replace_all(celltype,opts$rename.celltypes)]


table(sample_metadata$celltype2)

###############################################
## Rename cell types from class 2 to class 3 ##
###############################################

opts$rename.celltypes <- c(
  "Caudal_epiblast" = "Epiblast",
  "Notochord" = "Endoderm",
  "Def._endoderm" = "Endoderm",
  "Gut" = "Endoderm",
  "Nascent_mesoderm" = "Mesoderm",
  "Mixed_mesoderm" = "Mesoderm",
  "Intermediate_mesoderm" = "Mesoderm",
  "Caudal_Mesoderm" = "Mesoderm",
  "Paraxial_mesoderm" = "Mesoderm",
  "Somitic_mesoderm" = "Mesoderm",
  "Pharyngeal_mesoderm" = "Mesoderm",
  "Cardiomyocytes" = "Mesoderm",
  "Allantois" = "Mesoderm",
  "ExE_mesoderm" = "Mesoderm",
  "Mesenchyme" = "Mesoderm",
  "Haematoendothelial_progenitors" = "Mesoderm",
  "Endothelium" = "Mesoderm",
  "Blood_progenitors" = "Mesoderm",
  "NMP" = "NMP",
  "Neurectoderm" = "Ectoderm",
  "Neural_crest" = "Ectoderm",
  "Forebrain_Midbrain_Hindbrain" = "Ectoderm",
  "Spinal_cord" = "Ectoderm",
  "Surface_ectoderm" = "Ectoderm",
  "Visceral_endoderm" = "Endoderm"
)

sample_metadata %>% .[,celltype3:=stringr::str_replace_all(celltype2,opts$rename.celltypes)]

table(sample_metadata$celltype3)

#################
## Save output ##
#################

fwrite(sample_metadata, args$outfile, sep="\t", na="NA", quote=F)
