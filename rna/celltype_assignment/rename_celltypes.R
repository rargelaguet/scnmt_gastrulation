suppressPackageStartupMessages(library(argparse))

here::i_am("rna/mapping/run/parse_sample_metadata_after_mapping.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query samples')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
# p$add_argument('--mapping_seurat',    type="character", nargs="+", help='Results of the Seurat mapping')
p$add_argument('--mapping_mnn',    type="character",  nargs="+", help='Results of the MNN mapping')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))

## START TEST ##
# args$query_samples <- opts$samples
# args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$mapping_mnn <- file.path(io$basedir,sprintf("results/rna/mapping/mapping_mnn_%s.txt.gz",args$query_samples))
# args$mapping_seurat <- file.path(io$basedir,sprintf("results/rna/mapping/mapping_seurat_%s.txt.gz",args$query_samples))
# args$outfile <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

##########################
## Load mapping results ##
##########################

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
  "Allantois" = "ExE_mesoderm"
)

sample_metadata %>% 
	.[,celltype2:=stringr::str_replace_all(celltype,opts$rename.celltypes)]

###############################################
## Rename cell types from class 2 to class 3 ##
###############################################

# Germ layer: ectoderm, mesoderm, endoderm

# opts$rename.celltypes <- c(
#   "Erythroid1" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid3" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Rostral_neurectoderm" = "Neurectoderm",
#   "Caudal_neurectoderm" = "Neurectoderm",
#   "Anterior_Primitive_Streak" = "Primitive_Streak",
#   "Mixed_mesoderm" = "Nascent_mesoderm",
#   "Allantois" = "ExE_mesoderm"
# )

# sample_metadata %>% 
# 	.[,celltype3:=stringr::str_replace_all(celltype2,opts$rename.celltypes)]

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)
