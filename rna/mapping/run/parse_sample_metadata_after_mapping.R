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

# MNN
mapping_mnn.dt <- args$mapping_mnn %>% map(~ fread(.)) %>% rbindlist %>% setnames("cell","id_rna")
stopifnot(mapping_mnn.dt$id_rna%in%sample_metadata$id_rna)

# Seurat
# mapping_seurat.dt <- args$mapping_seurat %>% map(~ fread(.)) %>% rbindlist %>% setnames("cell","id_rna")
# stopifnot(mapping_seurat.dt$id_rna%in%sample_metadata$id_rna)

###########
## Merge ##
###########

# mapping.dt <- merge(mapping_mnn.dt, mapping_seurat.dt, by="id_rna", suffixes=c("_mnn","_seurat"))
# to.save <- sample_metadata %>% merge(mapping.dt, by="id_rna", all.x=TRUE)

to.save <- sample_metadata %>% merge(mapping_mnn.dt, by="id_rna", all.x=TRUE)

#######################
## Rename cell types ##
#######################

opts$rename.celltypes <- c(
  "Erythroid1" = "early_Erythroid",
  "Erythroid2" = "early_Erythroid",
  "Erythroid3" = "late_Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm"
)

to.save %>% 
	# .[,celltype.mapped_seurat:=stringr::str_replace_all(celltype.mapped_seurat,opts$rename.celltypes)] %>%
	# .[,celltype.mapped_mnn:=stringr::str_replace_all(celltype.mapped_mnn,opts$rename.celltypes)]
	.[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename.celltypes)]

to.save %>% .[,celltype_class:=sprintf("%s_%s",celltype.mapped,class)]

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)

######################
## Compare mappings ##
######################

# foo <- merge(
#   mapping_mnn.dt[,c("id_rna","celltype.mapped")] %>% setnames("celltype.mapped","celltype_mnn"),
#   mapping_seurat.dt[,c("id_rna","predicted.id")] %>% setnames("predicted.id","celltype_seurat"),
#   by = c("id_rna")
# )
