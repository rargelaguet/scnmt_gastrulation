here::i_am("mapping/trajectories/parse_sample_metadata_after_mapping.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
# p$add_argument('--mapping_seurat',    type="character", nargs="+", help='Results of the Seurat mapping')
p$add_argument('--mapping_mnn',    type="character",  nargs="+", help='Results of the MNN mapping')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
args$metadata <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$mapping_mnn <- file.path(io$basedir,"results_new/rna/mapping/trajectories/blood/mapping_mnn.txt.gz")
args$outfile <- file.path(io$basedir,"results_new/rna/mapping/trajectories/blood/sample_metadata_after_mapping.txt.gz")
## END TEST ##

stopifnot(file.exists(args$mapping_mnn))

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[,c("cell","id_rna","id_met","id_acc","sample","class","celltype.mapped")] %>%
  setnames("celltype.mapped","global_mapping")

##########################
## Load mapping results ##
##########################

mapping_mnn.dt <- args$mapping_mnn %>% map(~ fread(.)) %>% rbindlist
stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

###########
## Merge ##
###########

to.save <- sample_metadata %>% 
  merge(mapping_mnn.dt, by=c("cell"))

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)
