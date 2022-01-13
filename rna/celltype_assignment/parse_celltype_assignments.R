here::i_am("rna/celltype_assignment/parse_celltype_assignments.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Metadata file to use as input')
p$add_argument('--celltype_assignments',    type="character",  nargs="+", help='Results of the cell type assignments')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
args$celltype_assignments <- c(
	file.path(io$basedir,"results/rna/celltype_assignment/E.5/celltype_assignment_E35.txt.gz"),
	file.path(io$basedir,"results/rna/celltype_assignment/E4.5/celltype_assignment_E45.txt.gz"),
	file.path(io$basedir,"results/rna/celltype_assignment/E5.5/celltype_assignment_E55.txt.gz")
)
args$outfile <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

###############
## Load data ##s
###############

celltype_assignments.dt <- args$celltype_assignments %>% map(~ fread(.)) %>% rbindlist
stopifnot(mapping_mnn.dt$id_rna%in%sample_metadata$id_rna)

# .[,c("celltype.mapped","celltype.score","closest.cell"):=as.character(NA)] 

###########
## Merge ##
###########

to.save <- sample_metadata %>% merge(mapping_mnn.dt, by="id_rna")

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)
