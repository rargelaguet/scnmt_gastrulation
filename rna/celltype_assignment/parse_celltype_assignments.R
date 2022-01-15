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
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping_all_samples.txt.gz")
# args$celltype_assignments <- c(
# 	file.path(io$basedir,"results/rna/celltype_assignment/E3.5/celltype_assignment_E35.txt.gz"),
# 	file.path(io$basedir,"results/rna/celltype_assignment/E4.5/celltype_assignment_E45.txt.gz"),
# 	file.path(io$basedir,"results/rna/celltype_assignment/E5.5/celltype_assignment_E55.txt.gz")
# )
# args$outfile <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

###############
## Load data ##
###############

celltype_assignments.dt <- args$celltype_assignments %>% map(~ fread(.)) %>% rbindlist
stopifnot(celltype_assignments.dt$id_rna%in%sample_metadata$id_rna)

###########
## Merge ##
###########

sample_metadata.A <- sample_metadata[id_rna%in%celltype_assignments.dt$id_rna] %>%
  .[,celltype:=NULL] %>% merge(celltype_assignments.dt, by="id_rna")
sample_metadata.B <- sample_metadata[!id_rna%in%celltype_assignments.dt$id_rna]

to.save <- rbind(sample_metadata.A[,colnames(sample_metadata.B),with=F],sample_metadata.B) %>% setorder(stage)

# Sanity checks
stopifnot(!is.na(to.save[pass_rnaQC==TRUE,celltype]))
stopifnot(sort(sample_metadata$id_rna)==sort(to.save$id_rna))

#################
## Save output ##
#################

fwrite(to.save, args$outfile, sep="\t", na="NA", quote=F)
