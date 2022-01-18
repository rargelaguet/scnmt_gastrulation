here::here("metacc/qc/qc.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata_met',  type="character",              help='Cell metadata file')
p$add_argument('--metadata_acc',  type="character",              help='Cell metadata file')
p$add_argument('--outfile',  type="character",              help='Output file')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$metadata_met <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
# args$metadata_acc <- file.path(io$basedir,"results/acc/qc/sample_metadata_after_acc_qc.txt.gz")
# args$outfile  <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

sample_metadata_met.dt <- fread(args$metadata_met)
sample_metadata_acc.dt <- fread(args$metadata_acc)

###########
## Merge ##
###########

sample_metadata.dt <- merge(sample_metadata_met.dt[,c("cell","pass_metQC","nCG","met_rate")], sample_metadata_acc.dt, by="cell") %>% 
  .[,c("cell", "sample", "id_met", "id_acc", "plate", "id_rna", "method", "embryo", "stage", "nCount_RNA", "nFeature_RNA", 
  "mit_percent_RNA", "rib_percent_RNA", "celltype", "celltype2", "celltype3", "celltype.score", "closest.cell", "nCG","met_rate", "nGC", "acc_rate", "pass_rnaQC", "pass_metQC", "pass_accQC"
  )]

print(head(sample_metadata.dt))

##########
## Save ##
##########

fwrite(sample_metadata.dt, args$outfile, sep="\t", na = "NA", quote=F)
