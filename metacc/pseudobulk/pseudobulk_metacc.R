here::here("metacc/pseudobulk/pseudobulk_cpg_level.R")

suppressMessages(library(furrr))

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--outdir',  type="character",              help='Output directory')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--min_cells',     type="integer",    default=15,   help='Minimum number of cells')
p$add_argument('--ncores',     type="integer",    default=1,   help='Number of cores')
p$add_argument('--group_by',    type="character",  nargs="+",  help='Metadata column that represents the group')
p$add_argument('--test',    action="store_true",             help='Test mode? subset number of cells')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("metacc/pseudobulk/utils.R"))

## START TEST ##
# args <- list()
# args$indir <- file.path(io$basedir,"processed/met/cpg_level")
# args$outdir <- file.path(io$basedir,"processed/met/cpg_level/pseudobulk")
# args$featuresdir  <- file.path(io$basedir,"features/genomic_contexts")
# args$metadata <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
# args$context <- "CG"
# args$group_by <- "sample"
# args$min_cells <- 10
# args$ncores <- 1
# args$test <- TRUE
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

dir.create(args$outdir, showWarnings = F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[!is.na(celltype)]

if (args$context=="CG") {
  sample_metadata <- sample_metadata %>% .[pass_metQC==TRUE]
} else {
  sample_metadata <- sample_metadata %>% .[pass_accQC==TRUE]
}

stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

# Filter groups by minimum number of cells
sample_metadata <- sample_metadata[,N:=.N,by=c(args$group_by)] %>% .[N>=args$min_cells] %>% .[,N:=NULL]

table(sample_metadata[[args$group_by]])

##############################
## Load data and pseudobulk ##
##############################

# Parallel processing options
if (args$ncores>1) {
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}

# i <- unique(sample_metadata[[args$group_by]])[1]
for (i in unique(sample_metadata[[args$group_by]])) {
  outfile = sprintf("%s/%s.tsv.gz",args$outdir,i)
  if (file.exists(outfile)) {
    print(sprintf("%s already exists, skipping...",outfile))
  } else {
    
    # Define input files 
    if (args$context=="CG") {
      cells <- sample_metadata[eval(as.name(args$group_by))==i,id_met]
    } else {
      cells <- sample_metadata[eval(as.name(args$group_by))==i,id_acc]
    }
    
    if (args$test) cells <- head(cells,n=5)
    
    # cells <- head(cells,n=2)
    files <- paste0(args$indir, "/", cells, ".tsv.gz")
    
    # split into chunks for parallel processing
    if (args$ncores>1) {
      chunks <- ceiling(seq_along(files)/opts$chunk_size)
      file_list <- split(files, chunks)
    } else {
      file_list <- list(files)
    }
    
    # pseudobulk
    init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_sites=as.integer(NA), nonmet_sites=as.integer(NA))
    # data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
    data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
      purrr::reduce(merge_and_sum) %>%
      .[,rate:=round(100*met_sites/(met_sites+nonmet_sites))] %>%
      .[,.(chr,pos,met_sites,nonmet_sites,rate)]
    
    # filter
    data <- data %>% .[complete.cases(.)]
    
    # Save
    fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
  }
}

# Completion token
# file.create(file.path(args$outdir,"completed.txt"))

###########################
## Save group statistics ##
###########################

tmp <- table(sample_metadata[[args$group_by]])
to_save.dt <- data.table(group=names(tmp), N=tmp)
fwrite(to_save.dt, file=file.path(args$outdir,"stats.txt"), quote=F, col.names=T, sep="\t")
