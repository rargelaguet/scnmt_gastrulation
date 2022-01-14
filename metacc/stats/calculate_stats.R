here::here("metacc/stats/calculate_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--metadata',  type="character",              help='Cell metadata file')
p$add_argument('--outfile',  type="character",              help='Output file')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$indir <- file.path(io$basedir,"processed/acc/gpc_level")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping_all_samples.txt.gz")
# args$outfile <- file.path(io$basedir,"results/acc/stats/sample_metadata_after_acc_stats.txt.gz")
# args$context <- "GC"
# args$test <- TRUE
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# Define cells
opts$cells <- list.files(args$indir, pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)
if (args$test) opts$cells <- head(opts$cells,n=5)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

# if (args$context=="CG") {
#   stopifnot(opts$cells%in%sample_metadata$id_met)
# } else {
#   stopifnot(opts$cells%in%sample_metadata$id_acc)
# }
# mean(opts$cells%in%sample_metadata$id_acc)

########################################
## Load data and calculate statistics ##
########################################

stats <- data.table(cell=opts$cells) %>% 
  .[,c("N","rate"):=as.numeric(NA)]

for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",args$indir,i))) {
    print(i)

    # data <- fread(sprintf("%s/%s.tsv.gz",args$indir,i), sep="\t", verbose=F, showProgress=F, select=c(1,2,4)) %>%
    #   setnames(c("chr","pos","rate"))
    data <- fread(sprintf("%s/%s.tsv.gz",args$indir,i), sep="\t", verbose=F, showProgress=F) %>%
      .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
      .[chr%in%opts$chr]

    # Compute genome-wide statistics
    stats[cell==i, c("N","rate"):=list(nrow(data), round(100*mean(data$rate),2))]

  } else {
    print(sprintf("Sample %s not found",i))
  }
}

# Define column names
if (args$context=="CG") {
  to.save <- sample_metadata %>% 
    merge(stats %>% setnames(c("id_met","nCG","met_rate")), by="id_met", all.x=TRUE)
} else if (args$context=="GC") {
  to.save <- sample_metadata %>% 
    merge(stats %>% setnames(c("id_acc","nGC","acc_rate")), by="id_acc", all.x=TRUE)
}

##########
## Save ##
##########

fwrite(to.save, args$outfile, sep="\t", na = "NA", quote=F)
