here::here("metacc/stats/calculate_stats_per_chr.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--outfile',  type="character",              help='Output file')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$indir <- file.path(io$basedir,"processed/met/cpg_level")
# args$outfile <- file.path(io$basedir,"results/met/stats/met_stats_per_chr.txt.gz")
# args$context <- "CG"
# args$test <- TRUE
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# Define cells
opts$cells <- list.files(args$indir, pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)

if (args$test) opts$cells <- head(opts$cells,n=5)

###################
## Load metadata ##
###################

# sample_metadata <- fread(args$metadata)

# if (args$context=="CG") {
#   stopifnot(opts$cells%in%sample_metadata$id_met)
# } else {
#   stopifnot(opts$cells%in%sample_metadata$id_acc)
# }
# mean(opts$cells%in%sample_metadata$id_acc)

########################################
## Load data and calculate statistics ##
########################################

stats <- data.table(expand.grid(opts$cells,opts$chr)) %>% 
  setnames(c("cell","chr")) %>%
  .[,c("N","rate"):=as.numeric(NA)]

# i <- opts$cells[1]
for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",args$indir,i))) {
    print(i)

    # Load data
    # data <- fread(sprintf("%s/%s.tsv.gz",args$indir,i), sep="\t", verbose=F, showProgress=F, select=c(1,2,4)) %>%
    #   setnames(c("chr","pos","rate"))
    data <- fread(sprintf("%s/%s.tsv.gz",args$indir,i), sep="\t", verbose=F, showProgress=F) %>%
      .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
      .[chr%in%opts$chr]

    # Sanity checks
    stopifnot(unique(data$chr)%in%opts$chr)

    # Compute per-chr statistics
    for (j in opts$chr) {
      stats[cell==i & chr==j, c("N","rate"):=list(data[chr==j,.N], round(100*mean(data[chr==j,rate]),2))]
    }

  } else {
    print(sprintf("Sample %s not found",i))
  }
  
}

# Define column names
if (args$context=="CG") {
  stats %>% setnames(c("id_met","chr","nCG","met_rate"))
} else if (args$context=="GC") {
  stats %>% setnames(c("id_acc","chr","nGC","acc_rate"))
}

##########
## Save ##
##########

fwrite(stats, args$outfile, sep="\t", na = "NA", quote=F)
