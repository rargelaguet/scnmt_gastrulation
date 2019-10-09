library(ggplot2)
library(data.table)
library(purrr)
library(VGAM)
library(MASS)
library(RColorBrewer)
library(argparse)

round_df <- function(df, digits) {
  nums <- names(which(vapply(df, is.numeric, FUN.VALUE = logical(1))))
  df[,(nums) := round(.SD,digits), .SDcols=nums]
  return(df)
}

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('-a',  '--anno',           type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('-s',  '--stage_lineage',  type="character",  nargs='+',  help='stage_lineage (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-o',  '--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

# args <- list()
# args$anno <- "CGI"
# args$stage_lineage <- c("E4.5_EPI","E4.5_PE")
# args$outfile <- "/homes/ricard/gastrulation/met/variability/betabinomial_model/out/test.txt"

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  source("/Users/ricard/gastrulation/met/variability/bb_model.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  source("/homes/ricard/gastrulation/met/variability/bb_model.R")
}
io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata_scNMT.txt")
io$in.data <- paste0(io$basedir,"/met/parsed")
io$in.features  <- paste0(io$basedir, "/features/filt")
io$outfile <- args$outfile
dir.create(dirname(io$outfile), recursive = TRUE)

####################
## Define options ##
####################

opts <- list()

# Define genomic contexts
opts$anno <- args$anno

# Define stage and lineage
opts$stage_lineage <- args$stage_lineage

# Define colors for plotting
opts$colors <- c(E4.5="#F8766D", E5.5="#7CAE00", E6.5="#00BFC4", E7.5="#C77CFF")

# Filtering options
opts$min.GpCs <- 10
opts$min.cells <- 10

# Define which cells to use
opts$cells <- fread(io$in.sample_metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>% 
  .[pass_metQC==TRUE & KO_3b=="not" & stage_lineage%in%opts$stage_lineage,id_met]


#############################
## Load methylation data ##
#############################

data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$in.data,opts$anno), sep="\t", quote="") %>%
  .[V1%in%opts$cells]
colnames(data) <- c("sample","id","anno","rate","Nmet","N")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(file=io$in.sample_metadata, header=T, sep="\t", stringsAsFactors=FALSE) %>%
  .[,c("id_met","stage","lineage")] %>%
  .[id_met%in%opts$cells] %>% setnames("id_met","sample")

sample_metadata[stage=="E6.75",stage:="E6.5"] %>% .[,stage_lineage:=paste(stage,lineage,sep="_")]

###########################################
## Merge methylation data and metadata ##
###########################################

data <- merge(sample_metadata, data, by="sample")

###############################
## Filter methylation data ##
###############################

# Filter features by number of CpGs
data <- data[N>=opts$min.GpCs]

# Filter features by number of cells (by stage)
for (i in unique(data$stage)) {
  data[stage==i,Ntotal:=sample_metadata[stage==i,.N]]
}
keep_cov_sites <- data %>% split(.$stage) %>% map(~ .[, cov:=.N, by=c("id","anno")] %>% .[cov >= opts$min.cells] %>% .[,id_anno:=paste(id,anno,sep="_")] %>% .$id_anno)
data <- data %>% .[,id_anno:=paste(id,anno,sep="_")] %>% .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,"Ntotal":=NULL]


#############################
## Fit beta binomial model ##
#############################

# bb_res <- tmp[, bb_mle(cbind(N, Nmet))[c("rho", "mu", "z_test", "chi2_test")], by = c("id","stage","anno")]
bb_res <- data[, bb_mle(cbind(N, Nmet)), by = c("id","stage","anno")] %>%
  .[is_conv==TRUE]

##########
## Save ##
##########

fwrite(round_df(bb_res,5), io$outfile, sep="\t")
