###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

# This script overlaps the output bismark files (individual CpG sites) with genomic features such as promoters, gene bodies, etc.

# - Preprocessing of annotations: collect all CpG sites
# - Preprocessing of samples: collect all CpG sites from mm10 using the package "BSgenome.Mmusculus.UCSC.mm10"
# - Annotate samples with the preprocessed annotations

## Input ##
# (1) output bismark file with one of the two following formats
# input_format=1
#   chr     pos      met_reads    nonmet_reads
#   19    3152031     1       0
#   19    3152424     0       1
# input_format=2
#   chr     pos      rate
#   19    3152031     100
#   19    3152424     0

# (2) genomic feature annotation files in BED6 format
#   1	3531624	3531843	*	CGI_1	CGI
#   1	3670619	3671074	*	CGI_2	CGI
#   1	3671654	3672156	*	CGI_3	CGI

## Output ##
# (1) a tmp folder with a punch of tsv files with the preprocessed samples and annotations. 
#   For example:
#     tmp/[SAMPLE]_[FEATURE].tsv
# (2) all.tsv: all annotations and samples in one dataframe
# output format=1
#   sample	id	anno	rate	weight
#   3289STDY6312493	super_enhancers_100	super_enhancers	42	31
#   3289STDY6312493	super_enhancers_1001	super_enhancers	0	1
#   3289STDY6312493	super_enhancers_1002	super_enhancers	0	2

# """
# To-do:
# - Currently only input_format=2 is accepted
# - Currently only output_format=1 is accepted
# - 
# """

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

# args <- list()
# args$context <- "GC"
# args$cores <- 1

io <- list()
opts <- args

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))

io$basedir <- "/Users/ricard/data/NMT-seq_ESC"
io$anno.folder <- paste0(io$basedir,"/features/filt")
io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")

# GC
if (opts$context == "GC") {
  io$in.folder <- paste0(io$basedir,"/acc/raw")
  io$out.folder <- paste0(io$basedir,"/acc/parsed")
# CG
} else if (opts$context=="CG") {
  io$in.folder <- paste0(io$basedir,"/met/raw")
  io$out.folder <- paste0(io$basedir,"/met/parsed")
}



## Options ##

# Valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

# Annotations to analyse
opts$annos <- c("LTR","LINE","genebody","prom_2000_2000_cgi","prom_2000_2000_noncgi","CGI")
# opts$annos <- "all"
if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Valid chromosomes: %s\n",paste(opts$chr_list, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_metQC==T,sample]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_accQC==T,sample]
} else{
  stop()
}
stopifnot(all(!duplicated(samples_keep)))


############################
## Preprocess annotations ##
############################

anno_list <- list()
for (i in 1:length(opts$annos)) {
  
  # Read annotation file
  anno.file <- sprintf("%s/%s.bed",io$anno.folder,opts$anno[i])
  dat_anno <- fread(anno.file ,sep="\t", header=F, select=c(1,2,3,4,5), verbose=F) %>% 
    setnames(c("chr","start","end","strand","id"))
  
  # Check that there are no weird chromosomes
  anno_list[[i]] <- dat_anno %>% .[chr%in%opts$chr_list,] %>% setkey(chr,start,end)
}
names(anno_list) <- opts$anno



#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",sample,opts$anno) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("zcat < %s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) %>%
      .[,c(1,2,3)] %>% setnames(c("chr","pos","rate")) %>%
      .[,c("start","end") := list(pos,pos)] %>% # Add 'start' and 'end' columns to do the overlap
      .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
      setkey(chr,start,end)
    
    # Overlap data with annotations
    for (anno in opts$anno) {
      # fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
      fname.out <- sprintf("%s/tmp/%s_%s",io$out.folder,sample,anno)
      if (file.exists(paste0(fname.out,".gz"))) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
        
        # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
          .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
          .[,c("sample","anno") := list(sample,anno)] %>%
          # Compute number of methylated CpGs and the corresponding methylation rates
          # NOTICE: THIS ASSUMES THAT RATE IS BINARISED
          .[,.(rate=round(mean(rate)), Nmet=sum(rate==100), N=.N), keyby=.(sample,id,anno)] %>%
          # Reorder columns
          .[,c("sample","id","anno","Nmet","N","rate")]
        
        # Store and save results
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        system(sprintf("gzip -f %s",fname.out))
        
        rm(ov)
      }
    }
    rm(dat_sample)
  }
}



# Concatenate everything and save it
for (i in opts$anno) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    # files <- list.files(paste(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample' | gzip > %s", io$out.folder, i, outfile))
    cat("\n")
  }
}
