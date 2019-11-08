###########################################################################################
##  Script to quantify DNA methylation and chromatin accessibility over genomic features ##
###########################################################################################

# This script overlaps the individual CpG sites with genomic features such as promoters, gene bodies, etc.

## Input ##

# (1) For every cell, a long data.table with (at least ) columns ["chr", "pos", "rate"]
# Example:
#   chr     pos      rate
#   19    3152031     1
#   19    3152424     0

# IMPORTANT: For every CpG/GpC site, the the rate must be 0 or 1

# (2) genomic feature annotation files in BED6 format
#   chr start end strand id anno
#   1	3531624	3531843	*	CGI_1	CGI
#   1	3670619	3671074	*	CGI_2	CGI
#   1	3671654	3672156	*	CGI_3	CGI

suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

#####################
## Parse arguments ##
#####################

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

opts <- args

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
}

if (opts$context == "GC") {
  io$in.folder <- paste0(io$basedir,"/acc/gpc_level")
  io$out.folder <- paste0(io$basedir,"/acc/feature_level")
} else if (opts$context=="CG") {
  io$in.folder <- paste0(io$basedir,"/met/cpg_level")
  io$out.folder <- paste0(io$basedir,"/met/feature_level")
}
io$anno.folder <- paste0(io$basedir,"/features/genomic_contexts")
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

## Options ##

# Valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

# Annotations to analyse
opts$annos <- "all"
# opts$annos <- c(
#   "E12.5_intestine_H3K27ac_distal"
# )

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
  samples_keep <- fread(io$sample.metadata) %>% .[!is.na(id_met),id_met]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$sample.metadata) %>% .[!is.na(id_acc),id_acc]
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
  anno.file <- sprintf("%s/%s.bed",io$anno.folder,opts$annos[i])
  dat_anno <- fread(anno.file, select=c(1,2,3,4,5)) %>% 
    setnames(c("chr","start","end","strand","id"))
  
  # Check that there are no weird chromosomes
  anno_list[[i]] <- dat_anno %>% .[chr%in%opts$chr_list,] %>% setkey(chr,start,end)
}
names(anno_list) <- opts$annos


#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",sample,opts$annos) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,sample)) %>%
      .[,c("chr","pos","rate")] %>%
      .[,c("start","end") := list(pos,pos)] %>% # Add 'start' and 'end' columns to do the overlap
      .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
      setkey(chr,start,end)
    
    # Overlap data with annotations
    for (anno in opts$annos) {
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
          .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), N=.N), keyby=.(sample,id,anno)] %>%
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

#####################################
## Concatenate everything and save ##
#####################################

for (i in opts$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample|id_met|id_acc' | gzip > %s", io$out.folder, i, outfile))
    cat("\n")
  }
}
