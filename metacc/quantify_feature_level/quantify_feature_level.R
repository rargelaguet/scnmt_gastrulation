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

suppressMessages(library(stringr))
suppressMessages(library(argparse))

#####################
## Parse arguments ##
#####################

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('--context',  type="character",              required=TRUE,  help='cg/CG or gc/GC')
p$add_argument('--annos',    type="character",  nargs="+",  required=TRUE,  help='Genomic annotation')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$context <- "CG"
# args$anno <- c("prom_200_200")
## END TEST ##

# Define what context to look at: CG (MET) or GC (ACC)
args$context <- toupper(args$context)
stopifnot(args$context %in% c("CG","GC"))

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

#########
## I/O ##
#########

if (args$context == "GC") {
  io$in.folder <- io$acc_data_raw
  io$out.folder <- io$acc_data_parsed
} else if (args$context=="CG") {
  io$in.folder <- io$met_data_raw
  io$out.folder <- io$met_data_parsed
}

#############
## Options ##
#############

# Annotations to analyse
# args$annos <- "all"
# args$annos <- c(
#   "BG207_BG295_Tal1_Mesoderm_intersected_with_atac",
#   "BG251_SLX7049_Tal1_intersected_with_atac",
#   "D340004_Scl_intersected_with_atac"
# )

# if (is.null(args$anno)) {
#   args$anno <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
# }

# Define cells
if (args$context=="CG") {
  samples_keep <- fread(io$metadata) %>% .[!is.na(id_met),id_met]
} else if (args$context=="GC") {
  samples_keep <- fread(io$metadata) %>% .[!is.na(id_acc),id_acc]
}

# ###########
# ## Print ##
# ###########

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$features.dir))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(args$annos, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",args$context))

######################
## Load annotations ##
######################

anno_list <- list()
for (i in 1:length(args$anno)) {
  anno_list[[i]] <- fread(sprintf("%s/%s.bed.gz",io$features.dir,args$anno[i]), header=F, select=c(1,2,3,4,5)) %>%
    setnames(c("chr","start","end","strand","id")) %>%
    .[,chr:=as.factor(chr)] %>%
    setkey(chr,start,end)
}
names(anno_list) <- args$anno


#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",samples_keep[i],args$annos) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",samples_keep[i])) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",samples_keep[i]))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,samples_keep[i]), sep="\t", verbose=F, showProgress=F) %>%
      .[,c("chr","pos","rate")] %>%
      .[,c("start","end") := list(pos,pos)] %>% # Add 'start' and 'end' columns to do the overlap
      .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
      setkey(chr,start,end)
    
    # Sanity check
    stopifnot(all(dat_sample$rate %in% c(0,1)))
    
    # Overlap data with annotations
    for (anno in args$annos) {
      fname.out <- sprintf("%s/tmp/%s_%s.gz",io$out.folder,samples_keep[i],anno)
      if (file.exists(fname.out)) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",samples_keep[i],anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",samples_keep[i],anno))
        
        # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
          .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
          .[,c("sample") := list(samples_keep[i])] %>%
          # Compute number of methylated CpGs and the corresponding methylation rates
          .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), N=.N), keyby=.(samples_keep[i],id,anno)] %>%
          # Reorder columns
          .[,c("sample","id","Nmet","N","rate")]
        
        # Store and save results
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE)
      }
    }
  }
}

#####################################
## Concatenate everything and save ##
#####################################

for (i in args$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample|id_met|id_acc' | pigz > %s", io$out.folder, i, outfile))
  }
}
