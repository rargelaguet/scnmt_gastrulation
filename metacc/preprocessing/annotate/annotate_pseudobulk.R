###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(stringr))
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
# args$context <- "CG"

io <- list()
opts <- args

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))

## Own computer ##
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$anno.folder <- paste0(io$basedir,"/features/filt")
  io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata_scNMT_pseudobulk.txt")
  
  # GC
  if (opts$context == "GC") {
    io$in.folder <- paste0(io$basedir,"/acc/raw/pseudobulk")
    io$out.folder <- paste0(io$basedir,"/acc/parsed/pseudobulk")
  # CG
  } else if (opts$context=="CG") {
    io$in.folder <- paste0(io$basedir,"/met/raw/pseudobulk")
    io$out.folder <- paste0(io$basedir,"/met/parsed/pseudobulk")
  }
  
## Cluster ##
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$anno.folder <- paste0(io$basedir,"/features/filt")
  io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata_scNMT_pseudobulk.txt")
  
  # GC
  if (opts$context == "GC") {
    io$in.folder <- paste0(io$basedir,"/acc/raw/pseudobulk")
    io$out.folder <- paste0(io$basedir,"/acc/parsed/pseudobulk")
  # CG
  } else if (opts$context=="CG") {
    io$in.folder <- paste0(io$basedir,"/met/raw/pseudobulk")
    io$out.folder <- paste0(io$basedir,"/met/parsed/pseudobulk")
  }
}


## Options ##

# Valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

# Annotations to analyse
opts$annos <- "all"
# opts$annos <- c(
#   "H3K27ac_distal_E7.5_Ect_intersect12",
#   "H3K27ac_distal_E7.5_Ect_intersect12_500",
#   "H3K27ac_distal_E7.5_End_intersect12",
#   "H3K27ac_distal_E7.5_End_intersect12_500",
#   "H3K27ac_distal_E7.5_Mes_intersect12",
#   "H3K27ac_distal_E7.5_Mes_intersect12_500",
#   "genebody",
#   "prom_2000_2000",
#   "prom_2000_2000_cgi",
#   "prom_2000_2000_noncgi",
#   "window2000_step1000"
# )

if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

###############
## Load data ##
###############

# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[,id_met]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[,id_acc]
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
  # if (sum(str_detect(samples_processed,paste(sample,"_",opts$anno))) == length(opts$anno)) {
  if (all(sprintf("%s_%s.tsv.gz",sample,opts$anno) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F, stringsAsFactors=F, header=T) %>%
      .[,c("start","end") := list(pos,pos)] %>%
      .[,pos:=NULL] %>% setkey(chr,start,end)
    
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

          # Compute number of methylated CpGs/GpCs and the corresponding methylation rates
          #   Option1: weighting the confidence of each CpG or GpC site
          #   Option2: treating each CpG/GpC site independently and ignoring uncertainity

          # Option 1
          .[,.(rate=round(mean(rate)*100), Nmet=sum(met_reads), N=sum(met_reads)+sum(nonmet_reads)), by=c("sample","id","anno")] %>%
          .[,rate:=round(100*Nmet/N)] %>%

          # Option 2
          # .[,.(rate=round(mean(rate)), Nmet=sum(rate>=0.5), N=.N), by=c("sample","id","anno")] %>%

          # Reorder columns
          setcolorder(c("sample","id","anno","Nmet","N","rate"))
        # Store and save results
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
        system(sprintf("gzip -f %s",fname.out))
      }
    }
    rm(dat_sample)
  }
}




# Concatenate everything and save it
for (i in opts$anno) {
  cat(i)
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample|id_met' | gzip > %s.tsv.gz", io$out.folder, i, outfile))
    cat("\n")
  }
}
  