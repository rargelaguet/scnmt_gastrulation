###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-n','--cores', type="integer" ,help='Number of cores', default=1)

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

## I/O ##

io <- list()
## Own computer ##
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation/public_data/Dai_2016"
  io$anno.folder <- "/Users/ricard/data/gastrulation/features/genomic_contexts"
  io$sample.metadata <- str_c(io$basedir,"/sample_metadata.txt")
  io$in.folder <- str_c(io$basedir,"/met/cpg_level")
  io$out.folder <- str_c(io$basedir,"/met/feature_level")
  
## Cluster ##
} else {
  # Mouse
  # io$anno.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/mouse/features/genomic_contexts"
  # io$sample.metadata <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/mouse/sample_metadata.txt"
  # io$in.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/mouse/cpg_level_CCH"
  # io$out.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/mouse/feature_level_CCH"
  
  # Human
  # io$anno.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/human/features/genomic_contexts"
  # io$sample.metadata <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/human/sample_metadata.txt"
  # io$in.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/human/cpg_level_CCH"
  # io$out.folder <- "/hps/nobackup/stegle/users/ricard/Ecker_2017_feature_level/human/feature_level_CCH"
}


## Options ##

opts <- args

# opts <- list()
# opts$cores <- 1

# Valid chromosomes
opts$chr_list <- c(1:22,"X","Y")

# Annotations to analyse
opts$annos <- c(
  # "CGI",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Mes_intersect12"
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12_500",
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K4me3_E7.5_Ect",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Mes",
  # "genebody",
  # "prom_2000_2000"
  # "window2000_step1000"
)

if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Valid chromosomes: %s\n",str_c(opts$chr_list, collapse=" ")))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
samples_keep <- fread(io$sample.metadata, header=T) %>% .[,id_met]
# samples_keep <- fread(io$sample.metadata, header=T) %>% .[specie=="Homo_sapiens",sample]
stopifnot(all(!duplicated(samples_keep)))

############################
## Preprocess annotations ##
############################

# Run in parallel
# registerDoParallel(cores=args$cores)
# anno_list <- foreach(i=1:length(opts$anno)) %dopar% {
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

# Run in parallel (MEMORY LEAK!!!)
# registerDoParallel(cores=args$cores)
# invisible(foreach(i=1:length(samples_keep)) %dopar% {
for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  samples_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse cpg_level methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) %>%
      .[,c("chr","pos","met_reads","nonmet_reads","rate")] %>%
      .[,total_reads:=met_reads+nonmet_reads]
    
    # Add 'start' and 'end' columns to do the overlap
    dat_sample <- dat_sample[,c("start","end") := list(pos,pos)][,pos:=NULL] %>% .[,chr:=as.factor(chr)] %>% setkey(chr,start,end)
    
    # Overlap data with annotations
    for (anno in opts$anno) {
      fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
      if (file.exists(paste0(fname.out,".gz"))) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
        
        # Overlap and calculate methylation status for each region in the annotation by summarising over all CG sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
          .[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)*100)), keyby=.(sample,id,anno)]
        
        # Store and save results
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        system(sprintf("gzip -f %s",fname.out))
      }
    }
  }
}


# Concatenate everything and save it
for (i in opts$anno) {
  cat(i)
  # files <- list.files(str_c(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
  system(sprintf("cat %s/tmp/*_%s.tsv.gz | zgrep -v 'sample' | gzip > %s/%s.tsv.gz", io$out.folder, i, io$out.folder, i))
  # system(sprintf("cat %s/tmp/*_%s.tsv.gz > %s/%s.tsv.gz", io$out.folder, i, io$out.folder, i))
  cat("\n")
}
