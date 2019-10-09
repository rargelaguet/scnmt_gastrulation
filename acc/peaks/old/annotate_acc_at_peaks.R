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
suppressMessages(library(purrr))
# suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
#suppressMessages(library(argparse))

# Initialize argument parser
# p <- ArgumentParser(description='')
# p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
# p$add_argument('-n','--cores', type="integer",help='Number of cores')
# 
# # Read arguments
# args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

args <- list()
args$context <- "GC"
args$cores <- detectCores()

io <- list()
opts <- args

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))


io$basedir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
# io$anno.folder <- str_c(io$basedir,"/features/sjc/motifs/motifdb/filt")
# io$anno.folder <- str_c(io$basedir,"/features/sjc/motifs/motifdb/filt/motif_only")
io$anno.folder <- str_c(io$basedir,"/features/sjc/peaks")
io$in.sample_metadata <- str_c(io$basedir,"/sample_metadata.txt")
  
io$in.folder <- str_c(io$basedir,"/acc/raw/")
# io$out.folder <- str_c(io$basedir,"/acc/parsed/motifs")
io$out.folder <- str_c(io$basedir,"/acc/parsed/peaks")


## Options ##

# Valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

# Annotations to analyse
opts$annos <- "all"
if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Valid chromosomes: %s\n",str_c(opts$chr_list, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_metQC==T,id_met]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_accQC==T,id_acc]
} else{
  stop()
}
stopifnot(all(!duplicated(samples_keep)))

# cat(sprintf("- Processing samples: %s\n",str_c(samples_keep, collapse=" ")))


############################
## Preprocess annotations ##
############################

# Run in parallel
# registerDoParallel(cores=args$cores)
# anno_list <- foreach(i=1:length(opts$anno)) %dopar% {
# anno_list <- list()
# for (i in 1:length(opts$annos)) {
#   
#   # Read annotation file
#   anno.file <- sprintf("%s/%s.bed",io$anno.folder,opts$anno[i])
#   dat_anno <- fread(anno.file ,sep="\t", header=F, select=c(1,2,3,4,5), verbose=F) %>% 
#     setnames(c("chr","start","end","strand","id"))
#   
#   # Check that there are no weird chromosomes
#   anno_list[[i]] <- dat_anno %>% .[chr%in%opts$chr_list,] %>% setkey(chr,start,end)
# }
# names(anno_list) <- opts$anno

anno_list <- map(opts$annos, ~paste0(io$anno.folder, "/", ., ".bed") %>%
                   fread(select=1:5, verbose=FALSE) %>%
                   .[complete.cases(.)]
) 
# rename annotations with non-standard characters
opts$annos <- gsub("\\(var.2)", "", opts$annos)
names(anno_list) <- opts$annos

#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
temp_dir <- sprintf("%s/temp",io$out.folder)
dir.create(temp_dir, recursive=T)


anno_list <- map2(anno_list, names(anno_list), ~.x[, anno:=.y]) %>%
  rbindlist() %>%
  setkey(chr, start, end)

registerDoParallel(cores=2)
invisible(foreach(i=(samples_keep)) %dopar% {
for (i in rev(samples_keep)){
  samples_processed <- list.files(temp_dir)
  if (sum(str_detect(samples_processed, str_c(i,"_",opts$anno))) == length(opts$anno)) {
        cat(sprintf("Sample %s already processed for all required annotations...\n",i)) 
      } else {
        cat(sprintf("Sample %s has not been processed, annotating...\n",i)) 
        
        file <- dir(io$in.folder, pattern=i, full=TRUE) %>% paste("zcat ", .)
        dt <- fread(file, select=c(1,2,5)) %>%
          .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>%
          setkey(chr, start, end) %>%
          foverlaps(anno_list, nomatch=0L) %>%
          .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), .N), .(id, anno)] %>%
          .[, sample := i]
        setkey(dt, anno)
        walk(dt[, unique(anno)], ~{
          tmp_file <- paste0(io$out.folder, "/temp/", i, "_", ., ".tsv")
          fwrite(dt[anno==.], tmp_file, sep="\t")          
        })
      }
  
})

# concatenate files
walk(opts$anno, ~{
  out_file <- paste0(io$out.folder, "/", ., ".tsv")
  if (file.exists(paste0(out_file, ".gz"))) return(paste0(basename(.), " already exists. Not overwritten!"))
  files <- dir(temp_dir, pattern=paste0(., ".tsv"), full=TRUE)
  if (length(files)==0) return()
  dt <- map(files, fread) %>%
    rbindlist()
  fwrite(dt, out_file, sep="\t")
  system(paste("gzip", out_file))
})




# # Run in parallel
# registerDoParallel(cores=opts$cores)
# invisible(foreach(i=1:length(samples_keep)) %dopar% {
# #for (i in 1:length(samples_keep)) {
#   sample=samples_keep[i]
#   samples_processed <- list.files(sprintf("%s/tmp",io$out.folder))
#   if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
#     cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
#   } else {
#     cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
#     
#     # Read and parse raw methylation data
#     dat_sample <- fread(sprintf("zcat < %s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) %>%
#       .[,c("chr","pos","rate")]
#     
#     # Add 'start' and 'end' columns to do the overlap
#     dat_sample <- dat_sample[,c("start","end") := list(pos,pos)][,pos:=NULL] %>% .[,chr:=as.factor(chr)] %>% setkey(chr,start,end)
#     
#     # Overlap data with annotations
#     for (anno in opts$anno) {
#       fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
#       if (file.exists(paste0(fname.out,".gz"))) {
#         cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
#       } else {
#         cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
#         
#         # Overlap
#         setkey(anno_list[[anno]], chr, start, end)
#         ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% .[,"i.end":=NULL] %>% setnames("i.start","pos")
#         
#         # Calculate methylation status for each region in the annotation by summarising over all CG sites
#         out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), N=.N), keyby=.(sample,id,anno)]
#         
#         # Store and save results
#         fwrite(out, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
#         system(sprintf("gzip -f %s",fname.out))
#       }
#     }
#   }
# })
# 
# 
# 
# # Concatenate everything and save it
# for (i in opts$anno) {
#   cat(i)
#   # files <- list.files(str_c(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
#   system(sprintf("cat %s/tmp/*_%s.tsv.gz | zgrep -v 'sample' | gzip > %s/%s.tsv.gz", io$out.folder, i, io$out.folder, i))
#   # system(sprintf("cat %s/tmp/*_%s.tsv.gz > %s/%s.tsv.gz", io$out.folder, i, io$out.folder, i))
#   cat("\n")
# }
