###########################
## Script to extract mCH ##
###########################

SCRIPT IS UNFINISHED

# Input data:
# 2       3061990 3061990 0       0       1
# 2       3062007 3062007 0       0       1

# Output format:
# chr     pos     strand  context        met        total
# 12      3001502 -       CTT     0       1
# 12      3001504 -       CTC     0       1

# Load libraries
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

#####################
## Define settings ##
#####################

## I/0 ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  # io$basedir = 
} else {
  # io$basedir = 
}
dir.create(io$outdir, showWarnings = F)
  
## Options ##
opts <- list()
opts$cores <- 2


###############
## Load data ##
###############

# Load samples
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))
# samples <- sub(".tsv.gz","",list.files(io$indir,pattern="^(E7.5_Plate3).*(.tsv.gz)$"))

# Parallelise processing
# registerDoParallel(cores=opts$cores)
# invisible(foreach(i=1:length(samples)) %dopar% {
for (i in 1:length(samples)) {
  outfile <- sprintf("%s/%s.tsv",io$outdir,samples[i])
  if (file.exists(outfile)) {
    cat(sprintf("File %s already exists, skipping...\n",outfile))
  } else {
    
    ## Load data
    cat(sprintf("Processing %s...\n",samples[i]))
    data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F)
    
    ## Filter NMT trinucleotides

    # Get sequence
    seq <- getSeq(Mmusculus, paste0("chr",data$chr), start=data$pos-1, end=data$pos+1)
    data[,c("up","center","down") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]

    # Do sanity checks
    # stopifnot(all(data$center%in%c("C","G")))
    # data %>% split(.$center) %>% walk2(.,names(.), function(x,y) if (y=="C") { stopifnot(all(x$down1=="G")) } else if (y=="G") { stopifnot(all(x$up1=="C")) })
    # data[,trinucleotide:=ifelse(center=="C",paste0(up1,"C",down1),as.character(reverseComplement(DNAStringSet(paste0(up1,center,down1)))))]
    
    # if (opts$context=="CG") {
    #   data[,strand:=ifelse(base=="C","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
    # } else {
    #   data[,strand:=ifelse(base=="G","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
    # }
    # 
    # if (opts$context=="GC") {
    #   data <- data[trinucleotide%in%c("GCA","GCC","GCT")]
    # } else if (opts$context=="CG") {
    #   data <- data[trinucleotide%in%c("ACG","TCG")]
    # 
    # data[,c("up1","center","down1","trinucleotide"):=NULL]
    
    # Save results
    fwrite(data, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
    }
  }
}

# Compress output files
# system(sprintf("gzip -f %s/*.tsv",io$outdir))
system(sprintf("pigz -p %d -f %s/*.tsv", opts$cores, io$outdir))
