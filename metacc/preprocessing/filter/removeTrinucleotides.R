##################################################
## Script to filter NMT-specific trinucleotides ##
##################################################


"""
There are a few SNPs with respect to the reference genome and are therefore removed. Ideally one should do the filtering
directly at the read level using bismark with the option --NOMe
"""

# Load libraries
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

# Initialize argument parser
# p <- ArgumentParser(description='')
# p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
# p$add_argument('-n','--cores', type="integer",help='Number of cores')
# opts <- p$parse_args(commandArgs(TRUE))

opts <- list()
opts$context <- "CG"

# Define options
opts$context <- toupper(opts$context); stopifnot(opts$context %in% c("CG","GC"))

# Define I/0
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  if (opts$context=="CG") {
    io$indir <- "/Users/ricard/data/gastrulation/met/raw"
    io$outdir <- "/Users/ricard/data/gastrulation/met/raw/test/NMT_processed"
  } else {
    stop()
  }
} else {
  stop()
  if (opts$context=="CG") {
    # io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/met/raw"
    # io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/met/raw/filtered"
  } else {
    # io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw"
    # io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw/filtered"
  }
}
dir.create(io$outdir, showWarnings = F)

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
    data[,c("up1","center","down1") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]

    # Do sanity checks
    stopifnot(all(data$center%in%c("C","G")))
    # data %>% split(.$center) %>% walk2(.,names(.), function(x,y) if (y=="C") { stopifnot(all(x$down1=="G")) } else if (y=="G") { stopifnot(all(x$up1=="C")) })
    data[,trinucleotide:=ifelse(center=="C",paste0(up1,"C",down1),as.character(reverseComplement(DNAStringSet(paste0(up1,center,down1)))))]
      
    if (opts$context=="GC") {
      data <- data[trinucleotide%in%c("GCA","GCC","GCT")]
    } else if (opts$context=="CG") {
      data <- data[trinucleotide%in%c("ACG","TCG")]

    data[,c("up1","center","down1","trinucleotide"):=NULL]
    
    # Save results
    fwrite(data, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
    }
  }
}

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
# system(sprintf("pigz -p %d -f %s/*.tsv", opts$cores, io$outdir))
