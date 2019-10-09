##################################
## Script to binarise CpG sites ##
##################################

## Input: non-binarised methylation file, with rates that range from 0-100 or 0-1
# opts$input_format=1
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 50

# opts$input_format=2
# chr     pos     met_reads	nonmet_reads
# 1       3019021 1	1
# 1       3027398 0	1
# 1       3052955 5	10

## Output: binarised methylation files, rates are either 0 or 100.
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 0

suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
# p <- ArgumentParser(description='')
# p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
# p$add_argument('-n','--cores', type="integer",help='Number of cores')
# opts <- p$parse_args(commandArgs(TRUE))

opts <- list()
opts$context <- "CG"
opts$cores <- 4

opts$context <- toupper(opts$context); stopifnot(opts$context %in% c("CG","GC"))

# Define options
opts$input_format <- 2
opts$remove50 <- TRUE # if TRUE, sites with methylation rate of 50% are removed, otherwise they are rounded to 100%

# Define I/0
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  if (opts$context=="CG") {
    io$indir <- ""
    io$outdir <- ""
  } else {
    io$indir <- ""
    io$outdir <- ""
  }  
} else {
  if (opts$context=="CG") {
    io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/E5.5_scBS/2683_2684/merged"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/E5.5_scBS/2683_2684/merged/binarised"
  } else {
    io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw/filtered"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw/filtered/binarised"
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
    cat(sprintf("Sample %s already processed, skipping...\n",samples[i]))
  } else {
    print(sprintf("%s (%d/%d)", samples[i], i, length(samples)))
    
    # Load data
    data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F)
    
    # Input format 1 (chr,pos,rate)
    if (opts$input_format == 1) {
      colnames(data) <- c("chr","pos","rate")
      # data[,rate:=round(rate*100)] # if rate goes from 0 to 1...
    }
    
    # Input format 2 (chr,pos,met_reads,nonmet_reads)
    if (opts$input_format == 2) {
      colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
      data[,rate:=round((met_reads/(met_reads+nonmet_reads))*100)]# %>% .[,c("met_reads","nonmet_reads"):=NULL]
    }
    
    # Sanity check
    tmp <- sum((max(data$rate) > 100) | (min(data$rate) < 0))
    if (tmp>0) cat(sprintf("%s: There are %d CpG sites that have methylation rate higher than 100 or lower than 0\n",samples[i],tmp))
    
    # Deal with uncertain sites with rate=50
    if (opts$remove50) {
      data <- data[rate!=50,]
    } else {
      data[rate==50,rate:=100]
    }
    
    # Calculate binary methylation status
    cat(sprintf("%s: There are %0.03f%% of sites with non-binary methylation rate\n", samples[i], 100*mean(!data$rate %in% c(0,100))))
    data[,rate:=round(rate/100)]
    
    # Save results
    fwrite(data, file=outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
  }
}

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
# system(sprintf("pigz -p %d -f %s/*.tsv", opts$cores, io$outdir))