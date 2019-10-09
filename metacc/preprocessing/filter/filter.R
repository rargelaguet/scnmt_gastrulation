#####################################
## Script to filter specific sites ##
#####################################

# Current filters:
# - dinucleotides (non-cg methylation for example)
# - chromosomes

## Input: 
# single-cell methylation files output from Bismark. In either one of the two following formats:
# input_format=1:
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 100

# input_format=2:
# chr     pos     met_reads non nomet_reads
# 1       3019021 0 1
# 1       3027398 1 1
# 1       3052955 1 0

## Output:
# Same format as the input

# Load libraries
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')
opts <- p$parse_args(commandArgs(TRUE))


# Define options
opts$input_format <- 2
opts$chr_list <- c(1:19,"X","Y")
opts$context <- toupper(opts$context); stopifnot(opts$context %in% c("CG","GC"))

# opts$remove_dinucleotides can be a character vector of dinucleotides, NULL or 'non_cg'
# if (opts$context=="CG") {
#   opts$remove_dinucleotides <- "non_cg"
# } else {
#   opts$remove_dinucleotides <- "non_gc"
# }

# Define I/0
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  stop()
} else {
  if (opts$context=="CG") {
    io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/met/raw"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/met/raw/filtered"
  } else {
    io$indir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_NMT/acc/raw/filtered"
  }
}
dir.create(io$outdir)

cat("Options:\n")
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat(sprintf("- Valid chromosomes: %s\n",paste(opts$chr_list, collapse=" ")))
cat(sprintf("- Dinucleotides: %s\n",paste(opts$remove_dinucleotides,collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat("\n")

# Load samples
# samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="^(E7.5_Plate3).*(.tsv.gz)$"))

# Parallelise processing
# registerDoParallel(cores=opts$cores)
# invisible(foreach(i=1:length(samples)) %dopar% {
for (i in 1:length(samples)) {
  outfile <- sprintf("%s/%s.tsv",io$outdir,samples[i])
  if (file.exists(outfile)) {
    cat(sprintf("File %s already exists, skipping...\n",outfile))
  } else {
    # Load data
    cat(sprintf("Processing %s...\n",samples[i]))
    data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F)
    
    # Input format 1 (chr,pos,rate)
    if (opts$input_format == 1) {
      colnames(data) <- c("chr","pos","rate")
    # Input format 2 (chr,pos,met_reads,nonmet_reads)
    } else if (opts$input_format == 2) {
      colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
    }
    
    # Filter by chromosomes
    data <- data[chr %in% opts$chr_list,]
    
    # Filter by dinucleotides
    # if (length(opts$remove_dinucleotides)>0) {
    #     
    #   # Get sequence
    #   seq <- getSeq(Mmusculus, paste0("chr",data$chr), data$pos-1, data$pos+1)
    #   data[,c("base_up","base","base_down") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]
    #   
    #   # Do sanity checks
    #   stopifnot(unique(data$base) %in% c("G","C"))
    #   
    #   # Guess nucleotide taking into account strand info
    #   # data[,dinucleotide:=ifelse(base=="C",paste0(base,base_down),paste0(base_up,base))]
    #   data[,dinucleotide:=ifelse(base==substr(opts$context,1,1),paste0(base,base_down),paste0(base_up,base))]
    #     
    #   if (opts$remove_dinucleotides == "non_cg") {
    #     idx_keep <- which(data$dinucleotide=="CG")
    #     data <- data[idx_keep][,dinucleotide:=NULL]
    #     data %>% split(.$base) %>% walk2(.,names(.), function(x,y) if (y=="C") { stopifnot(all(x$base_down=="G")) } else if (y=="G") { stopifnot(all(x$base_up=="C")) })
    #     
    #   } else if (opts$remove_dinucleotides == "non_gc") {
    #     idx_keep <- which(data$dinucleotide=="GC")
    #     data <- data[idx_keep][,dinucleotide:=NULL]
    #     data %>% split(.$base) %>% walk2(.,names(.), function(x,y) if (y=="G") { stopifnot(all(x$base_down=="C")) } else if (y=="C") { stopifnot(all(x$base_up=="G")) })
    #     
    #   } else {
    #     idx_keep <- which(!data$dinucleotide %in% opts$remove_dinucleotides)
    #     data <- data[idx_keep][,dinucleotide:=NULL]
    #   }
    #   data[,c("base_up","base","base_down"):=NULL]
    # }
    
    # Save results
    fwrite(data, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=FALSE)
  }
}

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
# system(sprintf("pigz -p %d -f %s/*.tsv", opts$cores, io$outdir))