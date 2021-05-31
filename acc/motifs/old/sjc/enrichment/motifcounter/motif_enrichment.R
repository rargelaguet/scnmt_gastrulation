suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(motifcounter))
suppressPackageStartupMessages(library(MotifDb))
suppressPackageStartupMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--anno',       type="character",                 help='genomic contexts (i.e. genebody, promoters, etc.')
p$add_argument('--min.diff',   type="double",                    help='Minimum differential rate(%)')
p$add_argument('--order',      type="integer",     default=1L,   help='order of the markov model')
p$add_argument('--lineages_foreground',    type="character",                help='Lineages for the foreground set')
p$add_argument('--lineages_background',    type="character",                help='Lineages for the background set')
p$add_argument('--test',       action="store_true",              help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## TESTING ##
# args <- list()
# args$anno <- c("H3K4me3")
# args$test <- TRUE
# args$min.diff <- 25
# args$condition <- "D3"
# args$order <- 1L
## TESTING ##

# stopifnot(args$lineage%in%c("Mesoderm","Endoderm","Ectoderm"))

###################
## Sanity checks ##
###################

cat(sprintf("Running motif enrichment for: %s\n",paste(args$anno,collapse=" + ")))

#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/met/results/differential")
io$outdir <- paste0(io$basedir,"/met/results/differential/motif_enrichment")
dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts$threshold_fdr <- 0.01

##############################
## Load genomic coordinates ##
##############################

feature_metadata <- lapply(args$anno, function(i) 
  fread(sprintf("%s/%s.bed.gz",io$features.dir,i))[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno")) %>%
  .[,chr:=paste0("chr",chr)]

###########################################
## Load differential methylation results ##
###########################################

dt <- args$anno %>%
  map(~ fread(sprintf("%s/%s_vs_%s_%s.txt.gz",io$indir,args$lineages_foreground,args$lineages_background,.))) %>%
  rbindlist

if (all(c("chr","start","end")%in%colnames(dt))) {
  dt <- dt %>% .[,chr:=paste0("chr",chr)]
}  else {
  dt <- dt %>% .[,c("id","anno","diff","padj_fdr","sig")] %>%
    merge(feature_metadata, by=c("id","anno"))
}

###########################################
## Define foreground and background sets ##
###########################################

# Option 1
# if (args$condition =="UD") {
#   dt.foreground <- dt[sig==TRUE & diff<(-args$min.diff)] %>% setorder(chr,start,end)
# } else if (args$condition =="D3") {
#   dt.foreground <- dt[sig==TRUE & diff>args$min.diff] %>% setorder(chr,start,end)
# }
# dt.background <- dt[sig==FALSE & abs(diff)<args$min.diff] %>% setorder(chr,start,end)

# Option 2
if (args$condition =="UD") {
  dt.foreground <- dt[sig==TRUE & diff<(-args$min.diff)] %>% setorder(chr,start,end)
  dt.background <- dt[sig==TRUE & diff>args$min.diff] %>% setorder(chr,start,end)
} else if (args$condition =="D3") {
  dt.foreground <- dt[sig==TRUE & diff>args$min.diff] %>% setorder(chr,start,end)
  dt.background <- dt[sig==TRUE & diff<(-args$min.diff)] %>% setorder(chr,start,end)
}

#######################
## Extract sequences ##
#######################

foreground.seq <- getSeq(Mmusculus, dt.foreground$chr, dt.foreground$start, dt.foreground$end)
names(foreground.seq) <- dt.foreground$id

background.seq <- getSeq(Mmusculus, dt.background$chr, dt.background$start, dt.background$end)
names(background.seq) <- dt.background$id

# Remove sequences with Ns
foreground.seq <- foreground.seq[alphabetFrequency(foreground.seq)[,"N"]==0]
background.seq <- background.seq[alphabetFrequency(background.seq)[,"N"]==0]

# Save sequences in FASTA format
# outfile <- sprintf("%s/%s_seq_foreground.fa",io$outdir,paste(args$anno,collapse="_"))
# writeXStringSet(foreground.seq, outfile, compress=FALSE)
# outfile <- paste0(io$outdir,"/seq_background.fa.gz")
# writeXStringSet(background.seq, outfile, compress=TRUE)

######################
## Motif enrichment ##
######################

# Load data base
motifs <- query(query(MotifDb, 'mmusculus'),"jaspar2018")

# [TEST MODE] Subset data base
if (isTRUE(args$test)) { 
  print("Testing mode, reducing the number of motifs...")
  motifs <- head(motifs,n=5)
}

# define motifcounter settings
motifcounterOptions(alpha=1e-4)
motifcounterOptions(ncores=1)

# read background sequence
bg <- readBackground(background.seq, order=args$order)

# compute motif enrichment in foreground sequence 
enrich.dt <- names(motifs) %>% map(function(x) {
  out <- tryCatch(
    motifEnrichment(foreground.seq, normalizeMotif(motifs[[x]]), bg), 
    error=function(e) list(pvalue=NA,fold=NA)
  )
  data.table(motif = x, pvalue = out$pvalue, fold = out$fold)
}) %>% rbindlist %>%
  setorder(pvalue, na.last=T) %>%
  .[,padj_fdr:=p.adjust(pvalue, method="fdr")] %>%
  .[,sig:=padj_fdr<=opts$threshold_fdr]

##########
## Save ##
##########

outfile <- sprintf("%s/%s_%s_%d.txt.gz", io$outdir, args$condition, paste(args$anno,collapse="_"), args$order)
fwrite(enrich.dt, outfile)
