#######################################################################
## Script to compute differential accessibility at the feature level ##
#######################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(argparse)
# source("/homes/ricard/gastrulation/met/differential/utils.R")

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-a',  '--anno',           type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('-s1', '--stage_lineage1', type="character",  nargs='+',  help='stage_lineage 1 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-s2', '--stage_lineage2', type="character",  nargs='+',  help='stage_lineage 2 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-cells', '--min.cells',   type="integer",                help='Minimum number of cells per group')
p$add_argument('-o',  '--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata_scNMT.txt")
io$data.dir <- paste0(io$basedir,"/acc/parsed/fimo_motifs")
# io$annos_dir  <- paste0(io$basedir, "/features/filt")
io$outfile <- args$outfile
dir.create(dirname(io$outfile), recursive = TRUE, showWarnings=F)

## Define options ##
opts <- list()

# Define genomic contexts
opts$annos <- args$anno

# Define stage and lineage
opts$groupA <- args$stage_lineage1
opts$groupB <- args$stage_lineage2

# Filter by variability
opts$fraction.sites <- 1.0    # Fraction of sites to keep based on variance

# Filter by coverage
opts$min.GpCs <- 3            # Minimum number of GpC per feature in each cell
opts$min.cells <- args$min.cells  # Minimum number of cells per feature in each group

# Multiple testing correction
opts$threshold_fdr <- 0.10

# Define which cells to use
opts$cells <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[pass_accQC==T & stage_lineage%in%c(opts$groupA,opts$groupB),id_acc]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample.metadata) %>% 
  .[id_acc%in%opts$cells] %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")]

# Load accessibility data
data <- lapply(opts$annos, function(n) fread(sprintf("zcat < %s/%s.tsv.gz",io$data.dir,n), showProgress=F)) %>% rbindlist
colnames(data) <- c("id","anno","rate","N","sample")

##############################
## Parse accessibility data ##
##############################

# Calculate Nacc
data[,Nacc:=round((rate/100)*N)]

# Merge accessibility data and sample metadata
data <- data %>% merge(sample_metadata[,c("sample","stage","stage_lineage")], by="sample") %>% setkey(anno)

# Define the two exclusive groups
data[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]
sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]

# Convert beta value to M value
# data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

###############################
## Filter accessibility data ##
###############################

# Remove sex chromosomes because it might have inherently large variation due to gender
# data <- merge(data,feature_metadata[,c("id","chr","anno")], by=c("id","anno")) %>% .[!chr %in% c("Y")] %>% .[,chr:=NULL]

# Filter features by coverage
data <- data[N>=opts$min.GpCs]

# Filter features by minimum number of cells per group
remove_n_sites <- data %>% split(.$anno) %>% map(~ .[,.(N=min(.N)), by=c("id","group")] %>% .[N<opts$min.cells,id])
data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[!id %in% remove_n_sites[[y]]]) %>% rbindlist

# Filter by variance
keep_hv_sites <- data %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% setorder(-var)  %>% head(n = nrow(.) * opts$fraction.sites) %>% .$id)
data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

#########################################
## Differential accessibility analysis ##
#########################################

# Binomial assumption: test of equal proportions using Fisher exact test
diff <- data[, .(
  A_acc=sum(.SD[group=="A",Nacc]), A_unacc=sum(.SD[group=="A",N-Nacc]),
  B_acc=sum(.SD[group=="B",Nacc]), B_unacc=sum(.SD[group=="B",N-Nacc]),
  p.value = fisher.test(
    x = matrix( c(
      A_acc=sum(.SD[group=="A",Nacc]), A_unacc=sum(.SD[group=="A",N-Nacc]),
      B_acc=sum(.SD[group=="B",Nacc]), B_unacc=sum(.SD[group=="B",N-Nacc])
    ), nrow = 2, ncol = 2))[["p.value"]]
), by = c("id","anno")]

# Multiple testing correction and define significant hits
diff[,c("prop1","prop2"):=list(A_acc/(A_acc+A_unacc), B_acc/(B_acc+B_unacc))] %>% 
  .[,diff:=prop2-prop1]  %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr")), by="anno"] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>0.10)] %>% 
  .[,c("prop1","prop2","diff","p.value","padj_fdr","log_padj_fdr"):=list(round(prop1,4),round(prop2,4),round(diff,4),round(p.value,10),round(padj_fdr,10),round(log_padj_fdr,5))] %>%
  setorderv("padj_fdr")

##################
## Save results ##
##################

write.table(diff, file=io$outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="NA")
