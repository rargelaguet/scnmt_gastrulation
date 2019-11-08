###########################################################
## Script to do differential expression between lineages ##
###########################################################

suppressMessages(library(scater))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-s1', '--stage_lineage1', type="character",  nargs='+',  help='stage_lineage 1 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-s2', '--stage_lineage2', type="character",  nargs='+',  help='stage_lineage 2 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-o',  '--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))


## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/Users/ricard/gastrulation/rna/differential/utils.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene.metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/homes/ricard/gastrulation/rna/differential/utils.R")
}
io$sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$rna <- paste(io$basedir,"rna/SingleCellExperiment.rds",sep="/")
io$outfile <- args$outfile

## Options ##
opts <- list()

# Define stage and lineage
opts$groupA <- args$stage_lineage1
opts$groupB <- args$stage_lineage2
# opts$groupA <- c("E7.5_Ectoderm")
# opts$groupB <- c("E6.5_Epiblast")

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1.0

# Define which cells to use
opts$cells <- fread(io$sample_metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[pass_rnaQC==T & stage_lineage%in%c(opts$groupA,opts$groupB),id_rna]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample_metadata, stringsAsFactors=T) %>% 
  .[id_rna %in% opts$cells] %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]

# Load SingleCellExperiment object
sce <- readRDS(io$rna)[,opts$cells]
sce$stage_lineage <- sample_metadata$stage_lineage

# Define the two exclusive groups
sample_metadata[,group:=as.factor(as.numeric(stage_lineage%in%opts$groupB))]
sce$group <- as.factor(as.numeric(sce$stage_lineage%in%opts$groupB))

# Load gene metadata
gene_metadata <- rowData(sce) %>% as.data.frame(row.names=rownames(sce)) %>% 
  tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% 
  as.data.table %>% setnames("ens_id","id")
gene_metadata[,c("symbol","id"):=list(as.factor(symbol),as.factor(id))]


################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce) %>%
  merge(gene_metadata, all.y=T, by=c("id")) %>% setorderv("padj_fdr", na.last=T)

##################
## Save results ##
##################

fwrite(out, file=io$outfile, sep="\t", na="NA", quote=F)
