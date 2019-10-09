#####################################################################
## Script to compute differential methylation at the feature level ##
#####################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(argparse)

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
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/Users/ricard/gastrulation/met/differential/utils.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene.metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/homes/ricard/gastrulation/met/differential/utils.R")
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$data.dir <- paste0(io$basedir,"/met/parsed")
io$annos_dir  <- paste0(io$basedir, "/features/filt")
io$stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$outfile <- args$outfile

## Define options ##
opts <- list()

# Define genomic contexts
opts$annos <- args$anno

# Define stage and lineage
opts$groupA <- args$stage_lineage1
opts$groupB <- args$stage_lineage2

# Overlap genomic features with nearby genes?
opts$OverlapWithGenes <- FALSE
opts$gene_window <- 5e4       # window length for the overlap

# Subset top most variable sites
opts$number_features <- 5000

# Filter by coverage
opts$min.CpGs <- 1            # Minimum number of CpG per feature in each cell
opts$min.cells <- args$min.cells  # Minimum number of cells per feature in each group

# Regress out global methylation rate
opts$regress.mean <- FALSE

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"
if (opts$regress.mean) {
  warning("Binomial test does not work when regressing out the global methylation, switching to t.test")
  opts$statistical.test <- "t.test"
}

# Minimum differential methylation (%) for statistical significance
opts$min.diff <- 5

# Multiple testing correction
opts$threshold_fdr <- 0.10

## START TESTING ##
# opts$annos <- "H3K27ac_distal_E7.5_End_intersect12"
# opts$groupA <- "E7.5_Endoderm"; opts$groupB <- c("E7.5_Ectoderm","E7.5_Mesoderm")
# opts$min.cells <- 10
## END TESTING ##

# Define which cells to use
opts$cells <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_metQC==T & stage_lineage%in%c(opts$groupA,opts$groupB),id_met]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample.metadata) %>% 
  .[id_met%in%opts$cells] %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]

# Load gene metadata
# gene_metadata <- fread(io$gene.metadata) %>% 
#   .[,chr:=as.factor(sub("chr","",chr))] %>%
#   setnames(c("ens_id","symbol"),c("id","gene"))

# Load genomic context metadata
# feature_metadata <- lapply(opts$annos, function(n) fread(sprintf("%s/%s.bed",io$annos_dir,n), showProgress=F)) %>% rbindlist
# colnames(feature_metadata) <- c("chr","start","end","strand","id","anno")

# Load methylation data
data <- lapply(opts$annos, function(n) fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$data.dir,n), showProgress=F, header=F)) %>% 
  rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

############################
## Parse methylation data ##
############################

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata[,c("id_met","stage","stage_lineage")], by="id_met") %>% setkey(anno)

# Define the two exclusive groups
data[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]
sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]

# Convert beta value to M value
if (opts$statistical.test == "t.test")
  data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#############################
## Filter methylation data ##
#############################

# Filter features by coverage
data <- data[N>=opts$min.CpGs]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by=c("id","anno")] %>% .[Ngroup==2] %>% .[,Ngroup:=NULL]

# Filter features by minimum number of cells per group
remove_n_sites <- data %>% split(.$anno) %>% map(~ .[,.(N=min(.N)), by=c("id","group")] %>% .[N<opts$min.cells,id])
data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[!id %in% remove_n_sites[[y]]]) %>% rbindlist

# Filter by variance
if (!is.na(opts$number_features)) {
  keep_hv_sites <- data %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% setorder(-var)  %>% head(n=opts$number_features) %>% .$id)
  data <- data %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist
}

###################################
## Regress out global statistics ##
###################################

if (opts$regress.mean) {

  # Load global statistics
  stats <- fread(io$stats) 
  
  # QC
  stopifnot(all(stats$mean>0 & stats$mean<1))
  stopifnot(all(data$id_acc %in% stats$id_acc))
  
  # Calculate M-values from B-values
  stats %>%  .[,covariate:=log2(((mean)+0.01)/(1-(mean/100)+0.01))]

  # Merge data with global statistics
  data <- merge(data, stats[,c("id_met","covariate")], by="id_met")

  # Fit the linear model and regress out the covariate effect
  data[, c("m"):=.(lm(formula=m~covariate)[["residuals"]]), by=c("id","anno")]
}

###############################################
## Associate the genomic features with genes ##
###############################################

if (opts$OverlapWithGenes==TRUE) {
  
  # Prepare feature metadata and gene metadata for the overlap
  gene_metadata_filt <- gene_metadata[, c("chr","start","end","gene","id")] %>%
    .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
    setkey(chr,start,end)
  
  feature_metadata_filt <- feature_metadata %>% split(.$anno) %>% 
    map2(.,names(.), function(x,y) x[id %in% data[anno==y,id]] ) %>%
    rbindlist
  
  # Do the overlap  
  data_list <- list()
  for (ann in unique(data$anno)){
    data_tmp <- data[anno == ann, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(data_tmp$id)) == FALSE)) {
      ov <- foverlaps(
        feature_metadata_filt[anno==ann, c("chr","start","end","id")] %>% setkey(chr,start,end),
        gene_metadata_filt[, c("chr","start","end","gene")],
        nomatch = NA) %>% .[,c("gene", "id")]
      
      # ov1 <- ov[is.na(gene)]
      # ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      # ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      data_list[[ann]] <- merge(ov, data_tmp, by = "id", allow.cartesian=T) %>%
        .[,c("id","gene","anno","id_met","rate","Nmet","N","stage_lineage","group")]
    }
    
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(data_tmp$id)) == TRUE)) {
      data_list[[ann]] <- merge(data_tmp, gene_metadata[, c("id", "gene")], by="id") %>%
        .[,c("id","gene","anno","id_met","rate","Nmet","N","stage_lineage","group")]
    }
  }
  data <- rbindlist(data_list)
  
} else {
  data[,gene:="NA"]
}


#######################################
## Differential methylation analysis ##
#######################################

# Binomial assumption: test of equal proportions using Fisher exact test
if (opts$statistical.test == "binomial") {
  diff <- data[, .(
    A_met=sum(.SD[group=="A",Nmet]), A_unmet=sum(.SD[group=="A",N-Nmet]),
    B_met=sum(.SD[group=="B",Nmet]), B_unmet=sum(.SD[group=="B",N-Nmet])), by = c("id","anno")] %>%
    .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id","anno")] %>%
    .[,c("rateA","rateB"):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))]
  
  # T-test under normality assumption
} else if (opts$statistical.test == "t.test") {
  diff <- data[, .(
    N_A = .SD[group=="A",.N], N_B = .SD[group=="B",.N],
    rateA = mean(.SD[group=="A",rate]), rateB = mean(.SD[group=="B",rate]),
    p.value = t.test(x=.SD[group=="B",m], y=.SD[group=="A",m], var.equal=FALSE)[["p.value"]]), by = c("id","anno")]
}

# Multiple testing correction and define significant hits
diff %>%
  .[,diff:=rateB-rateA] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr")), by="anno"] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,c("rateA","rateB","diff"):=list(round(rateA,2),round(rateB,2),round(diff,2))] %>%
  setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff, file=io$outfile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t", na="NA")
