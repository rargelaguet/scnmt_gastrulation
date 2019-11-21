#######################################################################
## Script to compute differential accessibility at the feature level ##
#######################################################################

library(data.table)
library(purrr)
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
  source("/Users/ricard/gastrulation/met/results/differential/utils.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  source("/homes/ricard/gastrulation/met/results/differential/utils.R")
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$data.dir <- paste0(io$basedir,"/acc/feature_level")
io$annos_dir  <- paste0(io$basedir, "/features/genomic_contexts")
io$stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$outfile <- args$outfile

## Define options ##
opts <- list()

# Define genomic contexts
opts$annos <- args$anno

# Regress out global mean accessibility rate
# opts$regress.mean <- FALSE

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"
if (opts$regress.mean) {
  warning("Binomial test does not work when regressing out the global accessibility, switching to t.test")
  opts$statistical.test <- "t.test"
}

# Define stage and lineage
opts$groupA <- args$stage_lineage1
opts$groupB <- args$stage_lineage2

# Subset top most variable sites
opts$number_features <- 5000

# Filter by coverage
opts$min.GpCs <- 5                # Minimum number of GpC per feature in each cell
opts$min.cells <- args$min.cells  # Minimum number of cells per feature in each group

# Minimum differential accessibility (%) for statistical significance
opts$min.diff <- 5

# Multiple testing correction
opts$threshold_fdr <- 0.10

## START TESTING ##
# opts$annos <- "H3K27ac_distal_E7.5_Mes_intersect12"
# opts$groupA <- "E7.5_Mesoderm"; opts$groupB <- c("E7.5_Ectoderm")
# opts$min.cells <- 10
## END TESTING ##

# Define which cells to use
opts$cells <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_accQC==T & stage_lineage%in%c(opts$groupA,opts$groupB),id_acc]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample.metadata) %>% 
  .[id_acc%in%opts$cells] %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]

# Load accessibility data
data <- lapply(opts$annos, function(n) fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$data.dir,n), showProgress=F, header=F, sep="\t")) %>% 
  rbindlist %>% setnames(c("id_acc","id","anno","Nacc","N","rate"))

##############################
## Parse accessibility data ##
##############################

# Merge accessibility data and sample metadata
data <- data %>% merge(sample_metadata[,c("id_acc","stage","stage_lineage")], by="id_acc") %>% setkey(anno)

# Define the two exclusive groups
data[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]
sample_metadata[,group:=as.factor( c("A","B")[as.numeric(stage_lineage%in%opts$groupB)+1] )]

# Convert beta value to M value
if (opts$statistical.test=="t.test")
  data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

###############################
## Filter accessibility data ##
###############################

# Filter features by coverage
data <- data[N>=opts$min.GpCs]

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

# if (opts$regress.mean) {
  
#   # Load global statistics
#   stats <- fread(io$stats) %>% .[,mean:=mean/100]
  
#   # QC
#   stopifnot(all(stats$mean>0 & stats$mean<1))
#   stopifnot(all(data$id_acc %in% stats$id_acc))
  
#   # Calculate M-values from B-values
#   stats %>%  .[,covariate:=log2(((mean)+0.01)/(1-(mean/100)+0.01))]
  
#   # Merge data with global statistics
#   data <- merge(data, stats[,c("id_acc","covariate")], by="id_acc")
  
#   # Fit the linear model and regress out the covariate effect
#   data[, c("m"):=.(lm(formula=m~covariate)[["residuals"]]), by=c("id","anno")]
# }


#########################################
## Differential accessibility analysis ##
#########################################

# Binomial assumption: test of equal proportions using Fisher exact test
if (opts$statistical.test == "binomial") {
  diff <- data[, .(
    A_acc=sum(.SD[group=="A",Nacc]), A_unacc=sum(.SD[group=="A",N-Nacc]),
    B_acc=sum(.SD[group=="B",Nacc]), B_unacc=sum(.SD[group=="B",N-Nacc])), by = c("id","anno")] %>%
    .[,p.value := fisher.test(x = matrix( c(A_acc, A_unacc, B_acc, B_unacc), nrow=2, ncol=2))[["p.value"]], by=c("id","anno")] %>%
    .[,c("rateA","rateB"):=list(100*(A_acc/(A_acc+A_unacc)), 100*(B_acc/(B_acc+B_unacc)))]
  
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
