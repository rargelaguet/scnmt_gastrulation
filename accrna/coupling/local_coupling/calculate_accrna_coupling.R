here::here("accrna/coupling/local_coupling/calculate_accrna_coupling.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--sce',  type="character",              help='SingleCellExperiment object')
p$add_argument('--tss',  type="character",              help='TSS annotation file')
p$add_argument('--up',  type="integer", default=3000,              help='Window upstream')
p$add_argument('--down',  type="integer", default=3000,              help='Window downstream')
p$add_argument('--window',  type="integer", default=3000,              help='Window size')
p$add_argument('--tile',  type="integer", default=100,              help='Tile size')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outfile',  type="character",              help='Output file')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################


## START TEST ##
# args <- list()
# args$sce <- io$rna.sce
# args$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$tss <- "/Users/argelagr/data/mm10_regulation/genes/TSS_protein_coding.bed"
# args$up <- 3000
# args$down <- 3000
# args$window <- 100
# args$tile <- 50
# args$test <- TRUE
# args$outfile  <- file.path(io$basedir,"results/accrna/coupling/precomputed_accrna_coupling.txt.gz")
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_accQC==TRUE & pass_rnaQC==TRUE]

# [TESTING MODE] Subset cells
if (args$test) {
  sample_metadata <- sample_metadata %>% .[,head(.SD,n=3),by="sample"]
}

#########################
## Load gene metadata  ##
#########################

tss <- fread(args$tss, colClasses=c("character","integer","integer","factor","factor","factor")) %>%
  setnames(c("chr","start","end","id","score","strand")) %>%
  .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
  .[,score:=NULL] %>%
  .[,tss:=start] %>%
  .[,c("start","end"):=.(start-args$up,end+args$down)] %>%
  setkey(chr,start,end)

# Define genomic window around TSS
tmp <- seq(from=0-args$up, to=0+args$down-args$window, by=args$tile)
foo <- data.table(window_center=tmp+(args$window/2), rel_start=tmp, rel_end=tmp+args$window) %>%
  .[,(c("window_center","rel_start","rel_end")):=lapply(.SD, as.integer),.SDcols=(c("window_center","rel_start","rel_end"))] %>%
  setkey(rel_start,rel_end)

####################
## Load RNA data  ##
####################

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$id_rna, normalise = TRUE)

# Add sample metadata as colData
colnames(sce) <- sample_metadata$cell
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)[,c("ens_id","symbol")] %>% 
  .[symbol%in%rownames(sce)] 

# Rename genes from symbol to ens_id
genes <- intersect(gene_metadata.dt$symbol,rownames(sce))
gene_metadata.dt <- gene_metadata.dt %>% .[symbol%in%genes] %>% setkey(symbol) %>% .[genes]
tmp <- gene_metadata.dt$ens_id; names(tmp) <- gene_metadata.dt$symbol
rownames(sce) <- tmp[rownames(sce)]

################################
## Calculate Met/RNA coupling ##
################################

cells.to.plot <- sample_metadata$cell# %>% head(n=5)

accrna_coupling.dt <- cells.to.plot %>% map(function(i) {

  print(i)
  
  # Define RNA expression data
  rna_dt <- data.table(
    id = rownames(sce),
    expr = logcounts(sce[,sample_metadata[cell==i,cell]])[,1]
  )

  # Define DNA methylation data
  acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,sample_metadata[cell==i,id_acc]), showProgress = F, header = T,
        select = c("chr"="character", "pos"="integer", "rate"="integer")) %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    setnames("pos","bp") %>% .[,c("start","end"):=list(bp,bp)] %>%
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(args$tile*round(dist/args$tile))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%

    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,.(rate=mean(rate)), by=c("id","window_center")]

    # Merge with RNA expression and calculate correlation coefficient per window
    merge(acc_dt, rna_dt, by="id") %>%
      .[,.(r = cor(x=rate, y=expr, method="pearson") %>% round(2)), by = c("window_center")] %>%
      .[,cell:=i] %>%
      return
}) %>% rbindlist

# print(length(unique(accrna_coupling.dt$cell)))

##########
## Save ##
##########

fwrite(accrna_coupling.dt, args$outfile)
