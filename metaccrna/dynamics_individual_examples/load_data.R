suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(scater))
suppressMessages(library(purrr))

load_data <- function(rna.id, met.id, met.anno, acc.id, acc.anno, min.cpg=1, min.gpc=1) {
  
  ################
  ## Define I/O ##
  ################
  
  io <- list()
  io$met.dir <- "/Users/ricard/data/gastrulation/met/parsed"
  io$acc.dir <- "/Users/ricard/data/gastrulation/acc/parsed"
  io$rna.file <- "/Users/ricard/data/gastrulation/rna/parsed/SingleCellExperiment.rds"
  
  ##########
  ## Load ##
  ##########
  
  # Load DNA methylation data
  met_dt <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,met.anno), stringsAsFactors=F, quote="", showProgress=F) %>%
    setnames(c("id_met","id","anno","Nmet","N","rate"))
  
  # Load DNA accessibility data
  acc_dt <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,acc.anno), stringsAsFactors=F, quote="", showProgress=F) %>%
    setnames(c("id_acc","id","anno","Nmet","N","rate"))
  
  # Load RNA data
  sce <- readRDS(io$rna.file)
  rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
    melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
    merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))
  
  ############
  ## Filter ##
  ############
  
  # Select ID
  met_dt <- met_dt[id%in%met.id] %>% setnames("rate","value")
  acc_dt <- acc_dt[id%in%acc.id] %>% setnames("rate","value")
  rna_dt <- rna_dt[ens_id%in%rna.id] %>% setnames("expr","value")
  
  # Filter by coverage
  met_dt <- met_dt[N>=min.cpg]
  acc_dt <- acc_dt[N>=min.gpc]
  
  return(list("met"=met_dt, "acc"=acc_dt, "rna"=rna_dt))
}