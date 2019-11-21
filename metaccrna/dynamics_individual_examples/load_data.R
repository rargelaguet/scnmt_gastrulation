suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(scater))
suppressMessages(library(purrr))

load_data <- function(rna.id, met.id, met.anno, acc.id, acc.anno, min.cpg=1, min.gpc=1) {
  
  ################
  ## Define I/O ##
  ################
  
  io <- list()
  io$met.dir <- "/Users/ricard/data/gastrulation/met/feature_level"
  io$acc.dir <- "/Users/ricard/data/gastrulation/acc/feature_level"
  io$rna.file <- "/Users/ricard/data/gastrulation/rna/SingleCellExperiment.rds"
  
  ##########
  ## Load ##
  ##########
  
  # Load DNA methylation data
  met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met.dir,met.anno)) %>%
    setnames(c("id_met","id","anno","Nmet","N","value")) %>%
    .[id%in%met.id]
  
  # Load DNA accessibility data
  acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc.dir,acc.anno)) %>%
    setnames(c("id_acc","id","anno","Nmet","N","value")) %>%
    .[id%in%acc.id]
  
  # Load RNA data
  sce <- readRDS(io$rna.file)[rna.id,]
  rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
    melt(id.vars = "id_rna", value.name = "value", variable.name = "ens_id") %>%
    merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","id"))
  
  ############
  ## Filter ##
  ############
  
  # Filter by coverage
  met_dt <- met_dt[N>=min.cpg]
  acc_dt <- acc_dt[N>=min.gpc]
  
  return(list("met"=met_dt, "acc"=acc_dt, "rna"=rna_dt))
}