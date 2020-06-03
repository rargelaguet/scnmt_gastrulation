suppressMessages(library(SingleCellExperiment))

scale <- function(value, min.data, max.data, min.scaled, max.scaled) {
  stopifnot(is.numeric(value))
  stopifnot(value<=max.data & value>=min.data)
  return ((max.scaled - min.scaled) * (value - min.data) / (max.data - min.data)) + min.scaled
}

load_data <- function(io, rna.id, met.id, met.anno, min.cpg=1) {
  
  # Load DNA methylation data
  met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,met.anno)) %>%
    setnames(c("id_met","id","anno","Nmet","N","value")) %>%
    .[id%in%met.id] %>% .[N>=min.cpg]
  
  # Load RNA data
  sce <- readRDS(io$rna)
  sce <- sce[rowData(sce)$symbol == rna.id,]
  rna_dt <- logcounts(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
    melt(id.vars = "id_rna", value.name = "value", variable.name = "id")
  rna_dt$id <- rna.id
  
  
  return(list("met"=met_dt, "rna"=rna_dt))
}