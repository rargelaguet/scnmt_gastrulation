matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))

#############################
## Load methylation data ##
#############################

met_dt <- lapply(opts$met.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n), showProgress=F) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

###################
## Load RNA data ##
###################

sce <- readRDS(io$rna)[,opts$rna_cells]

# Filter genes with low counts
sce <- sce[rowMeans(logcounts(sce))>1e-2,]

# Filter genes by variance
sce <- sce[apply(logcounts(sce),1,var)>0.5,]

# Convert to data.table
rna_dt <-  as.matrix(logcounts(sce)) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id", variable.factor = FALSE)# %>%
  # merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

##############################
## Merge data with metadata ##
##############################

met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage")], by="id_met")
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage")], by="id_rna")

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by minimum number of cells
met_dt <- met_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$met_min.cells] %>% .[,ncells:=NULL]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Filter genes by variance
keep_hv_genes <- rna_dt %>% .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id 
rna_dt <- rna_dt[ens_id%in%keep_hv_genes]

###########################
## Prepare data for MOFA ##
###########################

data1 <- rna_dt %>% .[,c("sample","ens_id","expr")] %>%
  setnames(c("sample","feature","value")) %>% .[,view:="RNA"]

data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,view:=paste0("met_",view)]

data <- do.call("rbind",list(data1,data2)) %>%
  .[,value:=round(value,4)]

