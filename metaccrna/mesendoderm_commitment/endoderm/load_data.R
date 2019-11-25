
##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_rna %in% opts$rna_cells]

###############
## Load data ##
###############

# Load Methylation data
met_dt <- lapply(names(opts$met.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met.dir,n)) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Load Accessibility data
acc_dt <- lapply(names(opts$acc.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc.dir,n)) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))


##############################
## Merge data with metadata ##
##############################

acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage","lineage10x_2","plate")], by="id_acc")
met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage","lineage10x_2","plate")], by="id_met")

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by  minimum number of cells
# met_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by variance
# met_dt <- met_dt[,var:=var(rate), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
# acc_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by variance
# acc_dt <- acc_dt[,var:=var(rate), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()
