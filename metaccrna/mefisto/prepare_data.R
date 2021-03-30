###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/load_settings.R")
} else {
  stop()
}

###################
## Load RNA data ##
###################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = opts$rna_cells, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Feature selection
decomp <- scran::modelGeneVar(sce)
hvgs <- rownames(decomp)[decomp$mean > 0.1 & decomp$p.value <= 0.10]

# Convert to data.table
rna_dt <-  as.matrix(logcounts(sce[hvgs,])) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

# Regress out technical covariates: number of expressed genes
foo <- data.table(id_rna=colnames(sce), covariate=sce$total_features_by_counts/nrow(sce))
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]], by="gene"] %>%
  .[,covariate:=NULL]

# Regress out technical covariates: batch effect between plates
foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>% .[,plate:=as.factor(grepl("PS_VE",plate))]
rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]], by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)
foo <- sample_metadata[stage=="E6.5" & plate%in%c("E6.5_late_Plate1","E6.5_late_Plate2"),c("id_rna","plate")]
rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]], by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)

###########################
## Load methylation data ##
###########################

met_dt <- lapply(opts$met.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n), showProgress=F) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by minimum number of cells
met_dt <- met_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$met_min.cells] %>% .[,ncells:=NULL]

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Rename annotations
met_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

# Regress out technical covariates 
foo <- fread(io$met.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_met","mean")]
met_dt <- met_dt %>% merge(foo, by="id_met") %>%
  .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(opts$acc.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n), showProgress=F) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nacc","N","rate"))

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt <- acc_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$acc_min.cells] %>% .[,ncells:=NULL]

# Calculate M value from Beta value
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Rename annotations
acc_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

# Regress out technical covariates 
foo <- fread(io$acc.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_acc","mean")]
acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
  .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

##############################
## Merge data with metadata ##
##############################

met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage")], by="id_met")
acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage")], by="id_acc")
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage")], by="id_rna")

###########################
## Prepare data for MOFA ##
###########################

data1 <- rna_dt %>% .[,c("sample","gene","expr")] %>%
  setnames(c("sample","feature","value")) %>% .[,c("view"):="RNA"]

data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,view:=paste0("met_",view)] %>%
  .[,feature:=paste0("met_",feature)]

data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,view:=paste0("acc_",view)] %>%
  .[,feature:=paste0("acc_",feature)]

# data <- do.call("rbind",list(data1,data2,data3)) %>%
data <- do.call("rbind",list(data1,data2)) %>%
  .[,value:=round(value,4)]

############################
## load pre-computed UMAP ##
############################

umap.dt <- fread("/Users/ricard/data/gastrulation/metaccrna/mofa/all_stages/umap_coordinates.txt") %>%
  .[sample%in%unique(data$sample)]
  # merge(sample_metadata[,c("id_rna","sample")]) %>% .[,sample:=NULL] %>%
  # .[id_rna%in%colnames(sce_filt)]

sample_metadata <- sample_metadata %>% merge(umap.dt,by="sample")

################
## Downsample ##
################

set.seed(42)
cells <- sample(unique(data$sample), size=1000)
data.downsampled <- data[sample%in%cells]

##########
## Save ##
##########

# file <- paste0(io$outdir,"/data.txt.gz")
# fwrite(data, file, col.names=T, quote=F, sep="\t")
