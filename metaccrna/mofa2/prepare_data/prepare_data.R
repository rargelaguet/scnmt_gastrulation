
#############################################################
## Script to prepare scNMT-seq gastrulation data for MOFA+ ##
#############################################################

####################
## Load libraries ##
####################

library(data.table)
library(purrr)
library(scater)

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/mofa2_gastrulation/load_settings.R")
} else {
  source("/homes/ricard/mofa2_gastrulation/load_settings.R")
}


#############################
## Load methylation data ##
#############################

met_dt <- lapply(opts$met.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met.dir,n), showProgress=F) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Rename annotations
met_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(opts$acc.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc.dir,n), showProgress=F) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nacc","N","rate"))

# Rename annotations
acc_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

###################
## Load RNA data ##
###################

sce <- readRDS(io$rna.file)[,opts$rna_cells]

# Re-calculate QC metrics before filtering
sce <- scater::calculateQCMetrics(sce) 

# Filter genes with low counts
sce <- sce[rowData(sce)$pct_dropout_by_counts < 90,]

# Re-calculate QC metrics after filtering
sce <- scater::calculateQCMetrics(sce)

# Convert to data.table
rna_dt <-  as.matrix(logcounts(sce)) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

############################
## Parse methylation data ##
############################

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

##############################
## Parse accessibility data ##
##############################

# Calculate M value from Beta value
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

##############################
## Merge data with metadata ##
##############################

met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage")], by="id_met")
acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage")], by="id_acc")
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage")], by="id_rna")

############################################################
## Regress out technical covariates in the RNA expression ##
############################################################

# Number of expressed genes
foo <- rna_dt[,.(covariate=sum(expr>0)), by=c("id_rna")]
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]]+mean(expr), by=c("gene","stage")] %>%
  .[,covariate:=NULL]

# Regress out batch effects in the RNA expression at E7.5
foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>%
  .[,plate:=as.factor(grepl("PS_VE",plate))]

rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]]+mean(expr), by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)

# Regress out batch effects in the RNA expression at E6.5
foo <- sample_metadata[stage=="E6.5",c("id_rna","plate")] %>%
  .[,plate:=as.factor(grepl("E6.5_late",plate))]

rna_dt <- rbind(
  rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
    .[,expr:=lm(formula=expr~plate)[["residuals"]]+mean(expr), by="gene"] %>% .[,plate:=NULL],
  rna_dt[!id_rna%in%foo$id_rna]
)

# Mithocondrial content
foo <- rna_dt[grepl("mt-",gene)] %>% .[,.(mt=sum(expr)), by="id_rna"]
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~mt)[["residuals"]]+mean(expr), by=c("gene","stage")] %>%
  .[,mt:=NULL]

############################################################
## Regress out technical variation in the DNA methylation ##
############################################################

# foo <- fread(io$met.stats) %>%
#   .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
#   .[,c("id_met","mean")]
# met_dt <- met_dt %>% merge(foo, by="id_met") %>%
#   .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage")]

#########################################################3##########
## Regress out technical variation in the chromatin accessibility ##
####################################################################

# Global accessibility rate (linked to the activity of the GpC methyltransferase)
# foo <- fread(io$acc.stats) %>%
#   .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
#   .[,c("id_acc","mean")]
# acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
#   .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage")]

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by minimum number of cells
met_dt <- met_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$met_min.cells] %>% .[,ncells:=NULL]

# Filter features by variance
# Important: if you select highly variable features using all cells, you will enrich for differences across stages
#   This is not what we want. Our aim is to look at the variation within stages, not between stages.
#   Here we will select highly variable features at E7.5 (the most interesting stage). 
#   Alternatively, we could also regress out differences between stages and then select highly variable features
keep_hv_sites <- met_dt[stage=="E7.5"] %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt <- acc_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$acc_min.cells] %>% .[,ncells:=NULL]

# Filter features by variance
# Important: if you select highly variable features using all cells, you will enrich for differences across stages
#   This is not what we want. Our aim is to look at the variation within stages, not between stages.
#   Here we will select highly variable features at E7.5 (the most interesting stage). 
#   Alternatively, we could also regress out differences between stages and then select highly variable features
keep_hv_sites <- acc_dt[stage=="E7.5"] %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()


################################
## Filter RNA expression data ##
################################

# Extract top N highly variable genes (per stage)
# keep_hv_genes <- rna_dt %>% split(.$stage) %>% 
#   map(~ .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var)  %>% head(n=opts$rna_ngenes) %>% .$ens_id) %>%
#   unlist %>% unname %>% as.character %>% unique
# rna_dt <- rna_dt[ens_id%in%keep_hv_genes]

# Important: if you select HVGs using all cells, you will enrich for differences across stages
#   This is not what we want. Our aim is to look at the variation within stages, not between stages.
#   Here we will select HVGs at E7.5 (the most interesting stage). 
#   Alternatively, we could also regress out differences between stages and then select HVGs
keep_hv_genes <- rna_dt[stage=="E7.5"] %>% .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id 
rna_dt <- rna_dt[ens_id%in%keep_hv_genes]

###########################
## Prepare data for MOFA ##
###########################

data1 <- rna_dt %>% .[,c("sample","stage","gene","expr")] %>%
  setnames(c("sample","group","feature","value")) %>% .[,c("view"):="RNA"]

data2 <- met_dt %>% .[,c("sample","stage","id","m","anno")] %>%
  setnames(c("sample","group","feature","value","view")) %>%
  .[,view:=paste0("met_",view)] %>%
  .[,feature:=paste0("met_",feature)]

data3 <- acc_dt %>% .[,c("sample","stage","id","m","anno")] %>%
  setnames(c("sample","group","feature","value","view")) %>%
  .[,view:=paste0("acc_",view)] %>%
  .[,feature:=paste0("acc_",feature)]

data <- do.call("rbind",list(data1,data2,data3)) %>%
  .[,value:=round(value,4)]

##########
## Save ##
##########

file <- paste0(io$outdir,"/data.txt")
fwrite(data, file, col.names=T, quote=F, sep="\t")
system(sprintf("pigz -f %s", file))
