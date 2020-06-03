
###################
## Load settings ##
###################

source("/Users/ricard/scnmt_gastrulation/metaccrna/mofa2/endoderm/load_settings.R")

#############################
## Load methylation data ##
#############################

met_dt <- lapply(opts$met.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n), showProgress=F) %>% .[V1%in%sample_metadata$id_met]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Rename annotations
met_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(opts$acc.annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n), showProgress=F) %>% .[V1%in%sample_metadata$id_acc]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nacc","N","rate"))

# Rename annotations
acc_dt[,anno:=stringr::str_replace_all(anno,opts$rename.annos)]

# Calculate M value from Beta value
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

###################
## Load RNA data ##
###################

sce <- readRDS(io$rna)[,sample_metadata$id_rna]

# Re-calculate QC metrics before filtering
sce <- sce[rowMeans(logcounts(sce))>0.1,]

# Convert to data.table
rna_dt <-  as.matrix(logcounts(sce)) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

##############################
## Merge data with metadata ##
##############################

met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage")], by="id_met")
acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage")], by="id_acc")
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage")], by="id_rna")

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by minimum number of cells
met_dt <- met_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$met_min.cells] %>% .[,ncells:=NULL]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt <- acc_dt[,ncells:=.N, by=c("id","anno")] %>% .[ncells>=opts$acc_min.cells] %>% .[,ncells:=NULL]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()


################################
## Filter RNA expression data ##
################################

# Extract highly variable genes
keep_hv_genes <- rna_dt %>% .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id 
rna_dt <- rna_dt[ens_id%in%keep_hv_genes]

###########################
## Prepare data for MOFA ##
###########################

data1 <- rna_dt %>% .[,c("sample","gene","expr")] %>%
  setnames(c("sample","feature","value")) %>% .[,c("view"):="RNA expression"]

data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,view:=sprintf("%s methylation",view)] %>%
  .[,feature:=paste0("met_",feature)]

data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,view:=sprintf("%s accessibility",view)] %>%
  .[,feature:=paste0("acc_",feature)]

data <- do.call("rbind",list(data1,data2,data3)) %>%
  .[,value:=round(value,4)]

##############
## Run MOFA ##
##############

# Create MOFA object
MOFAobject <- create_mofa(data)

# Plot data overview
plot_data_overview(MOFAobject)

# prepare MOFA object
MOFAobject <- prepare_mofa(MOFAobject)

# Train the MOFA model
MOFAobject <- run_mofa(MOFAobject)

#########################
## Add sample metadata ##
#########################

cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))

# sample_metadata.mofa <- sample_metadata[,c("sample","id_rna","stage","lineage10x","pass_rnaQC","pass_metQC","pass_accQC")] %>%
#   setnames("lineage10x","lineage") %>%
#   .[sample%in%cells] %>% setkey(sample) %>% .[cells]
# samples_metadata(MOFAobject) <- sample_metadata.mofa

foo <- fread("/Users/ricard/data/gastrulation/rna/results/mapping/mapping.txt") %>%
  merge(sample_metadata.mofa[,c("sample","id_rna")], by="id_rna") %>%
  setnames("celltype.mapped","lineage") %>%
  setkey(sample) %>% .[cells]
samples_metadata(MOFAobject) <- foo


####################
## Subset factors ##
####################

r2 <- MOFAobject@cache$variance_explained$r2_per_factor
factors <- sapply(r2, function(x) x[,"RNA expression"]>1)
MOFAobject <- subset_factors(MOFAobject, which(apply(factors,1,sum)>=1))
factors_names(MOFAobject) <- paste("Factor",1:get_dimensions(MOFAobject)[["K"]], sep=" ")

#############################
## Plot variance explained ##
#############################

plot_variance_explained(MOFAobject, legend = T)

##################
## Plot factors ##
##################

plot_factors(MOFAobject, factors=c(1,2), color_by = "lineage", dot_size = 2.5)
plot_factors(MOFAobject, factors=c(1,2), color_by = "Rhox5", dot_size = 2.5)
# plot_factors(MOFAobject, factors=c(1,2), color_by = "lineage") + scale_fill_manual(values=opts$colors)


################
## Save model ##
################

saveRDS(MOFAobject, paste0(io$outdir,"/MOFAmodel.rds"))
