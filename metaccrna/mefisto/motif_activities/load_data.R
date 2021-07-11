
###################
## Load RNA data ##
###################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = opts$rna_cells, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)
rownames(sce) <- rowData(sce)$symbol

# Remove genes manually
sce <- sce[grep("mt-",rownames(sce),invert = T),]

# Feature selection
decomp <- scran::modelGeneVar(sce)
hvgs <- rownames(decomp)[decomp$mean > 0.1 & decomp$p.value <= 0.10]

# Convert to data.table
rna_dt <-  as.matrix(logcounts(sce[hvgs,])) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "gene") 

# Regress out technical covariates: number of expressed genes
foo <- data.table(id_rna=colnames(sce), covariate=colMeans(counts(sce)>0))
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

#################################
## Load motif methylation data ##
#################################

met_dt <- fread(paste0(io$basedir,"/met/results/motifs/motif_met.txt.gz")) %>%
  setnames("anno","id") %>%
  .[,id:=gsub("multiome_peaks_jaspar2020_","",id)]

# Calculate M value from Beta value
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] 

# Regress out technical covariates 
foo <- fread(io$met.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_met","mean")]
met_dt <- met_dt %>% merge(foo, by="id_met") %>%
  .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id")]

# Filter features by variance
keep_hv_motifs <- met_dt[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id
met_dt <- met_dt[id %in% keep_hv_motifs] %>% droplevels

#############################
## Load accessibility data ##
#############################

acc_dt <- fread(paste0(io$basedir,"/acc/results/motifs/motif_acc.txt.gz")) %>%
  setnames("anno","id") %>%
  .[,id:=gsub("multiome_peaks_jaspar2020_","",id)]

# Calculate M value from Beta value
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Regress out technical covariates
foo <- fread(io$acc.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_acc","mean")]
acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
  .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id")]

# Filter features by variance
keep_hv_motifs <- acc_dt[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id
acc_dt <- acc_dt[id %in% keep_hv_motifs] %>% droplevels

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

data2 <- met_dt %>% .[,c("sample","id","m")] %>%
  setnames(c("sample","feature","value")) %>%
  .[,feature:=paste0(feature,"_met")] %>%
  .[,view:="motif_met"]

data3 <- acc_dt %>% .[,c("sample","id","m")] %>%
  setnames(c("sample","feature","value")) %>%
  .[,feature:=paste0(feature,"_acc")] %>%
  .[,view:="motif_acc"]

data <- do.call("rbind",list(data1,data2,data3)) %>%
# data <- do.call("rbind",list(data1,data2)) %>%
  .[,value:=round(value,4)]

################
## Downsample ##
################

# set.seed(42)
# cells <- sample(unique(data$sample), size=1000)
# data.downsampled <- data[sample%in%cells]

##########
## Save ##
##########

# file <- paste0(io$outdir,"/vignette/data.txt.gz")
# fwrite(data, file, col.names=T, quote=F, sep="\t")
