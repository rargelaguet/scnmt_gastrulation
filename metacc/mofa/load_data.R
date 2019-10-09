library(data.table)
library(purrr)
library(ggplot2)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

###############
## Load data ##
###############

# Load Methylation data
met_dt <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Load Accessibility data
acc_dt <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))

# Load annotation metadata
feature_metadata <- lapply(unique(c(opts$met.annos,opts$acc.annos)), function(i) 
  fread(sprintf("%s/%s.bed",io$annos_dir,i), stringsAsFactors=T)[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

# Load gene metadata 
gene_metadata <- fread(io$gene_metadata,stringsAsFactors=T) %>% 
  setnames(c("ens_id","symbol"),c("id","gene")) %>% 
  .[,chr:=as.factor(sub("chr","",chr))]

################
## Parse data ##
################

# Parse gene and feature metadata
feature_metadata_filt.met <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% met_dt[anno==y,id]] ) %>%
  rbindlist
feature_metadata_filt.acc <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% acc_dt[anno==y,id]] ) %>%
  rbindlist

gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","gene")] %>% 
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)

## Parse accessibility data ##
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


##############################
## Merge data with metadata ##
##############################

acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()

###########################################################
## Associate the genomic features with overlapping genes ##
###########################################################
  
# Methylation
if (opts$overlapGenes) {
  met_list <- list()
  for (i in unique(met_dt$anno)){
    
    # Subset corresponding anno
    met_tmp <- met_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(met_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.met[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      met_list[[i]] <- merge(met_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(met_tmp$id)) == TRUE)) {
      met_list[[i]] <- merge(met_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  met_dt <- rbindlist(met_list)
  rm(met_list, met_tmp,feature_metadata_tmp,ov)
} else {
  met_dt[,gene:="NA"]
}

# Accessibility
if (opts$overlapGenes) {
  acc_list <- list()
  for (i in unique(acc_dt$anno)){
    
    # Subset corresponding anno
    acc_tmp <- acc_dt[anno==i, ]
    
    # Non gene-associated feature
    if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == FALSE)) {
      
      # Extract coordiantes for methylation sites and for genes
      feature_metadata_tmp <- feature_metadata_filt.acc[anno==i, c("chr","start","end","id")] %>% 
        .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
      
      # Do the overlap
      ov <- foverlaps(
        gene_metadata_filt, 
        feature_metadata_tmp, 
        nomatch=0) %>% .[,c("gene", "id")]
      
      # If a feature overlaps with multiple genes, collapse them
      ov1 <- ov[is.na(gene)]
      ov2 <- ov[!is.na(gene)] %>% .[,.(gene=paste(gene,collapse="_")), by="id"]
      ov <- rbind(ov1,ov2)
      
      # Merge with methylation data
      acc_list[[i]] <- merge(acc_tmp, ov, by="id", allow.cartesian=T) 
    }
    # Gene-associated feature
    else if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == TRUE)) {
      acc_list[[i]] <- merge(acc_tmp, gene_metadata[,c("id","gene")], by="id")
    }
  }
  acc_dt <- rbindlist(acc_list)
  rm(acc_list, acc_tmp,feature_metadata_tmp,ov)
} else {
  acc_dt[,gene:="NA"]
}

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by  minimum number of cells (by stage_lineage)
# for (i in unique(met_dt$stage_lineage)) {
#   met_dt[stage_lineage==i,Ntotal:=sample_metadata[id_met%in%opts$met_cells & stage_lineage==i,.N]]
# }
# keep_cov_sites <- met_dt %>% split(.$stage_lineage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$met_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# met_dt <- met_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]
met_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0.5] %>% setorder(-var) %>% head(n = opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells (by stage_lineage)
# for (i in unique(acc_dt$stage_lineage)) {
#   acc_dt[stage_lineage==i,Ntotal:=sample_metadata[id_acc%in%opts$acc_cells & stage_lineage==i,.N]]
# }
# keep_cov_sites <- acc_dt %>% split(.$stage_lineage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$acc_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# acc_dt <- acc_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]
acc_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0.5] %>% setorder(-var) %>% head(n = opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()


#############################
## Create joint data.frame ##
#############################

data.met <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("met_",feature), paste0("met_",feature_group), "gastrulation")]
data.acc <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group), "gastrulation")]

data <- rbind(data.met,data.acc)

outfile <- paste0(io$outdir,"/data.txt")
fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
system(sprintf("pigz -f %s",outfile))


# outfile <- paste0(io$outdir,"/data_acc.txt")
# fwrite(data.acc, file=outfile, col.names=T, quote=F, sep="\t")
# system(sprintf("pigz -f %s",outfile))
