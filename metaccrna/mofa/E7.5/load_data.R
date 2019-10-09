
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

# Load RNA data
sce <- readRDS(io$rna.file) %>% .[,opts$rna_cells]

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

## Parse RNA expression data ##

# Convert to data.table
rna_dt <- exprs(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))
# rna_dt[,c("id_rna","gene","ens_id"):=list(as.factor(id_rna),as.factor(gene),as.factor(ens_id))]

## Parse accessibility data ##
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


##############################
## Merge data with metadata ##
##############################

acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc","stage","stage_lineage")], by="id_acc") %>% droplevels()
met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met","stage","stage_lineage")], by="id_met") %>% droplevels()
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna","stage","stage_lineage")], by="id_rna") %>% droplevels()


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

# Filter features by  minimum number of cells
met_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by  minimum number of cells (by stage_lineage)
# for (i in unique(met_dt$stage_lineage)) {
#   met_dt[stage_lineage==i,Ntotal:=sample_metadata[id_met%in%opts$met_cells & stage_lineage==i,.N]]
# }
# keep_cov_sites <- met_dt %>% split(.$stage_lineage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$met_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# met_dt <- met_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
met_dt <- met_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells
acc_dt[,N:=.N,by=c("id","anno","gene")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by  minimum number of cells (by stage_lineage)
# for (i in unique(acc_dt$stage_lineage)) {
#   acc_dt[stage_lineage==i,Ntotal:=sample_metadata[id_acc%in%opts$acc_cells & stage_lineage==i,.N]]
# }
# keep_cov_sites <- acc_dt %>% split(.$stage_lineage) %>% map(~ .[, ncells:=.N, by=c("id","anno","gene")] %>% .[ncells >= opts$acc_min.cells] %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>% .$id_anno)
# acc_dt <- acc_dt %>% .[,id_anno:=paste(as.character(id),as.character(anno),sep="_")] %>%
#   .[id_anno%in%Reduce("intersect",keep_cov_sites)] %>% .[,c("Ntotal","id_anno"):=NULL]

# Filter features by variance
acc_dt <- acc_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Remove lowly expressed genes
rna_dt <- rna_dt[,mean:=mean(expr),by="ens_id"] %>% .[mean>0.1] %>% .[,mean:=NULL]

# Remove genes with constant expression levels
rna_dt <- rna_dt[,var:=var(expr),by="ens_id"] %>% .[var>0.1] %>% .[,var:=NULL]

# Filter genes with low cellular detection rate and sites with low coverage across samples
rna_dt <- rna_dt[,cdr:=sum(expr>0)/length(opts$rna_cells), by="ens_id"] %>% .[cdr>=opts$rna_min.cdr] %>% .[,cdr:=NULL]

# Extract top N highly variable genes
rna_dt <- rna_dt[,var:=var(expr), by="ens_id"] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

############################
## Regress out covariates ##
############################

# RNA: number of expressed genes
foo <- data.table(id_rna=colnames(sce), covariate=sce$total_features_by_counts/nrow(sce))
rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
  .[,expr:=lm(formula=expr~covariate)[["residuals"]], by=c("gene")] %>%
  .[,covariate:=NULL]

# RNA: strong batch effects between E7.5 plates
if (any(sample_metadata$stage=="E7.5")) {
  foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>% 
    .[,plate:=as.factor(grepl("PS_VE",plate))]
  rna_dt <- rbind(
    rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
      .[,expr:=lm(formula=expr~plate)[["residuals"]], by=c("gene")] %>% .[,c("plate"):=NULL],
    rna_dt[!id_rna%in%foo$id_rna]
  )
}
if (any(sample_metadata$stage=="E7.5")) {
  foo <- sample_metadata[stage=="E7.5",c("id_rna","plate")] %>%
    .[,plate:=as.factor(plate%in%c("E7.5_Plate1","E7.5_Plate2","E7.5_Plate3","E7.5_Plate4"))]
  rna_dt <- rbind(
    rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
      .[,expr:=lm(formula=expr~plate)[["residuals"]], by=c("gene")] %>% .[,c("plate"):=NULL],
    rna_dt[!id_rna%in%foo$id_rna]
  )
}

# RNA: strong batch effect between E6.5 plates
# if (any(sample_metadata$stage=="E6.5")) {
#   foo <- sample_metadata[stage=="E6.5" & plate%in%c("E6.5_late_Plate1","E6.5_late_Plate2"),c("id_rna","plate")]
#   rna_dt <- rbind(
#     rna_dt[id_rna%in%foo$id_rna] %>% merge(foo, by="id_rna") %>%
#       .[,expr:=lm(formula=expr~plate)[["residuals"]], by=c("gene")] %>% .[,plate:=NULL],
#     rna_dt[!id_rna%in%foo$id_rna]
#   )
# }

# Methylation: differences in mean methylation rate
# foo <- met_dt[,.(covariate=mean(m)),by=c("id_met")]
foo <- fread(io$met.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_met","mean")]
met_dt <- met_dt %>% merge(foo, by="id_met") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

# Accessibility: differences in global accessibility rate
foo <- fread(io$acc.stats) %>%
  .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
  .[,c("id_acc","mean")]
acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
  .[,m:=mean(m) + lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

#####################################
## Select highly variable features ##
#####################################

# RNA: Extract top N highly variable genes
keep_hv_genes <- rna_dt[,.(var=var(expr)), by="ens_id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id
rna_dt <- rna_dt[ens_id%in%as.character(keep_hv_genes)] %>% droplevels()

# Accessibility: Extract top N most variable features
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

# Methylation: Extract top N most variable features
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

# met_dt[,var(m),by=c("id","anno")]

#############################
## Create joint data.frame ##
#############################

data1 <- rna_dt %>% .[,c("sample","gene","expr")] %>%  
  setnames(c("sample","feature","value")) %>% .[,c("feature_group","sample_group"):=list("RNA","gastrulation")]
data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("met_",feature), paste0("met_",feature_group), "gastrulation")]
data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
  setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group), "gastrulation")]

data <- rbind(data1,data2,data3)

outfile <- paste0(io$outdir,"/data.txt")
fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
system(sprintf("pigz -f %s",outfile))


###########################
## Create input matrices ##
###########################

met_cells <- as.character(unique(met_dt$sample))
rna_cells <- as.character(unique(rna_dt$sample))
acc_cells <- as.character(unique(acc_dt$sample))

rna_matrix <- rna_dt[,c("gene","expr","sample")] %>%
  .[,c("sample","gene"):=list(as.character(sample),as.character(gene))] %>%
  .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
  dcast(sample~gene, value.var="expr", drop=F) %>% matrix.please() %>% t

met_matrix_list <- list()
for (n in unique(met_dt$anno)) {
  met_matrix_list[[paste("met",n,sep="_")]] <- met_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F) %>% matrix.please() %>% t

  cat(sprintf("%s methylation matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(met_matrix_list[[paste("met",n,sep="_")]]), ncol(met_matrix_list[[paste("met",n,sep="_")]]),
              100*mean(is.na(met_matrix_list[[paste("met",n,sep="_")]]))))
}

cat("\n")

acc_matrix_list <- list()
for (n in unique(acc_dt$anno)) {
  acc_matrix_list[[paste("acc",n,sep="_")]] <- acc_dt[anno==n,c("id","gene","m","sample")] %>%
    .[,c("sample","gene","id"):=list(as.character(sample),as.character(gene),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
    .[,id_gene:=paste(id,gene,sep="_")] %>%
    dcast(sample~id_gene, value.var="m", drop=F) %>% matrix.please() %>% t

  cat(sprintf("%s accessibility matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(acc_matrix_list[[paste("acc",n,sep="_")]]), ncol(acc_matrix_list[[paste("acc",n,sep="_")]]),
              100*mean(is.na(acc_matrix_list[[paste("acc",n,sep="_")]]))))
}

all_matrix_list <- c(rna=list(rna_matrix),met_matrix_list,acc_matrix_list)
