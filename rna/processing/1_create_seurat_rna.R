suppressPackageStartupMessages(library(Seurat))

here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                    help='Cell metadata file')
p$add_argument('--counts',        type="character",                    help='Counts file')
p$add_argument('--gene_metadata',        type="character",                    help='Gene metadata file')
p$add_argument('--outdir',       type="character",                    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$counts <- file.path(io$basedir,"processed/rna/counts_merged.txt.gz")
# args$metadata <- file.path(io$basedir,"processed/sample_metadata_merged.txt.gz")
# args$gene_metadata <- io$gene_metadata
# args$outdir <- file.path(io$basedir,"processed/rna")
## END TEST ##

#######################
## Load count matrix ##
#######################

rna_counts.mtx <- fread(args$count) %>% matrix.please

##########################
## Load sample metadata ##
##########################

metadata <- fread(args$metadata)

metadata <- metadata[id_rna%in%colnames(rna_counts.mtx)]
metadata$id_rna[!colnames(rna_counts.mtx)%in%metadata$id_rna]
metadata$id_rna[!metadata$id_rna %in% colnames(rna_counts.mtx)]

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(args$gene_metadata)[,c("ens_id","symbol")] %>% 
  .[ens_id%in%rownames(rna_counts.mtx) & symbol!=""] %>% .[!duplicated(symbol)]

########################################
## Rename genes from symbol to ens_id ##
########################################

genes <- intersect(gene_metadata.dt$ens_id,rownames(rna_counts.mtx))
gene_metadata.dt <- gene_metadata.dt %>% .[ens_id%in%genes] %>% setkey(ens_id) %>% .[genes]
rna_counts.mtx <- rna_counts.mtx[genes,] 
tmp <- gene_metadata.dt$symbol; names(tmp) <- gene_metadata.dt$ens_id
rownames(rna_counts.mtx) <- tmp[rownames(rna_counts.mtx)]

# Sanity checks
stopifnot(!is.na(rownames(rna_counts.mtx)))
stopifnot(!duplicated(rownames(rna_counts.mtx)))
stopifnot(sort(colnames(rna_counts.mtx))==sort(metadata$id_rna))

##################
## Filter genes ##
##################

# Remove duplicated genes
rna_counts.mtx <- rna_counts.mtx[!duplicated(rownames(rna_counts.mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(rna_counts.mtx)))==0)
stopifnot(sum(duplicated(colnames(rna_counts.mtx)))==0)

##########################
## Create Seurat object ##
##########################

metadata_to_seurat <- metadata %>% setkey(id_rna) %>% .[colnames(rna_counts.mtx)] %>% as.data.frame
rownames(metadata_to_seurat) <- metadata_to_seurat$id_rna
stopifnot(rownames(metadata_to_seurat)==colnames(rna_counts.mtx))

seurat <- CreateSeuratObject(rna_counts.mtx, meta.data = metadata_to_seurat)

# head(seurat@meta.data)

# Add mitochondrial percenatge
seurat[["mit_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["rib_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata.to.save <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL]
  
##########
## Save ##
##########

fwrite(metadata.to.save, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))

