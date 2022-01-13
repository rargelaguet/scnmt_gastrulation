suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))

here::i_am("rna/mapping/run/mapping_mnn.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_stages',   type="character",   nargs='+',  help='Query stage(s)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--use_marker_genes',action = "store_true",          help='Use marker genes?')
p$add_argument('--cosine_normalisation',      action = "store_true",          help='Use cosine normalisation?')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# I/O
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir


# Load mapping functions
source(here::here("rna/mapping/run/mapping_functions.R"))


## START TEST ##
# args$atlas_stages <- c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")
# args$query_stages <- c("E6.5","E7.5","E8.5")
# args$query_sce <- paste0(io$basedir,"/processed/rna_new/SingleCellExperiment.rds")
# args$atlas_sce <- io$atlas.sce
# args$query_metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
# args$atlas_metadata <- io$atlas.metadata
# args$test <- TRUE
# args$npcs <- 50
# args$n_neighbours <- 25
# args$use_marker_genes <- FALSE
# args$cosine_normalisation <- FALSE
# args$outfile <- paste0(io$basedir,"/results/rna/mapping/test.txt.gz")
## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% 
  .[pass_rnaQC==TRUE & stage%in%args$query_stages] %>% 
  .[,cell:=NULL] %>% setnames("id_rna","cell")
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_query %>% .[cell%in%colnames(sce_query)] %>% setkey(cell) %>% .[colnames(sce_query)]
stopifnot(tmp$cell == colnames(sce_query))
colData(sce_query) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_query),] %>% DataFrame()

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages] %>%
  .[,sample:=factor(sample)]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, normalise = TRUE, cells = meta_atlas$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_atlas %>% .[cell%in%colnames(sce_atlas)] %>% setkey(cell) %>% .[colnames(sce_atlas)]
stopifnot(tmp$cell == colnames(sce_atlas))
colData(sce_atlas) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_atlas),] %>% DataFrame()

############
## Filter ##
############

# Remove ExE cell types

#############
## Prepare ## 
#############

# Rename ensemble IDs to gene names in the atlas
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_atlas) <- foo[rownames(sce_atlas)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

#####################
## Define gene set ##
#####################

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Filter some genes manually
genes.intersect <- genes.intersect[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm",genes.intersect,invert=T)]
genes.intersect <- genes.intersect[!genes.intersect=="Xist"]
genes.intersect <- genes.intersect[!genes.intersect%in%gene_metadata[chr=="chrY",symbol]]

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#######################
## Feature selection ##
#######################

if (args$use_marker_genes) {
  # Load marker genes
  marker_genes.dt <- fread(io$rna.atlas.marker_genes)
  genes_to_use <- genes.intersect[genes.intersect%in%unique(marker_genes.dt$gene)]
} else {
  # Load gene statistics from the atlas
  gene_stats.dt <- fread(paste0(io$atlas.basedir,"/results/gene_statistics/gene_statistics.txt.gz")) %>%
    .[gene%in%genes.intersect]
  genes_to_use <- gene_stats.dt %>% setorder(-var_pseudobulk, na.last = T) %>% head(n=2500) %>% .$gene  
  
  # Calculate mean-variance relationship and extract HVGs
  # decomp <- modelGeneVar(sce_atlas, block=sce_atlas$sample)
  # genes_to_use <- rownames(decomp)[decomp$p.value<=0.01 & decomp$mean>0.1]
}

stopifnot(genes_to_use%in%rownames(sce_atlas))
stopifnot(genes_to_use%in%rownames(sce_query))

#########
## Map ##
#########

mapping  <- mapWrap(
  sce_atlas = sce_atlas,
  meta_atlas = meta_atlas,
  sce_query = sce_query,
  meta_query = meta_query,
  genes = genes_to_use,
  npcs = args$npcs,
  k = args$n_neighbours,
  cosineNorm = args$cosine_normalisation,
  order = NULL
)


##########
## Save ##
##########

mapping.dt <- mapping$mapping %>% 
  .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% 
  as.data.table

fwrite(mapping.dt, args$outfile, sep="\t")

##########
## TEST ##
##########

# previous_mapping.dt <- fread("/Users/argelagr/data/gastrulation_histones/results/rna/mapping/old/sample_metadata_after_mapping.txt.gz") %>%
#   .[,c("cell","celltype.mapped_mnn","celltype.mapped_seurat")] %>% setnames(c("cell","old_mapping_mnn","old_mapping_seurat"))
# foo <- merge(previous_mapping.dt, mapping.dt[,c("cell","celltype.mapped")], by="cell")
# 