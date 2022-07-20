suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))

here::i_am("rna/mapping/trajectories/mapping_mnn_trajectory.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--trajectory_name',  type="character",              help='') 
p$add_argument('--outfile',          type="character",               help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load mapping functions
source(here::here("rna/mapping/trajectories/mapping_functions.R"))

## START TEST ##
args$query_samples <- opts$samples
args$query_sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
args$query_metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_sce <- file.path(io$atlas.basedir,"results/trajectories/blood_scanpy/blood_SingleCellExperiment.rds")
args$atlas_metadata <- file.path(io$atlas.basedir,"results/trajectories/blood_scanpy/blood_sample_metadata.txt.gz")
args$npcs <- 10
args$n_neighbours <- 15
args$trajectory_name <- "blood"
args$outfile <- paste0(io$basedir,"/results/rna/mapping/trajectories/blood/mapping_mnn.txt.gz")
## END TEST ##

# parse arguments
dir.create(dirname(args$outfile), showWarnings = F) 

# Options
opts$celltype_trajectory_dic <- list(
  # "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3"),
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors", "early_Erythroid", "late_Erythroid"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "endoderm" = c("Epiblast", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  "NMP" = c("Epiblast","Primitive_Streak","Caudal_epiblast","NMP")
)
stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
opts$celltypes <- opts$celltype_trajectory_dic[[args$trajectory_name]]

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% 
  .[pass_rnaQC==TRUE & sample%in%args$query_samples & celltype.mapped%in%opts$celltypes]

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$id_rna, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_query %>% .[cell%in%colnames(sce_query)] %>% setkey(cell) %>% .[colnames(sce_query)]
stopifnot(tmp$cell == colnames(sce_query))
colData(sce_query) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_query),] %>% DataFrame()

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata)# %>% .[pass_rnaQC==TRUE]

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(args$atlas_sce, normalise = TRUE, cells = meta_atlas$cell, remove_non_expressed_genes = TRUE)

# Update colData
tmp <- meta_atlas %>% .[cell%in%colnames(sce_atlas)] %>% setkey(cell) %>% .[colnames(sce_atlas)]
stopifnot(tmp$cell == colnames(sce_atlas))
colData(sce_atlas) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce_atlas),] %>% DataFrame()

#############
## Prepare ## 
#############

if (any(grepl("^ENS",rownames(sce_atlas)))) {
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
}

#####################
## Define gene set ##
#####################

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Filter some genes manually
# genes.intersect <- genes.intersect[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm",genes.intersect,invert=T)]
# genes.intersect <- genes.intersect[!genes.intersect=="Xist"]
# genes.intersect <- genes.intersect[!genes.intersect%in%gene_metadata[chr=="chrY",symbol]]

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#######################
## Feature selection ##
#######################
  
# Calculate mean-variance relationship and extract HVGs
# decomp <- modelGeneVar(sce_atlas, block=sce_atlas$sample)
decomp <- modelGeneVar(sce_atlas)
genes_to_use <- rownames(decomp)[decomp$p.value<=0.01 & decomp$mean>0.1]

#########
## Map ##
#########

meta_query$block <- "query"
meta_atlas$block <- "atlas"

mapping  <- mapWrap(
  sce_atlas = sce_atlas,
  meta_atlas = meta_atlas,
  sce_query = sce_query,
  meta_query = meta_query,
  genes = genes_to_use,
  npcs = args$npcs,
  k = args$n_neighbours
)

##########
## Save ##
##########

mapping.dt <- mapping$mapping %>% 
  .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% 
  # merge(meta_query,by="cell") %>%
  as.data.table

# foo <- merge(meta_query[,c("cell","class","celltype.mapped")],mapping.dt[,c("cell","celltype.mapped")], by="cell", all=T)

fwrite(mapping.dt, args$outfile, sep="\t")

