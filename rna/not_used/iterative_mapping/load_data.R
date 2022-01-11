
################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")]

if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce.atlas  <- readRDS(io$atlas.sce)[,meta_atlas$cell]

# Add metadata to the SCE object
colData(sce.atlas) <- meta_atlas %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(io$metadata) %>% 
  .[pass_rnaQC==T & stage%in%c("E6.5","E7.5")] %>%
  setnames("id_rna","cell")

# Load SingleCellExperiment
sce.query <- readRDS(io$rna)[,meta_query$cell]

# Add metadata to the SCE object
meta_query.df <- meta_query %>% DataFrame
rownames(meta_query.df) <- meta_query.df$cell
colData(sce.query) <- meta_query.df

#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce.query <- sce.query[rowMeans(counts(sce.query))>1e-5,]
sce.atlas <- sce.atlas[rowMeans(counts(sce.atlas))>1e-5,]

# Intersect genes
genes.intersect <- intersect(rownames(sce.query), rownames(sce.atlas))
sce.query  <- sce.query[genes.intersect,]
sce.atlas <- sce.atlas[genes.intersect,]

# Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$atlas.marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# marker_genes <- marker_genes[marker_genes%in%genes.intersect]
# 
# stopifnot(all(marker_genes%in%rownames(sce.atlas)))
# stopifnot(all(marker_genes%in%rownames(sce.query)))
# 
# sce.query <- sce.query[marker_genes,]
# sce.atlas <- sce.atlas[marker_genes,]
