library(muscat)
library(DESeq2)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/rna")

# Options
# opts$stage_lineage <- c(
#   "E3.5_ICM",
#   "E4.5_Epiblast",
#   "E5.5_Epiblast",
#   "E6.5_Epiblast",
#   "E6.5_Primitive_Streak",
#   "E7.5_Ectoderm",
#   "E7.5_Mesoderm",
#   "E7.5_Endoderm"
# )

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(io$metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[pass_rnaQC==TRUE & stage_lineage%in%opts$stage_lineage]
table(sample_metadata$stage_lineage)

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$id_rna)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("id_rna") %>% DataFrame

###################################
## Aggregate counts per celltype ##
###################################

# assays(sce)$cpm <- edgeR::cpm(assay(sce), normalized.lib.sizes = FALSE, log = FALSE)

sce_pseudobulk <- aggregateData(
  sce,
  assay = "counts",
  by = c("stage_lineage"),
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###############
## Normalise ##
###############

# create DESeq object
dds <- DESeqDataSet(sce_pseudobulk, design=~1)

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) 
# and then transforms the count data (normalized by division by the size factors or normalization factors), 
# yielding a matrix of values which are now approximately homoskedastic 
dds <- varianceStabilizingTransformation(dds)

logcounts(sce_pseudobulk) <- assay(dds)

###################
## Sanity checks ##
###################

# cor(
#   colMeans(logcounts(sce_pseudobulk)),
#   metadata(sce_pseudobulk)$n_cells
# )

##########
## Save ##
##########

saveRDS(sce_pseudobulk, paste0(io$outdir,"/SingleCellExperiment_pseudobulk.rds"))
