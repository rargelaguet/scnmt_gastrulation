#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$sample.metadata <- "/Users/ricard/data/gastrulation/public_data/Dai_2016/sample_metadata.txt"
io$met.dir <- "/Users/ricard/data/gastrulation/public_data/Dai_2016/met/feature_level"
io$outdir <- "/Users/ricard/gastrulation/public_datasets/Dai2016/boxplots/out"

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/gastrulation/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/differential/feature_level"

## Define options ##
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E6.5_Epiblast"
)

# Define genomic contexts for methylation
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)


# How to select lineage-defining hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.acc.diff <- 5
opts$min.met.diff <- 5

# opts$colors <- c(
#   "TKO" = "blue",
#   "WT" = "red"
# )

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>% 
  .[stage_lineage%in%opts$stage_lineage] 
opts$samples <- tmp %>% .[,id_met]

