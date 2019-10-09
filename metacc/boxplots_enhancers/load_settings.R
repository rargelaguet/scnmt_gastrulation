#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$met.dir <- "/Users/ricard/data/gastrulation/met/parsed"
io$acc.dir <- "/Users/ricard/data/gastrulation/acc/parsed"
io$outdir <- "/Users/ricard/gastrulation/metaccrna/boxplots_enhancers/out"

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/gastrulation/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/differential/feature_level"

# Global statistics per cell
io$met.stats <- "/Users/ricard/data/gastrulation/met/stats/samples/sample_stats.txt"
io$acc.stats <- "/Users/ricard/data/gastrulation/acc/stats/samples/sample_stats.txt"

## Define options ##
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm"
)

# Define genomic contexts for methylation
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
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

opts$colors <- c(
  "E4.5 Epiblast"="#C1CDCD",
  "E5.5 Epiblast"="#C1CDCD",
  "E6.5 Epiblast"="#C1CDCD",
  "E7.5 Epiblast"="#C1CDCD",
  "E7.5 Ectoderm"="steelblue",
  "E6.5 Primitive Streak"="sandybrown",
  "E7.5 Primitive Streak"="sandybrown",
  "E7.5 Endoderm"="#43CD80",
  "E7.5 Mesoderm"="#CD3278"
)

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[stage_lineage%in%opts$stage_lineage] 
opts$met_cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc_cells <- tmp %>% .[pass_accQC==T,id_acc]