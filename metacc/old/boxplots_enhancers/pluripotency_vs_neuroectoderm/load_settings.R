#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
}

# Sample metadata file
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

# Data folders
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")

# File with global statistics per cell
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")

# Folders with the differential analysis results
io$diff.met <- paste0(io$basedir,"/met/results/differential")
io$diff.acc <- paste0(io$basedir,"/acc/results/differential")

# Output directory
io$outdir <- paste0(io$basedir,"/metacc/boxplots_neuroectoderm_enhancers")

## Define options ##
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E7.5_Ectoderm"
  
)

# Define genomic contexts for methylation
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  # "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  # "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers"
)

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[stage_lineage%in%opts$stage_lineage] 
opts$met_cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc_cells <- tmp %>% .[pass_accQC==T,id_acc]