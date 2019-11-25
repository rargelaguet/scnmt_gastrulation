#####################
## Define settings ##
#####################

## Define I/O ##

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  stop()
}
io$outdir <- paste0(io$basedir,"/metaccrna/mesendoderm_commitment/density_maps")

# Folders with the differential analysis results
io$diff.met <- paste0(io$basedir,"/met/results/differential")
io$diff.acc <- paste0(io$basedir,"/acc/results/differential")

## Define options ##
opts <- list()

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

opts$annos <- opts$met.annos <- opts$acc.annos

# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2

opts$min.fdr <- 0.10
opts$min.acc.diff <- 5
opts$min.met.diff <- 5