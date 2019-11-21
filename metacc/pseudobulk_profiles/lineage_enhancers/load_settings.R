
################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  stop()
}

# Sample metadata
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")

# Genomic contects
io$annos_dir <- paste0(io$basedir,"/features/genomic_contexts")

# DNA methylation and chromatin accessibility data
io$met.dir <- paste0(io$basedir,"/met/cpg_level")
io$acc.dir <- paste0(io$basedir,"/acc/gpc_level")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")

# Output directory
io$pdfdir <- paste0(io$basedir,"/metacc/pseudobulked_profiles/profiles_lineage_enhancers")

####################
## Define options ##
####################

opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Endoderm"
  
)

# Define genomic contexts
opts$annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)

# Define window positions and characteristics
opts$positions <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="center",
  "H3K27ac_distal_E7.5_End_intersect12"="center",
  "H3K27ac_distal_E7.5_Ect_intersect12"="center"
)
opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[stage_lineage%in%opts$stage_lineage] 
opts$met.cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T,id_acc]

rm(tmp)