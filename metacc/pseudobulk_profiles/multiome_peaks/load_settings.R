################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
} else {
  stop()
}

# Genomic contects
io$annos_dir <- io$features.dir

# DNA methylation and chromatin accessibility data
io$met.dir <- io$met_data_raw
io$acc.dir <- io$acc_data_raw

# Output directory
io$outdir <- paste0(io$basedir,"/metacc/exploration/multiome_peaks")

####################
## Define options ##
####################

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E3.5_ICM",
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
  "multiome_peaks" = "Multiome peaks"
)

# Define window positions and characteristics
opts$positions <- c(
  "multiome_peaks"="center"
)

opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150


# Define which cells to use
tmp <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
  .[stage_lineage%in%opts$stage_lineage] 
opts$met.cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T,id_acc]

rm(tmp)
