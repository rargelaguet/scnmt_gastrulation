
#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$outdir <- "/Users/ricard/gastrulation/metaccrna/mesendoderm_commitment/endoderm/out"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$outdir <- "/homes/ricard/gastrulation/metaccrna/mesendoderm_commitment/endoderm/out"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")

io$diff.met <- "/Users/ricard/data/gastrulation/met/results/differential"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/results/differential"

# Previously computed pseudotime estimates
io$pseudotime  <- paste0(io$outdir, "/destiny_endoderm.tsv")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
)

opts$acc.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
)

opts$views_names <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm Enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm Enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm Enhancers"
)


# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm"
)


# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature

# Define colors
opts$colors <- c(
  Epiblast="#63B8FF",
  Mesoderm="#CD3278",
  Primitive_Streak="#F4A460",
  Endoderm="#43CD80"
)


# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[stage_lineage%in%opts$stage_lineage]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]

