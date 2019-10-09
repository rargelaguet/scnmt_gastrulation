
#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$met.dir <- "/Users/ricard/data/gastrulation/met/raw"
io$acc.dir <- "/Users/ricard/data/gastrulation/acc/raw"
io$outdir <- "/Users/ricard/gastrulation/metaccrna/profiles_enhancers/out"
io$annos_dir <- "/Users/ricard/data/gastrulation/features/filt"
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"

# Folders with the global statistics per cell
io$met.stats <- "/Users/ricard/data/gastrulation/met/stats/samples/sample_stats.txt"
io$acc.stats <- "/Users/ricard/data/gastrulation/acc/stats/samples/sample_stats.txt"

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/gastrulation/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/differential/feature_level"

## Define options ##
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  "E4.5_Epiblast",
  # "E4.5_Primitive_endoderm",
  
  "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm",

  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  # "E6.5_Nascent_mesoderm",
  # "E6.5_Visceral_endoderm",

  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm"
  
)
# Define genomic contexts for methylation
opts$annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
)

# Define window positions
opts$positions <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="center",
  "H3K27ac_distal_E7.5_Mes_intersect12"="center",
  "H3K27ac_distal_E7.5_End_intersect12"="center"
)

opts$overlapGenes <- FALSE

opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.met.diff <- 5
opts$min.acc.diff <- 5

opts$window_size <- 2000
# opts$met.tile <- 150
# opts$acc.tile <- 100

opts$met.tile <- 200
opts$acc.tile <- 150


# opts$colors <- c(
#   E4.5_Epiblast="#C1CDCD",
#   # E5.5_Epiblast="#9AC0CD",
#   E6.5_Epiblast="#C1CDCD",
#   E7.5_Epiblast="#C1CDCD",
#   E7.5_Ectoderm="steelblue",
#   E6.5_Primitive_Streak="sandybrown",
#   # E7.5_Primitive_Streak="#FF7F24",
#   E7.5_Nascent_mesoderm="#FF82AB",
#   E7.5_Mature_mesoderm="#CD3278",
#   E7.5_Embryonic_endoderm="#43CD80"
#   # E7.5_Visceral_endoderm="#2E8B57"
# )

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] 

if (opts$stage_lineage[1]=="all") {
  stop()
} else {
  tmp <- tmp %>% .[stage_lineage%in%opts$stage_lineage] 
}
opts$met.cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T,id_acc]