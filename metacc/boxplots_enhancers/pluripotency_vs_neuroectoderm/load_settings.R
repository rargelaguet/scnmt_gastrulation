#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$met.dir <- "/Users/ricard/data/gastrulation/met/parsed"
io$acc.dir <- "/Users/ricard/data/gastrulation/acc/parsed"
io$outdir <- "/Users/ricard/data/gastrulation/metacc/boxplots_neuroectoderm_enhancers"
io$annos_dir <- "/Users/ricard/data/gastrulation/features/filt"
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/gastrulation/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/differential/feature_level"

# Folders with the global statistics per cell
io$met.stats <- "/Users/ricard/data/gastrulation/met/stats/samples/sample_stats.txt"
io$acc.stats <- "/Users/ricard/data/gastrulation/acc/stats/samples/sample_stats.txt"

## Define options ##
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  # "E7.5_Epiblast",
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

opts$colors <- c(
  "E4.5 Epiblast"="grey70",
  "E5.5 Epiblast"="grey70",
  "E6.5 Epiblast"="grey70",
  "E7.5 Epiblast"="grey70",
  "E7.5 Ectoderm"="steelblue"
)

opts$width = 10
opts$height = 5

# Define which cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[stage_lineage%in%opts$stage_lineage] 
opts$met_cells <- tmp %>% .[pass_metQC==T,id_met]
opts$acc_cells <- tmp %>% .[pass_accQC==T,id_acc]