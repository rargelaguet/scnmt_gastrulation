
#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/gastrulation/metaccrna/mesendoderm_commitment/mesoderm/out"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/homes/ricard/gastrulation/metaccrna/mesendoderm_commitment/mesoderm/out"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/parsed")
io$acc.dir <- paste0(io$basedir,"/acc/parsed")
io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/samples/sample_stats.txt")
io$rna.file <- paste0(io$basedir,"/rna/parsed/SingleCellExperiment.rds")
io$annos_dir  <- paste0(io$basedir, "/features/filt")

# Previously computed pseudotime estimates
io$pseudotime  <- paste0(io$outdir, "/destiny_mesoderm.tsv")

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


# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Mesoderm"
)


# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define colors
opts$colors <- c(
  Epiblast="#63B8FF",
  Mesoderm="#CD3278",
  Primitive_Streak="#F4A460",
  Endoderm="#43CD80"
)

opts$views_names <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm Enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm Enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm Enhancers"
)

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[stage_lineage%in%opts$stage_lineage]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]


