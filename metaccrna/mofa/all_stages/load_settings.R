library(data.table)
library(purrr)

## Define I/O ##
io <- list()
if (grepl("C02RF23NFVH8",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/gastrulation/metaccrna/mofa/all_stages"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/metaccrna/mofa/all_stages"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$rna.file <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$annos_dir  <- paste0(io$basedir, "/features/genomic_contexts")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)

opts$acc.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)


# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast",
  
  # E5.5
  "E5.5_Epiblast",
  
  # E6.5
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define colors
opts$colors <- c(
  "Epiblast"="grey70",
  "Mesoderm"="#CD3278",
  "Primitive Streak"="sandybrown",
  "Endoderm"="#43CD80",
  "Ectoderm"="steelblue"
)

# window length for the overlap between genes and features
opts$overlapGenes  <- FALSE
opts$gene_window  <- 5e4

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[stage_lineage%in%opts$stage_lineage]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ]
