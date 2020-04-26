suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene.metadata <- "/Users/ricard/data/ensembl/"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation"
  io$gene.metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met.stats <- paste0(io$basedir,"/met/stats/sample_stats.txt")

io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc.stats <- paste0(io$basedir,"/acc/stats/sample_stats.txt")

io$rna <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

io$features.dir <- paste0(io$basedir,"/features/genomic_contexts")
# io$cpg.density <- paste0(io$basedir,"/met/stats/features/cpg_density_perfeature.txt.gz")

io$scmet <- paste0(io$basedir,"/met/results/variability")

#############
## Options ##
#############

opts <- list()

# opts$stage_lineage <- c(

#   # E4.5
#   "E4.5_Epiblast",

#   # E5.5
#   "E5.5_Epiblast",
  
#   # E6.5
#   "E6.5_Epiblast",
#   "E6.5_Primitive_Streak",
  
#   # E7.5
#   "E7.5_Epiblast",
#   "E7.5_Ectoderm",
#   "E7.5_Primitive_Streak",
#   "E7.5_Endoderm",
#   "E7.5_Mesoderm"
# )

opts$colors_lineages <- c(
  "Epiblast"="grey70",
  "Mesoderm"="#CD3278",
  "Primitive_Streak"="sandybrown",
  "Endoderm"="#43CD80",
  "Ectoderm"="steelblue",
  "Epiblast/Ectoderm"="steelblue"
)

opts$colors_stages <- c(
  "E6.5"="grey70",
  "E7.5"="#CD3278"
)

##########################
## Load sample metadata ##
##########################

# factor.cols <- c("id_rna","id_met","id_acc","stage","lineage","lab","plate","embryo")

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]
  # %>% .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>% droplevels