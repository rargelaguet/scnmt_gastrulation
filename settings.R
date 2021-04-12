suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/Users/ricard/data/mm10_sequence/mm10.genome"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$gene.metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/hps/nobackup2/research/stegle/users/ricard/mm10_sequence/mm10.genome"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

# Methylation
io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met_data_motifs <- paste0(io$basedir,"/met/feature_level/motifs")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$met.stats_per_chr <- paste0(io$basedir,"/met/results/stats/sample_stats_per_chr.txt.gz")

# Accessibility
io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc_data_motifs <- paste0(io$basedir,"/acc/feature_level/motifs")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$acc.stats_per_chr <- paste0(io$basedir,"/acc/results/stats/sample_stats_per_chr.txt.gz")

# RNA
io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$rna.stats <- paste0(io$basedir,"/rna/results/stats/rna_stats.txt")

# Other
io$features.dir <- paste0(io$basedir,"/features/genomic_contexts")
io$motifs.dir <- paste0(io$basedir,"/features/motifs")
# io$cpg.density <- paste0(io$basedir,"/met/stats/features/cpg_density_perfeature.txt.gz")
io$scmet <- paste0(io$basedir,"/met/results/variability")
io$mae <- paste0(io$basedir,"/metaccrna/MultiAssayExperiment/scnmtseq_gastrulation_mae.rds")

# RNA atlas (PijuanSala2019)
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Visceral_endoderm",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

opts$celltype2.colors <- c(
  "Epiblast"="grey70",
  "Mesoderm"="#CD3278",
  "Primitive_Streak"="sandybrown",
  "Endoderm"="#43CD80",
  "Ectoderm"="steelblue",
  "Epiblast/Ectoderm"="steelblue",
  "Visceral_endoderm"="darkgreen"
)

opts$stagelineage.colors <- c(
  "E4.5_Epiblast" = "#C1CDCD",
  "E4.5_Primitive_endoderm" = "darkgreen",
  "E5.5_Epiblast" = "#C1CDCD",
  "E5.5_Visceral_endoderm" = "darkgreen",
  "E6.5_Epiblast" = "#C1CDCD",
  "E6.5_Visceral_endoderm" = "darkgreen",
  "E6.5_Primitive_Streak"="sandybrown",
  "E7.5_Epiblast" = "#C1CDCD",
  "E7.5_Primitive_Streak"="sandybrown",
  "E7.5_Ectoderm" = "steelblue",
  "E7.5_Endoderm" = "#43CD80",
  "E7.5_Mesoderm" = "#CD3278"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$celltype.colors = c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))
opts$stages <- c("E4.5", "E5.5", "E6.5", "E7.5")

##########################
## Load sample metadata ##
##########################

# factor.cols <- c("id_rna","id_met","id_acc","stage","lineage","lab","plate","embryo")

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]
  # %>% .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>% droplevels
