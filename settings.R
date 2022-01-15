suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  io$multiome.basedir <- "/Users/ricard/data/gastrulation_multiome_10x"
  io$gene.metadata <- io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/Users/ricard/data/mm10_sequence/mm10.genome"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$multiome.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x"
  io$gene.metadata <- io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/hps/nobackup2/research/stegle/users/ricard/mm10_sequence/mm10.genome"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  # io$basedir <- "/Users/argelagr/data/scnmt_gastrulation"
  io$basedir <- "/Users/argelagr/data/scnmt_gastrulation_argelaguet2019"
  io$atlas.basedir <- "/Users/argelagr/data/pijuansala2019_gastrulation10x"
  io$multiome.basedir <- "/Users/argelagr/data/gastrulation_multiome_10x"
  io$gene.metadata <- io$gene_metadata <- "/Users/argelagr/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/Users/argelagr/data/mm10_sequence/mm10.genome"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
    if (grepl("Clark", Sys.info()['effective_user'])) {
      stop()
    } else if (grepl("argelag", Sys.info()['effective_user'])) {
      io$basedir <- "/bi/group/reik/ricard/data/scnmt_gastrulation_argelaguet2019"
      io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
      io$multiome.basedir <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x"
      io$gene.metadata <- io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
      io$mm10.genome <- "/bi/group/reik/ricard/data/mm10_sequence/mm10.genome"
    }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$plate_metadata <- paste0(io$basedir,"/plate_metadata.txt")

# Methylation
io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_pseudobulk_raw <- paste0(io$basedir,"/met/cpg_level/pseudobulk")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met_data_motifs <- paste0(io$basedir,"/met/feature_level/motifs")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$met.stats_per_chr <- paste0(io$basedir,"/met/results/stats/sample_stats_per_chr.txt.gz")

# Accessibility
io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_pseudobulk_raw <- paste0(io$basedir,"/acc/gpc_level/pseudobulk")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc_data_motifs <- paste0(io$basedir,"/acc/feature_level/motifs")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$acc.stats_per_chr <- paste0(io$basedir,"/acc/results/stats/sample_stats_per_chr.txt.gz")

# RNA
io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$rna.stats <- paste0(io$basedir,"/rna/results/stats/rna_stats.txt")
io$cell.cycle <- paste0(io$basedir,"/rna/results/cell_cycle/cell_cycle_scran.txt.gz")

# Other
io$features.dir <- paste0(io$basedir,"/features/genomic_contexts")
io$motifs.dir <- paste0(io$basedir,"/features/motifs")
# io$cpg.density <- paste0(io$basedir,"/met/stats/features/cpg_density_perfeature.txt.gz")
io$scmet <- paste0(io$basedir,"/met/results/variability")
io$mae <- paste0(io$basedir,"/metaccrna/MultiAssayExperiment/scnmtseq_gastrulation_mae.rds")
io$umap <- paste0(io$basedir,"/metaccrna/mofa/all_stages/umap_coordinates.txt")

# RNA atlas (PijuanSala2019)
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")

# Multiome
io$multiome.metadata <- paste0(io$multiome.basedir,"/sample_metadata.txt.gz")
io$multiome.rna.sce <- paste0(io$multiome.basedir,"/processed/rna/SingleCellExperiment.rds")
io$multiome.rna.pseudobulk.sce <- paste0(io$multiome.basedir,"/results/rna/pseudobulk/SingleCellExperiment.rds")
io$multiome.rna.differential <- paste0(io$multiome.basedir,"/results/rna/differential")
io$multiome.atac.peaks.bed <- paste0(io$multiome.basedir,"/PeakCalls/bed/peaks_archR_macs2.bed.gz")
io$multiome.atac.differential.dir <- paste0(io$basedir,"/results/atac/archR/differential/PeakMatrix")
io$multiome.atac.marker_peaks <- paste0(io$multiome.basedir,"/results/atac/archR/differential/PeakMatrix/markers/marker_peaks.txt.gz")
io$multiome.atac.peak_metadata <- paste0(io$multiome.basedir,"/PeakCalls/peak_metadata.tsv.gz")
io$multiome.atac.peak_stats <- paste0(io$basedir,"/results/atac/archR/peak_calling/peak_stats.txt.gz")
io$multiome.atac.pseudobulk.peakMatrix.se <- paste0(io$multiome.basedir,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")

#############
## Options ##
#############

opts <- list()

# TO-DO: ADD E8.5 CELL TYPES
# opts$stage_lineage <- c(
#   "E3.5_ICM",
#   "E4.5_Epiblast",
#   "E4.5_Primitive_endoderm",
#   "E5.5_Epiblast",
#   "E5.5_Visceral_endoderm",
#   "E6.5_Epiblast",
#   "E6.5_Primitive_Streak",
#   "E6.5_Visceral_endoderm",
#   "E7.5_Epiblast",
#   "E7.5_Ectoderm",
#   "E7.5_Primitive_Streak",
#   "E7.5_Endoderm",
#   "E7.5_Mesoderm"
# )

# TO-DO: ADD E8.5 CELL TYPES
# opts$stagelineage.colors <- c(
#   "E3.5_ICM" = "#b63fba",
#   "E4.5_Epiblast" = "#C1CDCD",
#   "E4.5_Primitive_endoderm" = "darkgreen",
#   "E5.5_Epiblast" = "#C1CDCD",
#   "E5.5_Visceral_endoderm" = "darkgreen",
#   "E6.5_Epiblast" = "#C1CDCD",
#   "E6.5_Visceral_endoderm" = "darkgreen",
#   "E6.5_Primitive_Streak"="sandybrown",
#   "E6.5_ExE_ectoderm"="black",
#   "E7.5_Epiblast" = "#C1CDCD",
#   "E7.5_Primitive_Streak"="sandybrown",
#   "E7.5_Ectoderm" = "steelblue",
#   "E7.5_Endoderm" = "#43CD80",
#   "E7.5_Mesoderm" = "#CD3278"
# )

opts$celltypes = c(
  "ICM",
  "Primitive_endoderm",
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
  "ICM" = "#C6E2FF",
  "Primitive_endoderm" = "darkgreen",
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

opts$celltype2.colors = c(
  "ICM" = "#989898",
  "Primitive_endoderm" = "darkgreen",
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors" = "#c9a997",
  "Erythroid" = "#EF4E22",
  "NMP" = "#8EC792",
  "Neurectoderm" = "#65A83E",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#1A1A1A"
)


opts$celltype3.colors = c(
  "ICM" = "#b63fba",
  "Primitive_endoderm" = "#BC8F8F",
  "Epiblast" = "grey70",
  "Primitive_Streak" = "sandybrown",
  "PGC" = "#FACB12",
  "Endoderm" = "#43CD80",
  "Mesoderm" = "#CD3278",
  "Ectoderm" = "steelblue",
  "Erythroid" = "#EF4E22",
  "NMP" = "#8EC792",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#1A1A1A"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))
opts$stages <- c("E3.5", "E4.5", "E5.5", "E6.5", "E7.5", "E8.5")

opts$stage.colors <- viridis::viridis(n=length(opts$stages))
names(opts$stage.colors) <- rev(opts$stages)

##########################
## Load sample metadata ##
##########################

# factor.cols <- c("id_rna","id_met","id_acc","stage","lineage","lab","plate","embryo")

# sample_metadata <- fread(io$metadata) %>% 
#   .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]
  # %>% .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>% droplevels
