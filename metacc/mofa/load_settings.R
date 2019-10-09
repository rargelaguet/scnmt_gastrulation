library(data.table)
library(purrr)

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/gastrulation_norsync_stuff/mofa"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  # io$outdir <- "/homes/ricard/gastrulation/metaccrna/mofa/Primitive_Streak/out"
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata10x.txt")
io$met.dir <- paste0(io$basedir,"/met/parsed")
io$acc.dir <- paste0(io$basedir,"/acc/parsed")
io$annos_dir  <- paste0(io$basedir, "/features/filt")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  # "genebody",
  # "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)

opts$acc.annos <- c(
  # "genebody",
  # "prom_200_200",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)


# Define which stage and lineages to look at 
# opts$stage_lineage <- c(
#   # "E6.5_EPI", 
#   # "E6.5_PS",
#   # "E6.5_NOIDEA"
#   "E7.5_Ectoderm",
#   "E7.5_Mesoderm",
#   "E7.5_Endoderm"
# )

opts$stage_lineage10x <- c(
  # "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm",

  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Nascent_mesoderm",
  "E6.5_Anterior_Primitive_Streak",
  # "E6.5_Visceral_endoderm",
  # "E6.5_ExE_endoderm",
  
  "E7.5_Notochord",
  "E7.5_Embryonic_endoderm",
  # "E7.5_Visceral_endoderm",
  # "E7.5_ExE_endoderm",
  
  "E7.5_Nascent_mesoderm",
  "E7.5_Mature_mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Anterior_Primitive_Streak",
  
  "E7.5_Epiblast"
  # "E7.5_ExE_ectoderm",
  # "E7.5_Surface_ectoderm"
  # "E7.5_Rostral_neurectoderm"
  # "E7.5_Caudal_epiblast"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
# opts$met_min.cells <- 5      # minimum number of cells per feature (and per stage_lineage)
opts$met_min.cells <- 50      # minimum number of cells per feature
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 3        # minimum number of GpC sites per feature
# opts$acc_min.cells <- 5      # minimum number of cells per feature (and per stage_lineage)
opts$acc_min.cells <- 50      # minimum number of cells per feature
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Use only samples that passed QC for all omics (opts$include_all=FALSE)?
opts$include_all <- TRUE

# window length for the overlap between genes and features
opts$overlapGenes  <- FALSE
opts$gene_window  <- 5e4

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  # .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  # .[stage_lineage%in%opts$stage_lineage]
  .[stage_lineage%in%opts$stage_lineage10x]
  # .[!is.na(id_met) & !is.na(id_acc)]

if (opts$include_all==T) { 
  opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
  opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]
} else {
  tmp <- tmp %>% .[pass_metQC==T & pass_accQC==T]
  opts$met_cells <- tmp %>% .[,id_met]
  opts$acc_cells <- tmp %>% .[,id_acc]
}

sample_metadata <- fread(io$sample.metadata,stringsAsFactors=T) %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  # .[,stage_lineage10x:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[id_met%in%opts$met_cells | id_acc %in% opts$acc_cells ] %>%
  droplevels()


# sample_metadata[,sum(pass_metQC,na.rm=T),by=c("stage_lineage")]
