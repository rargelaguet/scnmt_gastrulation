#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/utils.R")
} else {
  stop()
}

# I/O
io$umap <- paste0(io$basedir,"/metaccrna/mofa/all_stages/umap_coordinates.txt")
io$outdir <- paste0(io$basedir,"/metaccrna/mefisto")
io$mofa.outfile <- paste0(io$outdir,"/mefisto_model_motifs.rds")

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast",
  # "E4.5_Primitive_endoderm",
  
  # E5.5
  "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm",
  
  # E6.5
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  # "E6.5_Visceral_endoderm",
  # "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
  # "E7.5_Visceral_endoderm"
)

# Filtering options for methylation
opts$met_min.cells <- 50      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 250    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.cells <- 50      # minimum number of cells per feature (per stage)
opts$acc_nfeatures <- 250    # maximum number of features per view (filter based on variance)

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage] %>%
  droplevels

opts$met_cells <- sample_metadata %>% .[pass_metQC==T,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T,id_rna]
opts$acc_cells <- sample_metadata %>% .[pass_accQC==T,id_acc]
