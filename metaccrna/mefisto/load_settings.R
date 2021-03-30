#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else {
  stop()
}

# I/O
io$outdir <- paste0(io$basedir,"/metaccrna/mefisto")
io$mofa.outfile <- paste0(io$outdir,"/mefisto_model.rds")

## Define options ##
# Define which annotations to look at
opts$met.annos <- c(
  "prom_2000_2000",
  # "genebody",
  # "E3.5_H3K27ac_distal",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "E10.5_midbrain_H3K27ac_distal",
  # "E10.5_heart_H3K27ac_distal",
  # "E12.5_intestine_H3K27ac_distal"
)

opts$acc.annos <- c(
  "prom_2000_2000",
  # "genebody",
  # "E3.5_H3K27ac_distal",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "E10.5_midbrain_H3K27ac_distal",
  # "E10.5_heart_H3K27ac_distal",
  # "E12.5_intestine_H3K27ac_distal"
)


opts$rename.annos <- c(
  "prom_2000_2000"="Promoters",
  "prom_200_200"="Promoters",
  "genebody"="Gene bodies",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Enhancers"
)


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
  # "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  # "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
  # "E7.5_Visceral_endoderm"
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 1500    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 25      # minimum number of cells per feature (per stage)
opts$acc_nfeatures <- 1500    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
# opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
# opts$rna_ngenes <- 5000       # maximum number of genes (filter based on variance)

# Deefine cell type colors
# opts$colors <- c(
#   "Epiblast"="grey70",
#   "Mesoderm"="#CD3278",
#   "Primitive Streak"="sandybrown",
#   "Endoderm"="#43CD80",
#   "Ectoderm"="steelblue",
#   "ExE Endoderm"="#E066FF"
# )

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>% 
  # .[stage=="E7.5" & lineage10x=="Visceral_endoderm",lineage10x_2:="Visceral_endoderm"] %>%
  # .[lineage10x_2=="Visceral endoderm",lineage10x_2:="ExE Endoderm"] %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage] %>%
  # .[pass_rnaQC==T & (pass_accQC==T | pass_metQC==T) & stage_lineage%in%opts$stage_lineage] %>%
  # .[,c("id_rna","stage","lineage10x_2","stage_lineage")] %>%
  droplevels

table(sample_metadata$stage_lineage)

opts$met_cells <- sample_metadata %>% .[pass_metQC==T,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T,id_rna]
opts$acc_cells <- sample_metadata %>% .[pass_accQC==T,id_acc]
