source("/Users/ricard/scnmt_gastrulation/settings.R")

# Define which annotations to use
opts$met.annos <- c(
  # "prom_2000_2000",
  # "genebody",
  # "E3.5_H3K27ac_distal",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "E10.5_midbrain_H3K27ac_distal",
  # "E10.5_heart_H3K27ac_distal",
  # "E12.5_intestine_H3K27ac_distal"
)



# Define which stage and lineages to look at 
opts$lineages <- c(

  # E5.5
  # "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm",
  
  # E6.5
  # "E6.5_Epiblast",
  # "E6.5_Primitive_Streak",
  # "E6.5_Visceral_endoderm",
  # "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm",
  "E7.5_Visceral_endoderm"
)
# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 25      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 2000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2000       # maximum number of genes (filter based on variance)

############################
## Update sample metadata ##
############################

# Define which cells to use
sample_metadata <- sample_metadata %>%
  .[stage=="E7.5" & lineage10x%in%c("Visceral_endoderm","ExE_endoderm"),lineage10x_2:="Visceral_endoderm"] %>%
  .[pass_rnaQC==TRUE & pass_metQC==TRUE] %>%
  .[stage_lineage%in%opts$lineages]
opts$met_cells <- sample_metadata %>% .[pass_metQC==T,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T,id_rna]

# table(sample_metadata$lineage10x)
table(sample_metadata$stage_lineage)
