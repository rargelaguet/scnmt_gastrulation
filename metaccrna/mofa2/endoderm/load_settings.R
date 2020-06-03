suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(MOFA2))

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

################
## Define I/O ##
################

source("/Users/ricard/scnmt_gastrulation/settings.R")

####################
## Define options ##
####################

# Define which annotations to look at
opts$met.annos <- c(
  # "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
)

opts$acc.annos <- c(
  # "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
)


opts$rename.annos <- c(
  "prom_2000_2000"="Promoters",
  "prom_200_200"="Promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Enhancers"
)


# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 10      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 10      # minimum number of cells per feature (per stage)
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 1000       # maximum number of genes (filter based on variance)

# Deefine cell type colors
opts$colors <- c(
  "Epiblast"="grey70",
  "Mesoderm"="#CD3278",
  "Primitive Streak"="sandybrown",
  "Endoderm"="#43CD80",
  "Ectoderm"="steelblue",
  "ExE Endoderm"="#E066FF"
)

# Define which cells to use
sample_metadata <- sample_metadata %>%
  .[pass_rnaQC==TRUE & !is.na(id_met)] %>%
  .[,lineage10x_2:=stringr::str_replace_all(lineage10x_2,"_"," ")] %>%
  .[stage_lineage%in%c("E7.5_Endoderm")] %>%
  .[lineage10x%in%c("Def._endoderm","Gut","Visceral_endoderm","Anterior_Primitive_Streak")]

table(sample_metadata$lineage10x)
table(sample_metadata$pass_metQC)
table(sample_metadata$pass_accQC)
