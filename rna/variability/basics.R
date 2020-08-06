library(SingleCellExperiment)
library(BASiCS)

#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/settings.R")
io$outdir <- paste0(io$basedir,"/rna/results/variability")

## Define options ##
# Define stage and lineage
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak"
  # "E7.5_Epiblast",
  # "E7.5_Ectoderm",
  # "E7.5_Primitive_Streak",  # PLOTS HAVE BEEN DONE WITH THIS UNCOMMENTED
  # "E7.5_Endoderm",
  # "E7.5_Mesoderm"
)

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  # .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage]
  .[pass_rnaQC==T & !is.na(id_met) & stage_lineage%in%opts$stage_lineage]
table(sample_metadata$stage)

###############
## Load data ##
###############

# SingleCellExperiment object
sce <- readRDS(io$rna)[,sample_metadata$id_rna]

############
## basics ##
############

# When technical spike-in genes are not available, BASiCS uses a horizontal integration strategy 
# which borrows information across multiple technical replicates

sce$BatchInfo = sce$plate
Chain <- BASiCS_MCMC(sce, N = 100, Thin = 2, Burn = 2, Regression = FALSE, WithSpikes = FALSE)

basics.object <- SingleCellExperiment(assays = list(counts = counts(sce)))

# The VarThreshold parameter sets a lower threshold for the proportion of variability that is assigned to the biological component 
HVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.6, Plot = FALSE)

# For each gene, these functions return posterior probabilities as a measure of HVG/LVG evidence. A cut-off value for these posterior probabilities is set by controlling the EFDR (as a default option, EFDR is set as 0.10).

format(HVG)



