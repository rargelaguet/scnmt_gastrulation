library(scater)
library(RColorBrewer)
library(scran)

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
  "E6.5_Primitive_Streak",
  # "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
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

####################################
## Calculate variability per gene ##
####################################

# Gaussian variance
var <- apply(logcounts(sce),1,var)
dt.var <- data.table(ens_id=names(var), var=var)

# Biological overdispersion using scran
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
dt.biodisp <- data.table(ens_id=rownames(decomp), bio.disp=decomp$bio)

##########
## Save ##
##########

dt <- merge(dt.var,dt.biodisp, by="ens_id")

# Save results
fwrite(dt, paste0(io$outdir,"/gene_variability.txt.gz"), sep="\t", quote=F)
