library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
}
io$data.dir <- paste0(io$basedir,"/met/parsed")
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$outdir <- paste0(io$basedir,"/met/dimensionality_reduction")

####################
## Define options ##
####################

opts <- list()

# Define which annotations to look at
opts$annos <- c(
  "ESC_DHS"="DHS"
)

# Define which cells to use
opts$stage_lineage <- c(

  # E4.5
  # "E4.5_Epiblast",
  # "E4.5_Primitive_endoderm"

  # E5.5
  # "E5.5_Epiblast",
  # "E5.5_Visceral_endoderm"

  # E6.5
  # "E6.5_Epiblast",
  # "E6.5_Primitive_Streak",
  # "E6.5_Mesoderm",
  # "E6.5_Visceral_endoderm"
  
  # E7.5
  "E7.5_Endoderm",
  "E7.5_Mesoderm",
  "E7.5_Ectoderm"
)

# Define which cells to  use
opts$cells <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_metQC==T & stage_lineage%in%opts$stage_lineage,id_met]

# Stage colors
# opts$colors <- c(
#   "E4.5"="#FDCC8A", 
#   "E5.5"="#FC8D59", 
#   "E6.5"="#E34A33", 
#   "E7.5"="#600707"
# )

# Lineage colors
opts$colors <- c(
  ExE="#fc8d62", 
  Embryonic="#8da0cb", 
  Epiblast="#63B8FF",
  Ectoderm="steelblue",
  Mesoderm="#CD3278",
  Primitive_Streak="sandybrown",
  Endoderm="#43CD80"
)


# Filtering options
opts$min.CpGs <- 1          # minimum number of CpG sites per feature and cell
opts$min.coverage <- 0.10   # minimum coverage (fraction of cells with at least min.CpG measurements)
opts$nfeatures <- 5000     # number of features per view (filter based on variance)
# opts$nfeatures <- 10000     # number of features per view (filter based on variance)

# Output file
io$outfile = sprintf("%s/hdf5/model_%s_%s.hdf5",io$outdir,paste(names(opts$annos), collapse="_"), paste(opts$stage_lineage, collapse="_"))


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata, stringsAsFactors=T, showProgress=F) %>% 
  .[id_met%in%opts$cells] %>%
  .[,c("id_met","stage","lineage10x_2")] %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>% droplevels()
