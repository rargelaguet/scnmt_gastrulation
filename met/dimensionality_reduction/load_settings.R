io <- list()
opts <- list()

## Define I/O ##
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
  io$data.dir <- "/Users/ricard/data/gastrulation/met/parsed"
  io$outdir <- "/Users/ricard/data/gastrulation_norsync_stuff/met/dimensionality_reduction"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$sample.metadata <- "/hps/nobackup/stegle/users/ricard/gastrulation/sample_metadata.txt"
  io$data.dir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/parsed"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/dimensionality_reduction"
}
io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")

## Define options ##

# Define which annotations to look at
opts$annos <- c(
  # "genebody"="Gene body"
  # "prom_2000_2000"="Promoters"
  "ESC_DHS"="DHS"
  # "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  # "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  # "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers"
  # "H3K4me3_E7.5_Mes" = "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End" = "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect" = "H3K4me3_E7.5_Ect"
  # "LINE"="LINE",
  # "LTR"="LTR"
  # "window2000_step1000"="window2000_step1000"
)

# Define which cells to use
opts$stage_lineage10x <- c(

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

opts$cells <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_metQC==T & stage_lineage%in%opts$stage_lineage,id_met]

# opts$colors <- c(
#   "E4.5"="#FDCC8A", 
#   "E5.5"="#FC8D59", 
#   "E6.5"="#E34A33", 
#   "E7.5"="#600707"
# )

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
io$outfile = sprintf("%s/hdf5/model_%s_%s.hdf5",io$outdir,paste(names(opts$annos), collapse="_"), paste(opts$stage_lineage10x, collapse="_"))


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$sample.metadata, stringsAsFactors=T, showProgress=F) %>% 
  .[id_met%in%opts$cells] %>%
  .[,c("id_met","stage","lineage10x_2")] %>%
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>% droplevels()
