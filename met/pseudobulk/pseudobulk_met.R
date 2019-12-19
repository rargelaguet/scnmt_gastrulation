library(data.table)
library(purrr)
library(furrr)

######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_cpgs.x), met_cpgs.x := 0L] %>%
    .[is.na(met_cpgs.y), met_cpgs.y := 0L] %>%
    .[is.na(nonmet_cpgs.x), nonmet_cpgs.x := 0L] %>%
    .[is.na(nonmet_cpgs.y), nonmet_cpgs.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_cpgs=met_cpgs.x+met_cpgs.y, nonmet_cpgs=nonmet_cpgs.x+nonmet_cpgs.y)]
}

fread_and_merge <- function(dt, file){
  fread(file, colClasses=list(factor=1L)) %>% 
    setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate")) %>%
    .[,rate:=NULL] %>%
    merge_and_sum(dt)
}

################
## Define I/O ##
################

io <- list()

if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard"
}
io$data <- paste0(io$basedir,"/met/cpg_level")
io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$outdir <- paste0(io$basedir,"/met/cpg_level/pseudobulk")


####################
## Define options ##
####################

opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast",
  
  # E5.5
  "E5.5_Epiblast",
  
  # E6.5
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

# Define which cells to use
opts$cells <- fread(io$metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_metQC==T & stage_lineage%in%opts$stage_lineage,id_met]

# Parallel processing
opts$parallel <- TRUE    # do parallel processing?
opts$ncores <- 2         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_met%in%opts$cells]

##############################
## Load data and pseudobulk ##
##############################

# Parallel processing options
if (opts$parallel){
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}

for (i in opts$stage_lineage) {
  print(i)
  
  # Define input files 
  cells <- sample_metadata[stage_lineage%in%i,id_met]
  
  # cells <- head(cells,n=2)
  files <- paste0(io$data, "/", cells, ".tsv.gz")
  
  # split into chunks for parallel processing
  if (opts$parallel) {
    chunks <- ceiling(seq_along(files)/opts$chunk_size)
    file_list <- split(files, chunks)
  } else {
    file_list <- list(files)
  }
  
  # pseudobulk
  init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_cpgs=as.integer(NA), nonmet_cpgs=as.integer(NA))
  data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
  # data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
    purrr::reduce(merge_and_sum) %>%
    .[,rate:=round(100*met_cpgs/(met_cpgs+nonmet_cpgs))] %>%
    .[,.(chr,pos,met_cpgs,nonmet_cpgs,rate)]
  
  # Save
  outfile = sprintf("%s/%s.tsv",io$outdir,i)
  fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
  system(sprintf("pigz -p %d -f %s",opts$ncores,outfile))
}
