library(furrr)

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
}

io$outdir <- paste0(io$basedir,"/met/bedgraph")


####################
## Define options ##
####################

# Define groups based on stage and lineage
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

# Parallel processing
opts$parallel <- TRUE    # do parallel processing?
opts$ncores <- 2         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[pass_metQC==TRUE & stage_lineage%in%opts$stage_lineage,id_met]


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
  
  # filter
  data <- data[complete.cases(.)]
  
  # Save
  outfile = sprintf("%s/%s.bgr",io$outdir,i)
  fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
  system(sprintf("pigz -p %d -f %s",opts$ncores,outfile))
}
