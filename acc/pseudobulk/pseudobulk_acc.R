library(furrr)


################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
  source("/Users/ricard/scnmt_gastrulation/acc/pseudobulk/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/utils.R")
  source("/homes/ricard/scnmt_gastrulation/acc/pseudobulk/utils.R")
} else {
  stop()
}

io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$outdir <- paste0(io$basedir,"/acc/gpc_level/pseudobulk")


####################
## Define options ##
####################

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  "E3.5_ICM",  
  # "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  # "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Visceral_endoderm",
  "E6.5_Mesoderm",
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

# Define which cells to use
opts$cells <- fread(io$metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_accQC==TRUE & stage_lineage%in%opts$stage_lineage,id_acc]

# Parallel processing
opts$parallel <- FALSE    # do parallel processing?
opts$ncores <- 1         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_acc%in%opts$cells]

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
  outfile = sprintf("%s/%s.tsv.gz",io$outdir,i)
  if (file.exists(outfile)) {
    print(sprintf("%s already exists, skipping...",outfile))
  } else {
    
    # Define input files 
    cells <- sample_metadata[stage_lineage%in%i,id_acc]
    
    # cells <- head(cells,n=2)
    files <- paste0(io$acc_data_raw, "/", cells, ".tsv.gz")
    
    # split into chunks for parallel processing
    if (opts$parallel) {
      chunks <- ceiling(seq_along(files)/opts$chunk_size)
      file_list <- split(files, chunks)
    } else {
      file_list <- list(files)
    }
    
    # pseudobulk
    init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_cpgs=as.integer(NA), nonmet_cpgs=as.integer(NA))
    # data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
    data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
      purrr::reduce(merge_and_sum) %>%
      .[,rate:=round(100*met_cpgs/(met_cpgs+nonmet_cpgs))] %>%
      .[,.(chr,pos,met_cpgs,nonmet_cpgs,rate)]
    
    # filter
    data <- data %>% .[complete.cases(.)]
    
    # Save
    fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
  }
}
