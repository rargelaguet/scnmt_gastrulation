library(data.table)
library(purrr)
library(furrr)

######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  walk(list(dt1, dt2), setkey, chr, pos)
  merge(dt1, dt2, all = TRUE) %>%
    .[is.na(met_reads.x), met_reads.x := 0L] %>%
    .[is.na(met_reads.y), met_reads.y := 0L] %>%
    .[is.na(nonmet_reads.x), nonmet_reads.x := 0L] %>%
    .[is.na(nonmet_reads.y), nonmet_reads.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_reads=met_reads.x+met_reads.y, nonmet_reads=nonmet_reads.x+nonmet_reads.y)]
}

fread_and_merge <- function(dt, file){
  setkey(dt, chr, pos)
  fread(cmd=paste0("zcat < ", file), select=1:5, colClasses=list(factor=1L)) %>% 
    setnames(c("chr","pos","met_reads","nonmet_reads","rate")) %>%
    .[rate==0.50,c("met_reads","nonmet_reads"):=.(1L,1L)] %>%
    .[rate>0.50,c("met_reads","nonmet_reads"):=.(1L,0L)] %>%
    .[rate<0.50,c("met_reads","nonmet_reads"):=.(0L,1L)] %>%
    .[,rate:=NULL] %>%
    merge_and_sum(dt)
}

################
## Define I/O ##
################

io <- list()

if (grepl("ricard",Sys.info()['nodename'])) {
  io$in.data <- "/Users/ricard/data/gastrulation/acc/raw"
  io$in.metadata <- "/Users/ricard/data/gastrulation/sample_metadata_scNMT_lenient.txt"
  io$out.dir <- "/Users/ricard/data/gastrulation/acc/raw/pseudobulk"
} else {
  io$in.data <- "/hps/nobackup/stegle/users/ricard/gastrulation/acc/raw"
  io$in.metadata <- "/hps/nobackup/stegle/users/ricard/gastrulation/sample_metadata_scNMT_lenient.txt"
  io$out.dir <- "/hps/nobackup/stegle/users/ricard/gastrulation/acc/raw/pseudobulk"
}

####################
## Define options ##
####################

opts <- list()

# (mouse) chromosomes
opts$chr <- c(1:19,"X")

# Define pseudobulk groups in terms of stage_lineage
opts$groups <- c(
  "E4.5_EPI","E4.5_PE",
  "E5.5_EPI","E5.5_PE",
  "E6.5_EPI","E6.5_PS","E6.5_VE",
  "E7.5_Mesoderm","E7.5_Ectoderm","E7.5_Endoderm"
)

# Define which cells to use
opts$cells <- fread(io$in.metadata) %>% 
  .[,group:=paste(stage,lineage,sep="_")] %>%
  .[pass_accQC==T & group%in%opts$groups, id_acc]

# Parallel processing
opts$parallel <- F     # do parallel processing?
opts$ncores <- 1       # number of cores
opts$chunk_size <- 10  # chunk_size: the higher the less memory it is required????

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$in.metadata,stringsAsFactors=T)[,c("stage","lineage","id_acc")] %>% 
  .[id_acc%in%opts$cells] %>% setnames("id_acc","sample") %>%
  .[,group:=as.factor(paste(stage,lineage,sep="_"))] %>% setkey(group)

##############################
## Load data and pseudobulk ##
##############################


# Parallel processing options
if (opts$parallel){
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}


for (group in opts$groups) {
  
  # Define input files 
  cells <- sample_metadata[group,sample]
  files <- paste0(io$in.data, "/", cells, ".tsv.gz")
  
  # split into chunks for parallel processing
  if (opts$parallel) {
    chunks <- ceiling(seq_along(files)/opts$chunk_size)
    file_list <- split(files, chunks)
  } else {
    file_list <- list(files)
  }
  
  # pseudobulk
  init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_reads=as.integer(NA), nonmet_reads=as.integer(NA))
  data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=T) %>%
    purrr::reduce(merge_and_sum) %>%
    .[!is.na(chr)] %>% .[chr%in%opts$chr] %>%
    .[,rate:=round(100*met_reads/(met_reads+nonmet_reads))] %>%
    .[,.(chr,pos,met_reads,nonmet_reads,rate)]
  
  # Save
  outfile = sprintf("%s/%s.tsv",io$out.dir,group)
  fwrite(data, file=outfile, quote=F, col.names=T, sep="\t")
  system(sprintf("pigz -p %d -f %s",opts$ncores,outfile))
}
