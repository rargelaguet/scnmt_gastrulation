library(data.table)
library(purrr)
library(furrr)


# generate pseudobulk accessibility data for use in e.g. peak_caller.R 
# only generates one file per go


## in/out ##
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$sample_meta <- paste0(io$data_dir, "/sample_metadata.txt")
io$raw_acc_dir <- paste0(io$data_dir, "/acc/raw/scNMT")
io$out_dir <- paste0(io$data_dir, "/acc/raw/pseudobulk/new_24-05-2019")

dir.create(io$out_dir, recursive = TRUE)

## options ##
opts <- list()

opts$stages <- "E7.5"#"all"#c("E4.5", "E5.5", "E6.5", "E7.5")
opts$lineages <- "Epiblast" #"Primitive_Streak" # "all"#"PE"#"Ectoderm" #"all" #"Endoderm" #"all"
opts$parallel <- TRUE
opts$chunk <- 20 # for parallel processing split input files into chunks of this many (smaller chunk size requires more memory)
opts$format <- "cov" # cov = <chr><start><end><rate><met_reads><nonmet_reads> defaults to <chr><pos><met_reads><total_reads>
opts$gzip <- TRUE


## functions ##

merge_and_sum <- function(dt1, dt2){
  walk(list(dt1, dt2), setkey, chr, pos)
  merge(dt1, dt2, all = TRUE) %>%
    .[is.na(met_reads.x), met_reads.x := 0] %>%
    .[is.na(met_reads.y), met_reads.y := 0] %>%
    .[is.na(nonmet_reads.x), nonmet_reads.x := 0] %>%
    .[is.na(nonmet_reads.y), nonmet_reads.y := 0] %>%
    .[, .(chr=chr, pos=pos,
          met_reads=met_reads.x+met_reads.y,
          nonmet_reads=nonmet_reads.x+nonmet_reads.y)]
}

fread_gz <- function(path, ...){
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  if (ext == "gzfile") {
    return(fread(cmd = paste("zcat", path), ...))
  }
  fread(path, ...)  
}

fread_and_merge <- function(dt, file){
  setkey(dt, chr, pos)
  fread_gz(file, select=1:4, colClasses=list(character=1L)) %>%
    merge_and_sum(dt)
}

fwrite_tsv <- partial(fwrite, sep = "\t")


# load metadata and find files


meta <- fread(io$sample_meta)  %>%
  .[pass_accQC==TRUE] %>%
  .[, lineage := lineage10x_2]

if (opts$stages != "all") meta <- meta[stage %in% toupper(opts$stages)]
if (opts$lineages != "all") meta <- meta[lineage %in% opts$lineages]

cells <- meta[, sample]
all_files <- dir(io$raw_acc_dir, pattern = ".tsv.gz$", full = TRUE, recursive = TRUE)
files <- all_files[basename(all_files) %in% paste0(cells, ".tsv.gz")]


# split into chunks for parallel processing
chunks <- ceiling(seq_along(files)/opts$chunk)
file_list <- split(files, chunks)





if (opts$parallel){
  plan(multiprocess)
} else {
  plan(sequential)
}

init <- data.table(
  chr=as.character(NA), 
  pos=as.integer(NA), 
  met_reads=as.integer(NA), 
  nonmet_reads=as.integer(NA)
)

data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
  purrr::reduce(merge_and_sum) %>%
  .[!is.na(chr)] %>%
  .[, total_reads := .(met_reads+nonmet_reads)] %>%
  .[, rate := round(100 * met_reads / total_reads)]


# save file

stages <- paste(opts$stages, collapse = "_")
lineages <- paste(opts$lineages, collapse = "_")
out_file <- paste0(io$out_dir, "/", stages, "_", lineages)


if (opts$format == "cov" || opts$format == "bismark") {
  data <- data[, .(chr, pos, pos, rate, met_reads, nonmet_reads)]
  out_file <- paste0(out_file, ".cov")
  fwrite_tsv(data, out_file, col.names = FALSE)
  
} else {
  data <- data[, .(chr, pos, met_reads, total_reads)]
  out_file <- paste0(out_file, ".tsv")
  fwrite_tsv(data, out_file)
}

if (opts$gzip) system(paste("gzip -f", out_file))









