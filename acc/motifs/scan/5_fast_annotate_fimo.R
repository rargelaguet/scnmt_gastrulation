##########################################################
## Script to quantify accessibility data at motif sites ##
##########################################################

library(data.table)
library(purrr)
library(furrr)

options(future.globals.maxSize =  2000 * 1024 ^ 2)

## Define I/O ##
io           <- list()
io$base_dir  <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir  <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/motifs"
io$raw_acc   <- paste0(io$base_dir, "/acc/raw/scNMT") # raw accessibility data
io$meta_data <- paste0(io$base_dir, "/sample_metadata.txt")
io$out_dir   <- paste0(io$base_dir, "/acc/parsed/fimo_motifs")



## Define options ##
opts                 <- list()
opts$stage           <- "all"           # select specific stages or "all"
opts$lineage         <- "all"           # select specific lineages or "all"
opts$anno_regex      <- "intersect12"   # pattern to match features 
opts$parallel        <- TRUE
opts$gzip            <- TRUE
opts$extend_anno_len <- FALSE
opts$anno_len        <- 100
opts$filter_fimo_p   <- 1e-6


## Functions ##
fread_gz <- function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") {
    filename <- paste("zcat", filename)
    return(fread(cmd = filename, ...))
  }
  fread(filename, ...)
}

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")

################################################################################


###############
## Load data ##
###############

# Load metadata and filter by QC
if (grepl("met", io$raw_acc)) {
  meta <- fread(io$meta_data) %>%
    .[pass_metQC == TRUE]
} else {
  meta <- fread(io$meta_data) %>%
    .[pass_accQC == TRUE]
}

# Some plates need re-naming for consistency ##

meta[, plate := gsub("Plate11", "E7.5_Plate11", plate) %>% 
       gsub("Plate12", "E7.5_Plate12", .) %>% 
       gsub("Plate13", "E7.5_Plate13", .) %>% 
       gsub("Plate14", "E7.5_Plate14", .) %>% 
       gsub("Plate15", "E7.5_Plate15", .) %>% 
       gsub("Plate16", "E7.5_Plate16", .) %>% 
       gsub("E6.5_late", "E6.75", .)]

## Filter by stage and lineage ##
if (opts$stage[1] != "all") meta <- meta[stage %in% opts$stage]
if (opts$lineage[1] != "all") meta <- meta[lineage %in% opts$lineage]


cells <- meta[, sample]

## Load annotation files ##
anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)]%>%
   map(fread, colClasses = list("character" = 1)) %>%
  rbindlist() %>%
  setkey(chr, start, end)

if (opts$extend_anno_len) {
  anno <- anno[, mid := start + (end-start)/2] 
  anno[end-start < opts$anno_len, c("start", "end") := .(as.integer(round(mid - opts$anno_len/2), round(mid + opts$anno_len/2)))]
  setkey(anno, chr, start, end)
}

## Filter for p-value ##
anno <- anno[fimo_p<=opts$filter_fimo_p]


## Find raw accessibility files ##
files <- meta[, sprintf("%s/%s/%s.tsv.gz", io$raw_acc, plate, sample)] %>%
  .[file.exists(.)]

cells <- basename(files) %>% gsub(".tsv.gz", "", .)

if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}

## Read in accessibility files and perform overlap and quantification ##
acc_dt <- future_map2(files, cells, ~{
#acc_dt <- map2(files, cells, ~{
  if (!file.exists(.x)) return("file not found")
  fread_gz(.x, select = c(1:2, 5), colClasses = list("factor" = 1L)) %>%
    setnames(c("chr", "start", "rate")) %>%
    .[, end := start] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>%
    .[, .(rate = round(100 * mean(rate)), .N, sample = .y), .(id_motif, tf)]
}) %>%
  rbindlist()

## Split by TF ##
setkey(acc_dt, tf)
acc_dt <- split(acc_dt, by = "tf")

#################
## Save output ##
#################

dir.create(io$out_dir, recursive = TRUE)

file_names <- paste0(io$out_dir, "/", names(acc_dt), ".tsv")
if (opts$extend_anno_len) file_names <- sub(".tsv", "_", opts$anno_len, "bp.tsv")

file_names <- gsub("\\(|\\)|\\.|::", "_", file_names) %>%
  gsub("_tsv", ".tsv", .) %>%
  gsub("_.tsv", ".tsv", .)

map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
