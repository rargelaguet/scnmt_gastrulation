library(data.table)
library(purrr)
library(furrr)


options(future.globals.maxSize = 1000 * 1024 ^ 2)

# simple script to quantify accessibility rates at overlapping windows in 
# pseudobulk data WITHIN a given set of annotations, then refine the coordinates
# of the annotation to only include the peak

## in/out ##
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$raw_acc_psuedobulk_file <- paste0(io$data_dir, "/acc/raw/pseudobulk/E7.5_all.tsv.gz")
io$annotation_in <- paste0(io$data_dir, "/features/filt")
io$annotation_out <- paste0(io$data_dir, "/features/peaks/H3K27ac")

## options ##
opts <- list()
opts$anno_regex <- "intersect12"
opts$window <- 100
opts$step <- 20
opts$min_weight <- 50
opts$merge_overlapping_windows <- TRUE # if adjacent windows have similar rates then merge
opts$max_diff <- 10 # maximum rate difference between adjacent sites to merge sites


## functions ##
fread_gz <- function(path, ...){
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  if (ext == "gzfile") {
    return(fread(cmd = paste("zcat", path), ...))
  }
  fread(path, ...)  
}

fwrite_tsv <- partial(fwrite, sep = "\t")

# load data and subset to annotation

anno <- dir(io$annotation_in, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno")) %>%
  setkey(chr, start, end)

data <- fread_gz(io$raw_acc_psuedobulk_file) 

# if loading bismark cov then set col names....
if (grepl(".cov", io$raw_acc_psuedobulk_file)) {
  setnames(data, c("chr", "pos", "end", "rate", "met_reads", "nonmet_reads"))
  data <- data[, .(chr, pos, rate, met_reads, total_reads = met_reads + nonmet_reads)]
}

data <- data[!chr %in% "MT"] %>%
  .[, .(chr, start = pos, end = pos, met_reads, total_reads)] %>%
  setkey(chr, start, end) %>%
  foverlaps(anno, nomatch = 0L) %>%
  .[, .(chr, id, anno, pos = i.start, met_reads, total_reads)]

windows <- data[, .(start = seq(min(pos), max(pos), by = opts$step) %>% round), .(chr ,id, anno)] %>%
  .[, end := start + opts$window] %>%
  setkey(chr, start, end)

data <- data[, .(chr, start = pos, end = pos, met_reads, total_reads)] %>%
  setkey(chr, start, end) 


# split by chr to process as chunks

data <- split(data, by = "chr")
windows <- split(windows, by = "chr")

plan(multiprocess)
data <- future_map2(data, windows, foverlaps, nomatch=0L) %>%
  rbindlist() %>%
  .[, .(met_reads = sum(met_reads), total_reads = sum(total_reads)), .(chr, start, end, id, anno)] %>%
  .[, rate := round(100*met_reads/total_reads)]


# combine overlapping windows if acc rate is the same 

if (opts$merge_overlapping_windows){
  data <- data[order(chr, start)] %>%
    .[, dif_r := rate - shift(rate, type = "lag")] %>%
    .[is.na(dif_r), dif_r := 100] %>%
    .[, dif_s := shift(start, type = "lead") - start] %>%
    .[, not_overlap := TRUE] %>%
    .[dif_s < opts$window, not_overlap := FALSE]  %>% 
    .[, group := as.numeric(abs(dif_r) > opts$max_diff | not_overlap) %>% cumsum] %>%
    .[, .(start = min(start), 
          end = max(end), 
          met_reads = sum(met_reads),
          total_reads = sum(total_reads)), .(chr, group, id, anno)] %>%
    .[, rate := round(100*met_reads/total_reads)]
}




setkey(data, total_reads)
anno <- data[total_reads>=opts$min_weight] %>%
  .[order(id, -rank(rate))] %>%
  .[, head(.SD, 1), id] %>%
  .[, .(chr, start, end, strand = "*", id, anno)] %>%
  split(by = "anno")


out_files <- paste0(io$annotation_out, 
                   "/", 
                   names(anno),
                   "_peak_", 
                   opts$window, 
                   "_by_", 
                   opts$step, 
                   ".bed")



# re-format and save


dir.create(io$annotation_out, recursive=TRUE)


walk2(anno, out_files, fwrite_tsv)

