library(data.table)
library(purrr)
library(furrr)


options(future.globals.maxSize = 1000 * 1024 ^ 2)

# simple script to quantify accessibility rates at overlapping windows in 
# pseudobulk data then pick the top x% to use as an annotation
# for small window sizes (<500) it is quite memory hungry

## in/out ##
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$raw_acc_psuedobulk_file <- paste0(io$data_dir, "/acc/raw/pseudobulk/E7.5_all.tsv.gz")
io$annotation_out <- paste0(io$data_dir, "/features/peaks/march19")

## options ##
opts <- list()
opts$window <- 100
opts$step <- 50
opts$min_weight <- 100
opts$merge_overlapping_windows <- TRUE # if adjacent windows have similar rates then merge
opts$max_diff <- 10 # maximum rate difference between adjacent sites to merge sites
opts$cutoff <- 0.05 # use the top (or bottom) x fraction of sites 
opts$top_or_bottom <- "top" # pick the most ("top") or least ("bottom") accessible sites

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

# load data

data <- fread_gz(io$raw_acc_psuedobulk_file)

# if loading bismark cov then set col names....
if (grepl(".cov", io$raw_acc_psuedobulk_file)) {
  setnames(data, c("chr", "pos", "end", "rate", "met_reads", "nonmet_reads"))
  data <- data[, .(chr, pos, rate, met_reads, total_reads = met_reads + nonmet_reads)]
}

data <- data[!chr %in% "MT"]

windows <- data[, .(start = seq(min(pos), max(pos), by = opts$step)), chr] %>%
  .[, end := start + opts$window] %>%
  setkey(chr, start, end)

data <- setnames(data, "pos", "start") %>%
  .[, end := start] %>%
  setkey(chr, start, end) 


# split by chr to process as chunks

data <- split(data, by = "chr")
windows <- split(windows, by = "chr")

plan(multiprocess)
data <- future_map2(data, windows, foverlaps, nomatch=0L) %>%
  rbindlist() %>%
  .[, .(met_reads = sum(met_reads), total_reads = sum(total_reads)), .(chr, start, end)] %>%
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
          total_reads = sum(total_reads)), .(chr, group)] %>%
    .[, rate := round(100*met_reads/total_reads)]
}




setkey(data, total_reads)
anno <- data[total_reads>=opts$min_weight] 

if (opts$top_or_bottom == "bottom") {
  anno <- anno[order(rate)]
} else {
  anno <- anno[order(-rank(rate))]
}
anno <- anno[1:(.N*opts$cutoff)]



out_file <- paste0(io$annotation_out, 
                   "/", 
                   basename(io$raw_acc_psuedobulk_file) %>% gsub(".tsv.gz|.cov.gz", "", .),
                   "_", 
                   opts$window, 
                   "_by_", 
                   opts$step, 
                   "_",
                   opts$top_or_bottom,
                   "_",
                   opts$cutoff,
                   "_filt_",
                   opts$min_weight,
                   ".bed")

if (opts$merge_overlapping_windows) out_file <- gsub(".bed", "_merged_overlap.bed", out_file)

# re-format and save

anno_name <- basename(out_file) %>% gsub(".bed", "", .) %>% paste0("peaks_", .)
anno <- anno[, c("strand", "id", "anno") := .("*", paste0(anno_name, "_", .I), anno_name)]
anno <- anno[, .(chr, start, end, strand, id, anno)]

dir.create(io$annotation_out, recursive=TRUE)
fwrite(anno, out_file, sep="\t")


