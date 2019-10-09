library(data.table)
library(purrr)
library(furrr)


# streamlined annotate script 

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$data_dir, "/features/peaks/H3K27ac")
io$raw_acc <- paste0(io$data_dir, "/acc/raw/scNMT") # raw accessibility (or methylation) data
io$meta_data <- paste0(io$data_dir, "/sample_metadata.txt")
io$out_dir <- paste0(io$data_dir, "/acc/parsed/peaks/H3K27ac")


dir.create(io$out_dir, recursive = TRUE)

opts <- list()
opts$stage <- "all"#"E7.5" #"all"#"E7.5" #"all"
opts$lineage <- "all"
#opts$KO_3b <- "not"
opts$anno_regex <- ""
opts$parallel <- TRUE
opts$gzip <- TRUE
opts$extend_anno_len <- FALSE
opts$anno_len <- 500


# fun
fread_gz = function(filename, ...){
  paste("zcat", filename) %>%
    fread(cmd = ., ...)
}


if (grepl("met", io$raw_acc)) {
  meta <- fread(io$meta_data) %>%
    .[pass_metQC == TRUE]
} else {
  meta <- fread(io$meta_data) %>%
    .[pass_accQC == TRUE]
}



if (opts$stage[1] != "all") meta <- meta[stage %in% opts$stage]
if (opts$lineage[1] != "all") meta <- meta[lineage %in% opts$lineage]
#if (opts$KO_3b[1] != "all") meta <- meta[KO_3b %in% opts$KO_3b]

cells <- meta[, sample]
plates <- meta[, plate]


anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)]%>%
#    map2(., sub(".bed", "", basename(.)), ~fread(.x, colClasses = list("factor" = 1L)) %>% 
#           setnames(c("chr", "start", "end", "strand", "id", "anno")) %>% 
#           .[, anno := paste0(anno, "_", .y)]) %>%
   map(fread) %>%
  rbindlist() %>%
  setkey(chr, start, end)

if (opts$extend_anno_len) {
  anno <- anno[, mid := start + (end-start)/2] 
  anno[end-start < opts$anno_len, c("start", "end") := .(as.integer(round(mid - opts$anno_len/2), round(mid + opts$anno_len/2)))]
  setkey(anno, chr, start, end)
}


files <- paste0(io$raw_acc, "/", plates, "/", cells, ".tsv.gz")
cells <- cells[file.exists(files)]
files <- files[file.exists(files)]




if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}

acc_dt <- future_map2(files, cells, ~{
#acc_dt <- map2(files, cells, ~{
  if (!file.exists(.x)) return("file not found")
  fread_gz(.x, select = c(1:2, 5), colClasses = list("factor" = 1L)) %>%
    setnames(c("chr", "start", "rate")) %>%
    .[, end := start] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>%
    .[, .(sample = .y, total_reads = .N, rate = round(100 * mean(rate))), .(id, anno)] %>%
    .[, .(sample, id, anno, met_reads = round(total_reads * rate/100), total_reads, rate)]
}) %>%
  rbindlist()

setkey(acc_dt, anno)
acc_dt <- split(acc_dt, by = "anno")

# save one file per anno
file_names <- paste0(io$out_dir, "/", opts$stage, "_" ,names(acc_dt), ".tsv")
if (opts$extend_anno_len) file_names <- sub(".tsv", "_", opts$anno_len, "bp.tsv")

# future_map2(acc_dt, file_names, fwrite, sep = "\t")
map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
