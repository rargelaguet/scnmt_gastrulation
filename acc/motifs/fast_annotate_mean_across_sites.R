library(data.table)
library(purrr)
library(furrr)
library(dplyr)

options(future.globals.maxSize = 100000*1024^2)

# streamlined annotate script - requires more memory than ricard's script so may crash

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb_100bp/filt_by_peaks")
io$raw_acc <- paste0(io$data_dir, "/acc/raw") # raw accessibility (or methylation) data
io$meta_data <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$tf_motifs <- paste0(io$data_dir, "/features/sjc/significant_sites/motifs") # only needed if opts$filter_tfs==TRUE
io$out_dir <- paste0(io$data_dir, "/acc/parsed/motifdb_100bp/filt_by_peaks/all_loci")

dir.create(io$out_dir, recursive = TRUE)

opts <- list()
opts$stage <- "all"#"E7.5" #"all"#"E7.5" #"all"
opts$lineage <- "all"
opts$KO_3b <- "not"
opts$filter_tfs <- FALSE
opts$parallel <- TRUE
opts$gzip <- TRUE
opts$extend_anno_len <- FALSE
opts$anno_len <- 100
opts$mean_across_sites <- FALSE
opts$n_sites_per_anno <- FALSE#1e4


# fun
fread_gz = function(filename, ...){
  paste("zcat", filename) %>%
    fread(...)
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
if (opts$KO_3b[1] != "all") meta <- meta[KO_3b %in% opts$KO_3b]

cells <- meta[, sample]


# anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
#   .[grep(opts$anno_regex, .)]%>%
#    map2(., sub(".bed", "", basename(.)), ~fread(.x, colClasses = list("factor" = 1L)) %>% 
#           setnames(c("chr", "start", "end", "strand", "id", "anno")) %>% 
#           .[, anno := paste0(anno, "_", .y)]) %>%
#   #  map(fread, colClasses = list("character" = 1)) %>%
#   rbindlist() %>%
#   setkey(chr, start, end)

if (opts$filter_tfs) {
  sig_tfs <- dir(io$tf_motifs, pattern = "ame.txt", full = TRUE, recursive = TRUE) %>%
    map2(basename(dirname(.)), ~{
      fread(.x) %>%
        .[, .(type = .y, tf = V9)]
    }) %>%
    rbindlist()
  
  anno_regex <- sig_tfs[, unique(tf)] %>%
    gsub("::", "|", .) %>%
    gsub("\\(.*", "", .) %>%
    paste(collapse = "|")
} else {
  anno_regex <-  ""
}



anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(toupper(anno_regex), toupper(.))]%>%
  map(~{
    dt <- fread(., colClasses = list("factor" = 1L))
    if (opts$n_sites_per_anno != FALSE) {
      dt <- dt[order(-rank(score))] %>%
        .[1:min(opts$n_sites_per_anno, .N)]
    }
    dt
  }) %>%
  rbindlist() %>%
  setkey(chr, start, end)

if (opts$extend_anno_len) {
  anno <- anno[, mid := start + (end-start)/2] 
  anno[end-start < opts$anno_len, c("start", "end") := .(as.integer(round(mid - opts$anno_len/2), round(mid + opts$anno_len/2)))]
  setkey(anno, chr, start, end)
}


files <- paste0(io$raw_acc, "/", cells, ".tsv.gz")
cells <- cells[file.exists(files)]
files <- files[file.exists(files)]




if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}

acc_dt <- future_map2(files, cells, ~{
    if (!file.exists(.x)) return("file not found")
    
    dt <- fread_gz(.x, select = c(1:2, 5), colClasses = list("factor" = 1L)) %>%
      setnames(c("chr", "start", "rate")) %>%
      .[, end := start] %>%
      setkey(chr, start, end) %>%
      foverlaps(anno, nomatch = 0L) 
    
    if (opts$mean_across_sites){
      dt[, .(rate = round(100 * mean(rate)), .N, sample = .y), .(anno)]
    } else {
      dt[, .(rate = round(100 * mean(rate)), .N, sample = .y), .(id, anno)]
    }
    
    
}) %>%
    rbindlist()
  
setkey(acc_dt, anno)
acc_dt <- split(acc_dt, by = "anno")

# save one file per anno
file_names <- paste0(io$out_dir, "/", opts$stage, "_" ,names(acc_dt), ".tsv")
if (opts$extend_anno_len) file_names <- sub(".tsv", paste0("_", opts$anno_len, "bp.tsv"), file_names)
if (opts$mean_across_sites) file_names <- sub(".tsv", "_mean_across_sites.tsv", file_names)

file_names <- gsub("\\(|\\)", "_", file_names) %>%
  gsub(" ", "_", .)

# future_map2(acc_dt, file_names, fwrite, sep = "\t")
map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
