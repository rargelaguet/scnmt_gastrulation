library(data.table)
library(purrr)
library(furrr)

options(future.globals.maxSize =  2000 * 1024 ^ 2)


# streamlined annotate script - requires more memory than ricard's script so may crash

io <- list()
io$base_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$base_dir,"/features/motifs/020519/200bp")
io$raw_acc <- paste0(io$base_dir, "/acc/raw/scNMT") # raw accessibility (or methylation) data
io$meta_data <- paste0(io$base_dir, "/sample_metadata.txt")
io$out_dir <- paste0(io$base_dir, "/acc/parsed/fimo_motifs/020519/200bp/select_few")


dir.create(io$out_dir, recursive = TRUE)

opts <- list()
opts$stage <- "all"#"E7.5" #"all"#"E7.5" #"all"
opts$lineage <- "all"
opts$anno_regex <- "Ect_intersect12|End_intersect12|Mes_intersect12"
opts$parallel <- T
opts$gzip <- TRUE
opts$extend_anno_len <- FALSE
opts$anno_len <- 100
opts$filter_fimo_p <- 1e-4


# fun
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


if (grepl("met", io$raw_acc)) {
  meta <- fread(io$meta_data) %>%
    .[pass_metQC == TRUE]
} else {
  meta <- fread(io$meta_data) %>%
    .[pass_accQC == TRUE]
}



meta[, plate := gsub("Plate11", "E7.5_Plate11", plate) %>% 
       gsub("Plate12", "E7.5_Plate12", .) %>% 
       gsub("Plate13", "E7.5_Plate13", .) %>% 
       gsub("Plate14", "E7.5_Plate14", .) %>% 
       gsub("Plate15", "E7.5_Plate15", .) %>% 
       gsub("Plate16", "E7.5_Plate16", .) %>% 
       gsub("E6.5_late", "E6.75", .)]


if (opts$stage[1] != "all") meta <- meta[stage %in% opts$stage]
if (opts$lineage[1] != "all") meta <- meta[lineage %in% opts$lineage]


cells <- meta[, sample]


anno <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)]%>%
#    map2(., sub(".bed", "", basename(.)), ~fread(.x, colClasses = list("factor" = 1L)) %>% 
#           setnames(c("chr", "start", "end", "strand", "id", "anno")) %>% 
#           .[, anno := paste0(anno, "_", .y)]) %>%
   map(fread, colClasses = list("character" = 1)) %>%
  rbindlist() %>%
#   setnames("tf", "anno") %>%
  setkey(chr, start, end)

if (opts$extend_anno_len) {
  anno <- anno[, mid := start + (end-start)/2] 
  anno[end-start < opts$anno_len, c("start", "end") := .(as.integer(round(mid - opts$anno_len/2), round(mid + opts$anno_len/2)))]
  setkey(anno, chr, start, end)
}

anno <- anno[fimo_p<=opts$filter_fimo_p]



files <- meta[, sprintf("%s/%s/%s.tsv.gz", io$raw_acc, plate, sample)] %>%
  .[file.exists(.)]

cells <- basename(files) %>% gsub(".tsv.gz", "", .)

anno=anno[tf %in% c("Sox2", "Pou5f1::Sox2", "SP8", "FOXA1",   "Foxa2", "Sox17", "Gata1",  "Gata4", "Hand1::Tcf3", "MEIS1" , "TWIST1" )]


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
    .[, .(rate = round(100 * mean(rate)), .N, sample = .y), .(id_motif, tf)]
}) %>%
  rbindlist()

setkey(acc_dt, tf)
acc_dt <- split(acc_dt, by = "tf")

# save one file per anno
file_names <- paste0(io$out_dir, "/", names(acc_dt), ".tsv")
if (opts$extend_anno_len) file_names <- sub(".tsv", "_", opts$anno_len, "bp.tsv")


file_names <- gsub("\\(|\\)|\\.|::", "_", file_names) %>%
  gsub("_tsv", ".tsv", .) %>%
  gsub("_.tsv", ".tsv", .)

# future_map2(acc_dt, file_names, fwrite, sep = "\t")
map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
