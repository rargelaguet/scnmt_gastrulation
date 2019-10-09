library(data.table)
library(purrr)
library(furrr)

options(future.globals.maxSize =  2000 * 1024 ^ 2)

# find pairwise combinations of motifs


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$motif_bed <- paste0(io$data_dir,"/features/motifs/020519/200bp")
io$raw_acc <- paste0(io$data_dir, "/acc/raw/scNMT") # raw accessibility (or methylation) data
io$meta_data <- paste0(io$data_dir, "/sample_metadata.txt")
io$out_dir <- paste0(io$data_dir, "/acc/parsed/fimo_motifs/motif_pairs")

dir.create(io$out_dir, recursive = TRUE)
# io$out <- 

opts<-list()
opts$motfis <- c("SP8", "SOX3", "TWIST1", "GATA4", "FOXA1", "HNF1B", "POU3F1", "FOXA2", "SOX17", "HAND1")
opts$max_dist <- 500
opts$extend_bp <- 50

opts$stage <- "all"#"E7.5" #"all"#"E7.5" #"all"
opts$lineage <- "all"
opts$parallel <- T
opts$gzip <- TRUE

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

find_anno <- function(id){
  sp <- strsplit(id, "_")
  map_chr(sp, ~.[c(1:4)] %>% paste(collapse = "_"))
}

################################################################################
# find pairs

pairs <- expand.grid(opts$motfis, opts$motfis) %>%
  setDT() %>%
  .[Var1 != Var2]

motifs <- dir(io$motif_bed, full = TRUE, pattern = ".bed$") %>%
  map(fread) %>%
  rbindlist() %>%
  .[, tf := toupper(tf)]

id_chr <- motifs[, .(id, chr)] %>%
  setkey() %>%
  unique()

loci_of_pairs <- map(1:nrow(pairs), ~{
  
  motif_pair <- unlist(pairs[.])
  dt <- motifs[tf %in% motif_pair, .(tf, id, mid)] %>%
    split(by = "tf")
  
  if (length(dt)!= 2) return(NULL)
  
    #     map2(names(.), ~setnames(.x, "mid", paste0("mid_", .y))) %>%
  purrr::reduce(dt, merge, by = "id", allow.cartesian = TRUE) %>%
    .[, dist := mid.x - mid.y] %>%
    .[abs(dist) < opts$max_dist] %>%
    .[, c("start", "end") := .(apply(.SD, 1, min), apply(.SD, 1, max)), .SDcol = c("mid.x", "mid.y")] %>%
    .[, c("start", "end") := .(start - opts$extend_bp, end + opts$extend_bp)] %>%
    merge(id_chr, by = "id") %>%
    .[, .(id, tf = paste(tf.x, tf.y, sep = "_"), chr, start, end)]
}) %>%
  compact() %>%
  rbindlist() %>%
  .[, anno := find_anno(id)]



################################################################################
# now quantify acc data



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

files <- meta[, sprintf("%s/%s/%s.tsv.gz", io$raw_acc, plate, sample)] %>%
  .[file.exists(.)]

cells <- basename(files) %>% gsub(".tsv.gz", "", .)

anno <- loci_of_pairs %>%
  setkey(chr, start, end)


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
    .[, .(rate = round(100 * mean(rate)), .N, sample = .y), .(id, anno, tf)]
}) %>%
  rbindlist()

setkey(acc_dt, tf)
acc_dt <- split(acc_dt, by = "tf")

# save one file per anno
file_names <- paste0(io$out_dir, "/", names(acc_dt), ".tsv")



file_names <- gsub("\\(|\\)|\\.|::", "_", file_names) %>%
  gsub("_tsv", ".tsv", .) %>%
  gsub("_.tsv", ".tsv", .)

# future_map2(acc_dt, file_names, fwrite, sep = "\t")
map2(acc_dt, file_names, fwrite, sep = "\t")

if (opts$gzip) walk(file_names, ~system(paste("gzip -f", .)))
