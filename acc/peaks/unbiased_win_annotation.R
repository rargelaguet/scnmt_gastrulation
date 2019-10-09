library(data.table)
library(purrr)
library(furrr)

options(future.globals.maxSize = 2000 * 1024 ^ 2)

# generate an annotation of overlapping windows then quantify accessibility data

## in/out ##
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$sample_meta <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$raw_acc <- paste0(io$data_dir, "/acc/raw")
io$acc_pseudobulk <- paste0(io$data_dir, "/acc/raw/pseudobulk/e4.5.tsv.gz") # provides min and max bp of each chr 
io$out_dir <- paste0(io$data_dir, "/acc/parsed/unbiased/e4.5")

## options ##
opts <- list()
opts$window <- 100
opts$step <- 50
opts$stage <- "E4.5"#c("E4.5", "E5.5", "E6.5", "E7.5")
opts$lineages <- c("EPI", "PE") # c("Mesoderm", "Endoderm", "Ectoderm", "EPI", "PS")
opts$KO_3b <- "not"
opts$chunks <- 100 # files will be split into x chunks so that later processing can be done in chunks
opts$parallel <- FALSE # parallel processing by cells?
opts$gzip <- TRUE

# functions

fread_gz <- function(filename, ...){
  fread(cmd = paste("zcat", filename), ...)
}


## load metadata and find files

meta <- fread(io$sample_meta) %>% 
  .[stage %in% opts$stage & 
      lineage %in% opts$lineage & 
      KO_3b %in% opts$KO_3b & 
      pass_accQC == TRUE ]


stages <- paste0("_stages_", paste(opts$stage, collapse = "_"))
out_dir <- paste0(io$out_dir, "/win", opts$window, "_step", opts$step )
dir.create(paste0(out_dir, "/tmp/"), recursive = TRUE)
anno_name <- paste0(opts$window, "win_step", opts$step)

# choose cells
cells <- meta[, sample]

# check if files exists before regenrating...
file_base <- paste0(out_dir, "/tmp/chunk", 0:opts$chunks) 
file_check <- map(cells, ~paste0(file_base, "_", ., ".tsv")) %>%
  map_lgl(~file.exists(.) %>% all)

cells <- cells[!file_check]


# find raw files
files <- paste0(io$raw_acc, "/", cells, ".tsv.gz")

cells <- cells[file.exists(files)]
files <- files[file.exists(files)]

# generate annotation
anno <- fread_gz(io$acc_pseudobulk, select = 1:2) %>%
  .[, .(min = min(pos), max = max(pos)), chr] %>%
  .[, .(start = seq(min, max, by = opts$step)), chr] %>%
  .[, end := start + opts$window] %>%
  .[, id := .I] %>%
  .[, chunk := round(id/(.N/opts$chunk))] %>%
  setkey(chr, start, end)

# save annotation
fwrite(anno[, .(chr, start, end, strand = "*", id, anno = anno_name)], paste0(out_dir, "/annotation.bed"), sep = "\t")



# parallel processing
if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}


# #testing
# anno <- sample_n(anno, 1e6) %>%
#   setkey(chr, start, end)
# files <- files[1:4]
# cells <- cells[1:4]


tmp_files <- future_map2(files, cells, ~{
  dt <- fread_gz(.x, select = c(1:2, 5)) %>%
    .[, .(chr, start = pos, end = pos, rate = 100 * rate)] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch = 0L) %>%
    .[, .(rate = round(mean(rate)), weight = .N, sample = .y, anno = anno_name), .(id, chunk)] %>%
    setcolorder(c("sample", "id", "anno", "rate", "weight", "chunk")) %>%
    split(by = "chunk", keep.by = FALSE)
  
  out_files <- paste0(out_dir, "/tmp/chunk", names(dt), "_",.y, ".tsv")
  walk2(dt, out_files, fwrite, sep = "\t", col.names = FALSE)
  
  out_files
}) 

## append files (keeping them split by chunk)

tmp_files <- unlist(tmp_files)

#tmp_files <- dir(paste0(out_dir, "/tmp"), pattern = "chunk", full = TRUE)

chunks <- paste0("chunk", 0:opts$chunks, "_")

chunk_files <- map(chunks, ~tmp_files[grep(., tmp_files)])





stages <- paste(opts$stage, collapse = "_")

out_files <- paste0(out_dir, "/parsed_stage", stages, "/" ,chunks, ".tsv")
dir.create(dirname(out_files)[1], recursive = TRUE)

walk2(chunk_files, out_files, ~{
  file.copy(.x[1], .y)
  walk(.x[2:length(.x)], function(x) file.append(.y, x))
})

file.remove(tmp_files)

if (opts$gzip) walk(out_files, ~system(paste("gzip -f ", .)))

