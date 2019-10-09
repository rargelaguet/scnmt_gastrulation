library(data.table)
library(purrr)



io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$fimo_dir <- paste0(io$data_dir, "/features/fimo")
io$pwm_file <- paste0(io$data_dir, "/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt")
io$anno_dir <- paste0(io$data_dir, "/features/filt")
io$out_dir <- paste0(io$data_dir, "/features/sjc/motifs/fimo/bed/q_filt")

opts <- list()
opts$fimo_regex <- "intersect12"
opts$win_size <- 100
opts$fimo_q_threshold <- 0.1 # note that fimo results are already filtered on p-value


fread_gz <- function(file, ...){
  paste("zcat", file) %>%
    fread(...)
}

anno <- dir(io$anno_dir, pattern = "K27ac", full = TRUE) %>%
  #map2(basename(.) %>% gsub(".bed", "", .), ~fread(.x)[, anno := .y] %>% .[, id := paste0(anno, "_", .I)]) %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

#anno <- anno[anno %like% "intersect12"]

motif_names <- readLines(io$pwm_file) %>%
  .[grep("MOTIF", .)] %>%
  as.data.table() %>%
  tidyr::separate(".", c("null", "motif", "tf"), sep = " ") %>%
  .[, .(motif, tf)]

motif_anno <- dir(io$fimo_dir, pattern = opts$fimo_regex, full = TRUE) %>%
  dir(pattern = "fimo.txt", full = TRUE) %>%
  map(fread, select = c(1:4, 7:8)) %>%
  rbindlist() %>%
  setnames(c("motif", "id", "min", "max", "p", "q")) %>%
  .[q < opts$fimo_q_threshold] %>%
  merge(anno, by = "id") %>%
  .[, c("start", "end") := .(start + min, start + max)] %>%
  .[, mid := start + (start-end)/2] %>%
  .[, c("start", "end") := .(mid - opts$win_size/2, mid + opts$win_size/2)] %>%
  merge(motif_names, by = "motif") %>%
  setkey(id) %>%
  .[, id := make.names(id, unique = TRUE) %>% gsub("\\.", "_", .)] %>%
  .[, .(chr, start, end, tf, id)] %>%
  setkey(chr, start, end)

file_name <- dir(io$fimo_dir, pattern = opts$fimo_regex) %>%
  paste(collapse = "_") %>%
  gsub(".fa|.bed", "", .) %>%
  paste0(io$out_dir, "/", ., ".bed")

dir.create(io$out_dir, recursive = TRUE)

fwrite(motif_anno, file_name, sep = "\t")
