library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/EB_data_new_May2019"
io$in_dir <- paste0(io$data_dir, "/acc/differential/feature_level/")
# io$in_dir <- "/bi/home/clarks/gastrulation_out/acc/differential/feature_level/H3K27ac_distal_E7.5_Ect_intersect12/"

io$anno <- paste0(io$data_dir, "/features/filt")
io$out_dir <- paste0(io$data_dir, "/acc/differential/feature_level/fasta")

opts <- list()
opts$regex <- "E5.5E6.5Epi"
opts$up_or_down <- "down" # "up" "down" or "both"
opts$min_sites <- 100 # discard any sets of features with less than X rows
# opts$padj_cutoff <- 0.05

### functions ###

fread_gz <- function(path, ...){
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  if (ext == "gzfile") {
    return(fread(cmd = paste("zcat", path), ...))
  }
  fread(path, ...)  
}

get_seqs <- function(bed_dt, BSgenome){
  seq <- bed_dt[, getSeq(mm10, names = paste0("chr", chr), start = start, end = end)]
  seq@ranges@NAMES <- bed_dt[, as.character(id)]
  seq  
}

####################

anno <- dir(io$anno, full = TRUE, pattern = ".bed$") %>%
  .[grep(opts$regex, .)] %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno"))


dif_results <- dir(io$in_dir, full = TRUE, pattern = ".txt") %>%
  .[grep(opts$regex, .)] %>%
  set_names(basename(.) %>% gsub(".txt.gz|.txt", "", .)) %>%
  map(fread_gz) 

# insert name of annotation

names(dif_results) <- map2(dif_results, names(dif_results), ~{
  anno <- .x[, unique(anno)]
  nm <- .y
  if (!grepl(anno, .y)) nm <- paste(.y, anno, sep = "_")
  nm
})

dif_bed <- map(dif_results, ~.[sig == TRUE])

# select up or down accessibility (else include both)
if (toupper(opts$up_or_down) == "UP"){
  dif_bed <- map(dif_bed, ~.[diff < 0])
} else if (toupper(opts$up_or_down) == "DOWN") {
  dif_bed <- map(dif_bed, ~.[diff > 0])
}

dif_bed <- purrr::discard(dif_bed, ~nrow(.) < opts$min_sites) %>%
  map(~.[, .(anno, id)]) %>%
  map(merge, anno, by = c("anno", "id")) %>%
  map(~.[, .(chr, start, end, strand, anno, id)])



#load genome and find sequences
mm10 <- BSgenome.Mmusculus.UCSC.mm10


seqs <- map(dif_bed, get_seqs, mm10)

background_seqs <- split(anno, by = "anno") %>%
  map(get_seqs, mm10)

bed_anno_name <- map_chr(dif_bed, ~.[, unique(anno)]) 

out_dirs <- paste0(io$out_dir, "/", bed_anno_name, "/", opts$up_or_down)

out_files <- paste0(out_dirs, "/", names(seqs), ".fa")

walk(out_dirs, dir.create, recursive = TRUE)

# save fasta files

walk2(seqs, out_files, writeXStringSet)

# now save background fasta files

out_files <- paste0(io$out_dir, "/", names(background_seqs), "/background/", names(background_seqs), ".fa")

walk(out_files, ~dir.create(dirname(.), recursive = TRUE))

walk2(background_seqs, out_files, writeXStringSet)
