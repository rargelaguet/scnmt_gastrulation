library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno <- paste0(io$data_dir, "/features/filt")
io$out_dir <- paste0(io$data_dir, "/acc/differential/feature_level/fasta")

opts <- list()
opts$regex <- "E7.5_[A-z][a-z][a-z]_intersect12.bed$"


### functions ###

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



#load genome and find sequences
mm10 <- BSgenome.Mmusculus.UCSC.mm10


seqs <- get_seqs(anno, mm10)

annos <- anno[, unique(anno)] %>%
  paste(collapse="_")

out_file <- paste0(io$out_dir, "/", annos, ".fa")

writeXStringSet(seqs, out_file)
