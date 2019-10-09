library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$in_dir <- paste0(io$data_dir, "/features/filt")
io$out_dir <- paste0(io$data_dir, "/features/fasta")


opts <- list()


dir.create(io$out_dir, recursive = TRUE)

bed <- dir(io$in_dir, pattern = ".bed$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(~fread(.)[, file := basename(.)]) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno", "file")) %>%
  split(by = c("file"), keep.by = FALSE)

mm10 <- BSgenome.Mmusculus.UCSC.mm10

seqs <- map(bed, ~{
  seq <- .[, getSeq(mm10, names = paste0("chr", chr), start = start, end = end)]
  seq@ranges@NAMES <- .[, as.character(id)]
  seq
})

out_files <- paste0(io$out_dir, "/", names(seqs) %>% sub(".bed", ".fa", .))
walk2(seqs, out_files, writeXStringSet)
