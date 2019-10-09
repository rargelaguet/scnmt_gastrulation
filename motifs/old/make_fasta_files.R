library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)

# converts bed files to fasta files

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
#io$in_dir <- paste0(io$data_dir, "/features/sjc/significant_sites/", c("met", "acc"))
# io$in_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/bed")
# io$out_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/fasta")

# io$in_dir <- paste0(io$data_dir, "/features/sjc/sc_acc")
# io$out_dir <- paste0(io$data_dir, "/features/sjc/sc_acc/fasta")

# io$in_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/bed/met_low_acc_high")
# io$out_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/fasta/met_low_acc_high")

io$in_dir <- paste0(io$data_dir, "/acc/differential_results/bed")

#io$in_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/bed/epi_pe")
# io$in_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/bed")
# io$out_dir <- paste0(io$data_dir, "/features/ricard_significant_sites/fasta")
#io$out_dir <- paste0(io$data_dir, "/features/ChIP/fasta")
io$out_dir <- paste0(io$data_dir, "/acc/differential_results/fasta")

opts <- list()
opts$anno_regex <- ""

dir.create(io$out_dir, recursive = TRUE)

bed <- dir(io$in_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(~fread(.)[, file := basename(.)]) %>%
  rbindlist() %>%
#   setnames(c("chr", "start", "end", "file")) %>%
#   .[, anno := sub(".bed", "", file)] %>%
#   .[, id := paste(anno, .I, sep = "_")] %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno", "file")) %>%
  #.[, anno := gsub("E7.5", "", anno) %>% gsub("__", "_", .) %>% strsplit("_") %>% map_chr(~paste(.[1], .[2], sep = "_"))] %>%
  split(by = c("file"), keep.by = FALSE)

mm10 <- BSgenome.Mmusculus.UCSC.mm10

seqs <- map(bed, ~{
  seq <- .[, getSeq(mm10, names = paste0("chr", chr), start = start, end = end)]
  seq@ranges@NAMES <- .[, as.character(id)]
  seq
})

out_files <- paste0(io$out_dir, "/", names(seqs) %>% sub("\\.", "_", .), ".fa")
walk2(seqs, out_files, writeXStringSet)
