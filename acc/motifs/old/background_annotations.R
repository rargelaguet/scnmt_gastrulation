library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)

# generate random annotations to match a given set of annotations (matched by size and chromosome)

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb/filt")
io$out_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb/filt/background")

dir.create(io$out_dir, recursive=TRUE)

background_anno <- dir(io$anno_dir, full=TRUE, pattern=".bed") %>% 
  map(fread) %>%
  rbindlist() %>% 
  .[, size := end-start] %>%
  .[, c("min", "max") := .(min(start), max(end)), chr] %>% 
  .[, .(start = runif(1, min, max)), .(id, anno, chr, size)] %>% 
  .[, end := start+size] %>%
  split(by="anno")

walk2(background_anno, names(background_anno), ~{
  file_name <- paste0(io$out_dir, "/", .y, ".bed")
  bed_file <- .x[, .(chr, start, end, strand="*", id, anno)]
  fwrite(bed_file, file_name, sep="\t")
})





