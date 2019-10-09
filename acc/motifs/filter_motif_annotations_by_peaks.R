library(data.table)
library(purrr)
library(furrr)


# in/out
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb_short")
io$filter_bed_file <- paste0(io$data_dir, "/features/sjc/peaks/E4.5_all_100_by_50_top_0.05_filt_100_merged_overlap.bed")
io$anno_out <- paste0(io$data_dir, "/features/sjc/motifs/motifdb_100bp/filt_by_peaks")

dir.create(io$anno_out, recursive = TRUE)

opts <- list()
opts$anno_regex <- ""
opts$parallel <- TRUE
opts$size <- FALSE #100 # FALSE or integer to adjust size of regions to
opts$peak_or_motif <- "peak" # should the output annotation use the peaks or the motifs?

extend <- round(opts$size/2)

in_files <- dir(io$anno_dir, pattern = ".bed", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] 

filter <- fread(io$filter_bed_file, colClasses = list("factor" = 1L)) %>%
  .[, .(chr, start, end)] %>%
  setkey()

if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}

bed <- future_map(in_files, ~{
  dt <- fread(., colClasses = list("factor" = 1L)) 
  if (opts$size != FALSE){
    dt <- dt[, len := end - start] %>%
      .[, mid := round(start + len/2)] %>%
      .[len != opts$size, c("start", "end") := .(mid - extend, mid + extend)]
  }
  
  setkey(dt, chr, start, end)
  filt <- foverlaps(filter, dt, nomatch = 0L)
  
  if (opts$peak_or_motif == "motif") {
    filt <- filt[, .(chr, start, end, strand, id, anno)]
  } else {
    filt <- filt[, .(chr, start = i.start, end = i.end, strand, id, anno)]
  }
  filt
}) %>%
  rbindlist() %>%
  split(by = "anno")


 

out_files <- paste0(io$anno_out, "/", names(bed), ".bed")

walk2(bed, out_files, fwrite, sep = "\t")
