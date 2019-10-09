library(data.table)
library(dplyr)
library(purrr)
library(furrr)


# this script converts raw bisulfite files of format 
# chr<>pos<>met_reads<>nonmet_reads into bismark cov format 
# chr<>start<>end<>rate<>met_reads<>nonmet_reads


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$data_in <- paste0(io$data_dir, "/acc/raw/pseudobulk")
io$data_out <- paste0(io$data_dir, "/acc/raw/pseudobulk")

opts <- list()
opts$files_to_convert <- "all" # "all" or vector of file names
opts$delete_old_file <- FALSE
opts$compress_new_file <- TRUE



# find files

if (opts$files_to_convert == "all"){
  files <- dir(io$data_in, pattern = ".tsv", full = TRUE)
} else {
  files <- dir(io$data_in, pattern = ".tsv", full = TRUE, pattern = paste(opts$files_to_convert, collapse = "|"))
}

# convert
#plan(sequential)
plan(multiprocess)

future_map(files, ~{
  print(paste0("processing: ", .))
  if (tools::file_ext(.) == "gz") . <- paste("zcat ", .)
  dt <- fread(.) %>%
    .[, rate := round(100 * met_reads / (met_reads + nonmet_reads))] %>%
    .[, .(chr, pos, pos, rate, met_reads, nonmet_reads)]
  out_file <- paste0(io$data_out, "/", basename(.) , ".cov") %>%
    gsub(".tsv.gz|.tsv", "", .)
  fwrite(dt, out_file, sep = "\t", col.names = FALSE)
  if (opts$compress_new_file) {
    system(paste0("gzip -f ", out_file))
  }
})

if (opts$delete_old_file) {
  file.remove(files)
}


