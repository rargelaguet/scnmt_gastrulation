library(purrr)
library(data.table)
library(rtracklayer)
library(R.utils)




#functions

fread_gz <- function(path, ...){
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  if (ext == "gzfile") {
    return(fread(cmd = paste("zcat", path), ...))
  }
  fread(path, ...)  
}


fwrite_tsv <- partial(fwrite, sep = "\t")

#io


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$liftover_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz"
io$tsv_files <- paste0(io$data_dir, "/public_data/Dai_2016/raw")



tmp_file <- paste0(io$tsv_files, "/tmp/liftover.chain.gz")
chain_file <- gsub(".gz", "", tmp_file)
dir.create(dirname(tmp_file), recursive = TRUE)

download.file(io$liftover_url, tmp_file)
gunzip(tmp_file, chain_file)
chain = import.chain(chain_file)
file.remove(chain_file)

tsv_files <- dir(io$tsv_files, pattern = ".tsv.gz", full = TRUE)

walk(tsv_files, ~{
  out_file <- paste0(dirname(.), "/mm10/", basename(.)) %>%
    gsub(".gz", "", .)
  
  dir <- dirname(out_file)
  if (!file.exists(dir)) dir.create(dir, recursive = TRUE)
  
  fread_gz(.) %>%
    .[, .(chr = paste0("chr", chr), start = pos, end = pos, met_reads, nonmet_reads, rate)] %>% 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE, ignore.strand=TRUE) %>% 
    liftOver(chain) %>% 
    as.data.frame() %>% 
    setDT() %>%
    .[, .(chr = gsub("chr", "", seqnames), pos = start, met_reads, nonmet_reads, rate)] %>%
    fwrite_tsv(out_file)
  system(paste("gzip -f", out_file))
})


