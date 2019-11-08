library(data.table)
library(purrr)

# script to download PBAT BSseq from Dai 2016 (TET TKO embryos)
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

convert_to_tsv <- function(dt) {
  setkey(dt, Type) %>%
    .["CpG"] %>%
    .[, .(chr = gsub("chr", "", `#Chr`), 
          pos = Pos, 
          met_reads = Met, 
          nonmet_reads = UnMet, 
          rate = round(MetRate))]  
}



io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$data_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE76261&format=file"

io$data_out <- paste0(io$data_dir, "/public_data/Dai_2016")



cpg_level_dir <- paste0(io$data_out, "/cpg_level")
dir.create(cpg_level_dir, recursive = TRUE)
tar <- paste0(cpg_level_dir, "/tar.tar")

download.file(io$data_url, destfile = tar)

untar(tar, exdir = paste0(cpg_level_dir, "/tar"))
file.remove(tar)


# convert to tsv format

files <- dir(paste0(cpg_level_dir, "/tar"), pattern = "SingleCmet.txt.gz", full = TRUE)

walk(files, ~{
  outfile <- gsub(".SingleCmet.txt.gz", ".tsv", .) %>% sub("/tar/", "/", .)
  fread_gz(.) %>%
    convert_to_tsv() %>%
    fwrite(outfile)
  system(paste("gzip -f", outfile))
})




