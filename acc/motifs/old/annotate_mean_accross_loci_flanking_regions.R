library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(doParallel)

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_folder <- str_c(io$data_dir,"/features/sjc/motifs/motifdb/filt/motif_only")
io$sample_metadata <- str_c(io$data_dir,"/sample_metadata.txt")
io$in_folder <- str_c(io$data_dir,"/acc/raw/")
io$out_folder <- str_c(io$data_dir,"/acc/parsed/motifs/mean_accross_loci")



opts <- list()
opts$stages <- c("E4.5", "E5.5", "E6.5", "E7.5")
opts$KO_3b <- "not"
opts$flank_size <- 50

meta <- fread(io$sample_metadata) %>%
  .[method=="scNMT" & stage %in% opts$stages & KO_3b %in% opts$KO_3b & pass_accQC==TRUE ]

cells <- meta[, sample]
files <- paste0(io$in_folder, "/", cells, ".tsv.gz") %>% .[file.exists(.)]

anno <- dir(io$anno_folder, full=TRUE, pattern="bed")[1:2] %>%
  map(fread, select=c(1:3, 6), colClasses=list(character=1L)) %>%
  bind_rows() %>%
  .[, c("chr", "anno") := map(.(chr, anno), as.factor)]
  
anno <- list(copy(anno)[, c("start", "end") := .(start-opts$flank_size, start)],
             copy(anno)[, c("start", "end") := .(end, end+opts$flank_size)]) %>%
  bind_rows() %>%
  setkey(chr, start, end)

registerDoParallel(cores=detectCores())
data <- foreach(i=files) %dopar% {
  cell <- basename(i) %>% gsub(".tsv.gz", "", .)
  paste0("zcat ", i) %>%
    fread(select=c(1:2, 5), colClasses=list(factor=1L)) %>%
    setnames("pos", "start") %>%
    .[, end:=start] %>%
    setkey(chr, start, end) %>%
    foverlaps(anno, nomatch=0L) %>%
    .[, .(rate=mean(rate*100), sample=cell), anno]
}
data <- bind_rows(data)

dir.create(io$out_folder, recursive=TRUE)
file_name <- paste0(io$out_folder, "/flank_", opts$flank_size, "bp.tsv")
fwrite(data, file_name)
system(paste0("gzip -f ", file_name))







