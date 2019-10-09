library(data.table)
library(purrr)
library(furrr)

fread_gz <- function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") {
    filename <- paste("zcat", filename)
    return(fread(cmd = filename, ...))
  }
  fread(filename, ...)
}

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$indir <- paste0(io$data_dir, "/features/motifs/bed/1e-4/200bp")
io$meta <-  paste0(io$data_dir, "/sample_metadata.txt")
io$raw_acc <- paste0(io$data_dir, "/acc/raw")

opts <- list()
opts$lineage <- c("Ectoderm", "Mesoderm", "Endoderm")
opts$anno_regex <- "intersect12"
opts$tfs <- c("Sox17", "POU5F1", "T", "TWIST1", "Foxa2", "Gata4", "GATA6", "Pou5f1::Sox2", "HES7"  ,"Sox1", "FOXA1" , "EOMES", "Foxj3", "Pax6", "FOXC1"  )
opts$adjust_mean_rate <- FALSE

anno <- dir(io$indir, pattern = "bed$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(fread) %>%
  rbindlist() %>%
  .[tf %in% opts$tfs] %>%
  setkey(chr, start, end)

meta <- fread(io$meta) %>%
  .[pass_accQC == TRUE & pass_rnaQC == TRUE & lineage %in% opts$lineage] %>%
  .[, stage_lineage := paste0(stage, "_", lineage)]

cells <- meta[, sample]

files <- paste0(io$raw_acc, "/", cells, ".tsv.gz") %>%
  .[file.exists(.)]

motif_acc <- map(files, ~{
  
  cell <- basename(.) %>%
    gsub(".tsv.gz", "", .)
  
  dt <- fread_gz(.) %>%
    .[, .(chr, start = pos, end = pos, rate)] %>%
    setkey(chr, start, end)
  
  if (opts$adjust_mean_rate) mean_rate <- dt[, mean(rate)]
  
  dt <- foverlaps(dt, anno, nomatch = 0L) %>%
    .[, .(rate = mean(rate), .N, sample = cell), tf] %>%
    merge(meta[, .(sample, stage_lineage)], by = "sample")
  
  if (opts$adjust_mean_rate) dt[, rate := rate/mean_rate]
  
  dt
    
}) %>%
  rbindlist()


ggplot(motif_acc, aes(stage_lineage, rate)) +
  geom_boxplot(aes(fill = tf)) +
  facet_wrap(~tf) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
