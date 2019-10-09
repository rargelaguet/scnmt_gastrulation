library(data.table)
library(purrr)

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"

io$in.data <- paste0(io$data_dir, "/acc/raw")
io$out_file <- paste0(io$data_dir, "/acc/output/global_rates/global_rates.csv")


files <- dir(io$in.data, pattern=".tsv.gz", full=TRUE)
global_rates <- map(files, ~{
  cell <- basename(.) %>%
    gsub(".tsv.gz", "", .)
  paste("zcat", .) %>% 
    fread() %>%
    .[, .(sample=cell, global_rate=100*mean(rate), .N)]
    
}) %>%
  rbindlist()

fwrite(global_rates, io$out_file)