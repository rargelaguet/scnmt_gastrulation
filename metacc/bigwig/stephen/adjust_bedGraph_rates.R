library(data.table)
library(purrr)


# normalises nmt-seq rates to genome-average

io <- list()



io$indir      <- "/bi/scratch/Stephen_Clark/gastrulation_data/bedGraph/smooth200_step50"
io$outdir     <- file.path(io$indir, "adjusted")

dir.create(io$outdir)

opts <- list()
opts$transform_rates <- function(x){x - 10} # reverse the function that was used to adjust the rates when generating the bdgraph

infiles <- dir(io$indir, pattern = ".bedGraph", full = TRUE) %>% .[. %like% "GpC"]
outfiles <- file.path(io$outdir, basename(infiles))
.x=files[1]

walk2(infiles, outfiles, ~{
  dt <- fread(.x) %>% 
    .[, V4 := opts$transform_rates(V4)]
  mean <- dt[, mean(V4)]
  dt[, V4 := V4/mean + 0.1]
  
  fwrite(dt, .y, sep = "\t", col.names = FALSE, quote = FALSE)
})





