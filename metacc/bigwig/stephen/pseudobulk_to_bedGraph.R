library(data.table)
library(purrr)
#library(furrr)

# reads in nmt-seq data (optionally split by celltype) then pseudobulks and saves as a bedGraph file
# optionally smoothes the data by quantifying over sliding windows

io <- list()


io$nmt_meta    <- "/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata.txt"
io$nmt_acc     <- "/bi/scratch/Stephen_Clark/gastrulation_data/acc/raw/scNMT/"
io$nmt_met     <- "/bi/scratch/Stephen_Clark/gastrulation_data/met/raw/scNMT/"

io$chrom_sizes <- "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"

io$outdir      <- "/bi/scratch/Stephen_Clark/gastrulation_data/bedGraph/smooth200_step50"

dir.create(io$outdir, recursive = TRUE)

opts <- list()
opts$split_by        <- "lineage10x_2"
opts$transform_rates <- function(x){x * 100 + 10} # function to modify rates to enable better viewing in IGV
opts$smooth          <- TRUE # smooth data by quantifying over sliding windows
opts$win             <- 200
opts$step            <- 50 



# function to read in data one file at a time (in case of memory issues, otherwise just load all then do the sum)
fread_and_sum <- function(dt, file){
  dt2 <- fread(file)
  
  rbindlist(list(dt, dt2)) %>% 
    .[, .(met = sum(met_reads), non = sum(nonmet_reads)), .(chr, pos)]
    
}

# ncores <- parallel::detectCores()
# ncores
# if (ncores > 1) plan(multisession, workers = ncores)


meta <- fread(io$nmt_meta) %>% 
  .[, stage_lineage := paste0(stage, "_", get(opts$split_by))] %>% 
  .[! stage_lineage %like% "NA"] %>% 
  .[pass_accQC == TRUE & pass_metQC == TRUE] %>% 
  split(by = "stage_lineage")

names(meta)



meta <- meta[names(meta) %like% "E4.5|E5.5|E6.5"]
names(meta)
#meta <- meta[1]


if (opts$smooth){
  anno <- fread(io$chrom_sizes) %>% 
    setnames(c("chr", "max")) %>% 
    .[chr %in% paste0("chr", c(1:19, "X", "Y", "M"))] %>% 
    .[, .(start = seq(2, max - opts$win, by = opts$step)), chr] %>%  # start at 2 as ucsc_tools is getting confused with multiples of 10
    .[, end := start +opts$win - 1] %>% 
    setkey(chr, start, end)
    
}

files <- map2(meta, names(meta), ~{
  cells <- paste0(.x[, id_acc], ".tsv.gz")
  
  files <- list(met = io$nmt_met, acc = io$nmt_acc) %>% 
    map(dir, recursive = TRUE, pattern= ".tsv", full = TRUE) %>% 
    map(~.[basename(.) %in% cells])
  
  
  met <- map(files$met, fread, select = c(1:2,5)) %>% 
    rbindlist() %>% 
    .[, .(rate = mean(rate)), .(chr, pos)] %>% 
    .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>% 
    setcolorder(c("chr", "start", "end", "rate")) %>% 
    setkey(chr, start, end)
  
  if (met[1, !grepl("chr", chr)]) {
    met[, chr := paste0("chr", chr)]
    
  }
  
   met[, rate := opts$transform_rates(rate)]
  
  met[chr=="chrMT", chr := "chrM"]
  
  met_outfile <- paste0(io$outdir, "/CpG_met_", .y, ".bedGraph")
  
  if (opts$smooth){
    setkey(met, chr, start, end)
    met <- foverlaps(met, anno, nomatch = 0L) %>% 
      .[, .(rate = round(mean(rate))), .(chr, start ,end)] %>% 
      .[, end := start + opts$step-1]
    
    
    
    met_outfile <- gsub(".bedGraph$", "_smooth.bedGraph", met_outfile)
    
  }
  
  setkey(met, chr, start, end)
  
  
  print(paste("saving", met_outfile))
  fwrite(met, met_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  rm(met)
  
  
  acc <- map(files$acc, fread, select = c(1:2,5)) %>% 
    rbindlist() %>% 
    .[, .(rate = mean(rate)), .(chr, pos)] %>% 
    .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>% 
    setcolorder(c("chr", "start", "end", "rate")) %>% 
    setkey(chr, start, end)
  
  if (acc[1, !grepl("chr", chr)]) {
    acc[, chr := paste0("chr", chr)]
    
  }
  
   acc[, rate := opts$transform_rates(rate)]
  
  acc[chr=="chrMT", chr := "chrM"]
  
  acc_outfile <- paste0(io$outdir, "/GpC_acc_", .y, ".bedGraph")
  
  if (opts$smooth){
    setkey(acc, chr, start, end)
    acc <- foverlaps(acc, anno, nomatch = 0L) %>% 
      .[, .(rate = round(mean(rate))), .(chr, start ,end)]%>% 
      .[, end := start + opts$step-1]
    
    acc_outfile <- gsub(".bedGraph$", "_smooth.bedGraph", acc_outfile)
    
  }
  
  
  
  setkey(acc, chr, start, end)
  
  
  print(paste("saving", acc_outfile))
  fwrite(acc, acc_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  rm(acc)
  
  c(met_outfile, acc_outfile)
})


