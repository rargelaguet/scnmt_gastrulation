##################################################################################
## Script to generate fasta sequence files from differentially accessible sites ##
##################################################################################


library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)

## Define I/O ##
io            <- list()
io$data_dir   <- "/bi/scratch/Stephen_Clark/gastrulation_data"                # base directory       
io$in_dir     <- paste0(io$data_dir, "/acc/differential/feature_level/")      # results of differential accessibility
io$anno       <- paste0(io$data_dir, "/features/filt")                        # feature bed files
io$out_dir    <- paste0(io$data_dir, "/acc/differential/feature_level/fasta") # where to save output

## Define options ##
opts            <- list()
opts$regex      <- "intersect12" # pattern to select feature (bed) files
opts$up_or_down <- "up"          # "up" "down" or "both"
opts$min_sites  <- 100           # discard any sets of features with less than X rows


## Functions ##

# Read gz compressed file
fread_gz <- function(path, ...){
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  if (ext == "gzfile") {
    return(fread(cmd = paste("zcat", path), ...))
  }
  fread(path, ...)  
}

# Extract sequences
get_seqs <- function(bed_dt, BSgenome){
  seq <- bed_dt[, getSeq(mm10, 
                         names = paste0("chr", chr), 
                         start = start, 
                         end = end)]
  seq@ranges@NAMES <- bed_dt[, as.character(id)]
  seq  
}

################################################################################

###############
## Load data ##
###############

# Load feature annoations
anno <- dir(io$anno, full = TRUE, pattern = ".bed$") %>%
  .[grep(opts$regex, .)] %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

# Load differential results
dif_results <- dir(io$in_dir, full = TRUE, pattern = ".txt") %>%
  .[grep(opts$regex, .)] %>%
  set_names(basename(.) %>% gsub(".txt.gz|.txt", "", .)) %>%
  map(fread_gz) 

# Insert name of annotation
names(dif_results) <- map2(dif_results, names(dif_results), ~{
  anno <- .x[, unique(anno)]
  nm <- .y
  if (!grepl(anno, .y)) nm <- paste(.y, anno, sep = "_")
  nm
})

# Load genome sequence
mm10 <- BSgenome.Mmusculus.UCSC.mm10

###########################
## Parse and filter data ##
###########################

# Select only significant hits
dif_bed <- map(dif_results, ~.[sig == TRUE]) 

# Select up or down accessibility (else include both)
if (toupper(opts$up_or_down) == "UP"){
  dif_bed <- map(dif_bed, ~.[diff < 0])
} else if (toupper(opts$up_or_down) == "DOWN") {
  dif_bed <- map(dif_bed, ~.[diff > 0])
}

# Filter out low covered features and combine with annotation file for genomic coordinates
dif_bed <- purrr::discard(dif_bed, ~nrow(.) < opts$min_sites) %>%
  map(~.[, .(anno, id)]) %>%
  map(merge, anno, by = c("anno", "id")) %>%
  map(~.[, .(chr, start, end, strand, anno, id)])


####################
## Find sequences ##
####################

seqs <- map(dif_bed, get_seqs, mm10)

# Also find sequences of background (i.e. all K27ac sites)
background_seqs <- split(anno, by = "anno") %>%
  map(get_seqs, mm10)

##################
## Save results ##
##################

# Generate directory and file names
bed_anno_name <- map_chr(dif_bed, ~.[, unique(anno)]) 
out_dirs <- paste0(io$out_dir, "/", bed_anno_name, "/", opts$up_or_down)
out_files <- paste0(out_dirs, "/", names(seqs), ".fa")
# Generate background file names
b_out_files <- paste0(io$out_dir, "/", names(background_seqs), "/background/", names(background_seqs), ".fa")

# Save fasta files
walk(out_dirs, dir.create, recursive = TRUE)
walk2(seqs, out_files, writeXStringSet)

# Save background fasta file
walk(b_bout_files, ~dir.create(dirname(.), recursive = TRUE))
walk2(background_seqs, b_out_files, writeXStringSet)
