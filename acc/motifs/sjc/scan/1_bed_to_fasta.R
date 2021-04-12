############################################################
## Script to generate fasta sequence files from bed files ##
############################################################


library(data.table)
library(purrr)
library(BSgenome.Mmusculus.UCSC.mm10)

## Define I/O ##
io            <- list()
io$data_dir   <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$in_dir     <- paste0(io$data_dir, "/features/filt")
io$out_dir    <- paste0(io$data_dir, "/features/fasta")


###############
## Load data ##
###############

# Feature files
bed <- dir(io$in_dir, pattern = ".bed$", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(~fread(.)[, file := basename(.)]) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno", "file")) %>%
  split(by = c("file"), keep.by = FALSE)

# Genome sequence
mm10 <- BSgenome.Mmusculus.UCSC.mm10

####################
## Find sequences ##
####################

seqs <- map(bed, ~{
  seq <- .[, getSeq(mm10, 
                    names = paste0("chr", chr), 
                    start = start, 
                    end = end)]
  seq@ranges@NAMES <- .[, as.character(id)]
  seq
})

##################
## Save results ##
##################

dir.create(io$out_dir, recursive = TRUE)
out_files <- paste0(io$out_dir, "/", names(seqs) %>% sub(".bed", ".fa", .))
walk2(seqs, out_files, writeXStringSet)
