##############################################################
## Script to subset motif PWM file to only include sig hits ##
##############################################################

library(data.table)
library(purrr)


## Define I/O ##
io             <- list()
io$data_dir    <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$pwm         <- paste0(io$data_dir, "/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt")
io$ame_dir     <- paste0(io$data_dir, "/acc/results/differential/motif_enrichment/background_all_k27ac/")
io$out_dir     <- paste0(io$data_dir, "/acc/results/differential/motif_scan")


## Define options ##
opts <- list()
opts$ame_regex <- "" # optionally select ame file with this pattern



#########################
## Load and parse data ##
#########################

# Load position weight matrices (PWM)
pwm <- readLines(io$pwm)

# Load significantly enriched motifs from ame output
ame <- list.dirs(io$ame_dir) %>%
  dir(pattern = "ame.txt", full = TRUE) %>%
  .[!grepl("/background/", .)] %>%
  .[grep(opts$ame_regex, .)] %>%
  set_names(dirname(.) %>% basename) %>%
  map(fread, select = c(8, 9, 13, 16)) %>%
#   map(fread) %>%
  map2(names(.), ~.x[, comparision := .y]) %>%
  rbindlist() %>%
  setnames(c("motif", "tf", "p", "q", "comparison")) %>%
  .[, q := as.numeric(gsub(")", "", q))]


# Filter for significant hits
motifs <- ame[q < 0.05, unique(motif)] 

# Subset PWMs to only significant motifs
sub <- map(motifs, ~{
  i <- grep(., pwm)
  c(pwm[i[1]:i[2]], "")
}) %>%
  unlist() %>%
  c(pwm[1:9], .)

#################
## Save output ##
#################

dir.create(io$out_dir, recursive = TRUE)
writeLines(sub, paste0(io$out_dir, "/pwm_subset.txt"))
writeLines(pwm, paste0(io$out_dir, "/pwm_all.txt"))





