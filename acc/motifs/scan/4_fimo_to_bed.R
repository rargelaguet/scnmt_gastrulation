###########################################################################
## Script to convert the output from fimo motif scanning into a bed file ##
###########################################################################

library(data.table)
library(purrr)


## Define I/O ##
io           <- list()
io$data_dir  <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$fimo_dir  <- paste0(io$data_dir, "/acc/differential/feature_level/motif_scan/")
io$pwm_file  <- paste0(io$data_dir, "/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt")
io$anno_dir  <- paste0(io$data_dir, "/features/filt")
io$out_dir   <- paste0(io$data_dir, "/features/motifs")

## Define options ##
opts             <- list()
opts$fimo_regex  <- "intersect12"
opts$extend_win  <- TRUE # should motif positions be extended to include extra bp?
opts$win_size    <- 50
opts$p_threshold <- 1e-5 # note that fimo results are already filtered on p<1e-4

## Functions ##
fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA") # save as .tsv


#########################
## Load and parse data ##
#########################

anno <- dir(io$anno_dir, pattern = "bed$", full = TRUE) %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

motif_names <- readLines(io$pwm_file) %>%
  .[grep("MOTIF", .)] %>%
  as.data.table() %>%
  tidyr::separate(".", c("null", "motif", "tf"), sep = " ") %>%
  .[, .(motif, tf)]

motif_anno <- dir(io$fimo_dir, pattern = opts$fimo_regex, full = TRUE) %>%
  dir(pattern = "fimo.txt", full = TRUE) %>%
  set_names(basename(dirname(.)) %>% gsub(".fa", "", .)) %>%
  map(fread, select = c(1:4, 7)) %>%
  map(setnames, c("motif", "id", "min", "max", "p")) %>%
  map2(names(.), ~.x[, anno := gsub(".bed", "", .y) %>% gsub("E7_5", "E7.5", .)]) %>%
  rbindlist() %>%
  .[p < opts$p_threshold] %>%
  merge(anno, by = c("id", "anno")) %>%
  .[, c("start", "end") := .(start + min, start + max)] %>%
  split(by = "motif") %>%
  map(setkey, id) %>%
  map(~.[, id_motif := make.unique(id, sep = "_")]) %>%
  rbindlist() %>%
  setkey(motif) %>%
  merge(motif_names, by = "motif") %>%
  .[, .(chr, start, end, id_motif = paste0(tf, "_", id_motif), anno, tf, id, fimo_p = p)] 

if (opts$extend_win) {
  motif_anno <- motif_anno[, mid := start + (start-end)/2] %>%
    .[, c("start", "end") := .(mid - opts$win_size/2, mid + opts$win_size/2)]
}
  
motif_anno <- split(motif_anno, by = "anno", keep.by = FALSE)
  
  
#################
## Save output ##
#################

files <- paste0(io$out_dir, "/motifs_in_", names(motif_anno), ".bed")

if (opts$extend_win){
  files <- gsub(".bed$", paste0("_", opts$win_size, "bp.bed"), files) %>%
    gsub("/motifs_in_", paste0("/", opts$win_size, "bp/motifs_in_"), .)
}

dir.create(dirname(files[1]), recursive = TRUE)

walk2(motif_anno, files, fwrite_tsv)

