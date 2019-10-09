library(data.table)
library(purrr)



io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$pwm <- paste0(io$data_dir, "/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt")
io$ame_dir <- paste0(io$data_dir, "/acc/differential/feature_level/motif_enrichment/background_all_k27ac/")
io$out_dir <- paste0(io$data_dir, "/acc/differential/feature_level/motif_scan/")

dir.create(io$out_dir, recursive = TRUE)

opts <- list()
opts$ame_regex <- ""

pwm <- readLines(io$pwm)

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


motifs <- ame[q < 0.05, unique(motif)] 


sub <- map(motifs, ~{
  i <- grep(., pwm)
  c(pwm[i[1]:i[2]], "")
}) %>%
  unlist() %>%
  c(pwm[1:9], .)

writeLines(sub, paste0(io$out_dir, "/pwm_subset.txt"))
writeLines(pwm, paste0(io$out_dir, "/pwm_all.txt"))





