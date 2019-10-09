library(data.table)
library(purrr)
library(furrr)
library(ggplot2)
library(cowplot)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data/"
io$meta <- paste0(io$data_dir, "/sample_metadata.txt")
io$raw_acc <- paste0(io$data_dir, "/acc/raw")
io$fimo <- paste0(io$data_dir, "/features/fimo")
io$anno <- paste0(io$data_dir, "/features/filt")
io$pwm_file <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt"

opts <- list()
opts$stage <- c("E7.5", "E6.5")
opts$lineage <- c("EPI", "Ectoderm")#c("Mesoderm", "Ectoderm", "Endoderm")
opts$anno_regex <- "intersect12"
opts$motifs <- c("Pou5f1::Sox2", "Foxa2", "TWIST1", "Gata4", "MEIS1", "CTCF", "Sox17", "Gata1", "POU3F1")

opts$motifs <- c("MEIS2", "FOXB1" )

meta <- fread(io$meta) %>%
  .[stage %in% opts$stage & lineage %in% opts$lineage & KO_3b == "not"]

cells <- split(meta, by = "lineage") %>%
  map(~.[sample(1:.N, 100)]) %>%
  map(~.[, sample]) 

(files <- map(cells, ~paste0(io$raw_acc, "/", ., ".tsv.gz") %>% .[file.exists(.)]))

motif_names <- readLines(io$pwm_file) %>%
  .[grep("MOTIF", .)] %>%
  as.data.table() %>%
  tidyr::separate(".", c("null", "motif", "tf"), sep = " ") %>%
  .[, .(motif, tf)]

anno <- dir(io$anno, pattern = opts$anno_regex, full = TRUE) %>%
  map(fread) %>%
  rbindlist() %>%
  setnames(c("chr", "start_pos", "end_pos", "strand", "id", "anno"))

fimo <- dir(io$fimo, full = TRUE, pattern = opts$anno_regex) %>%
  set_names(basename(.) %>% gsub(".bed.fa", "", .)) %>%
  dir(pattern = "fimo.txt$", full = TRUE) %>%
  map(fread) %>%
  map(~.[`p-value`<0.05, .(motif = `#pattern name`, start, end = `stop`, id = `sequence name`)]) %>%
  rbindlist() %>%
  merge(motif_names, by = "motif") %>%
  merge(anno[, .(id, chr, start_pos)], by = "id") %>%
  .[, .(id, tf, chr, start = start_pos + start, end = start_pos + end)] %>%
  .[, mid := start + round((end-start)/2)] %>%
  .[, c("start", "end") := .(start - 1000, end + 1000)] 

fimo_select <- fimo[tf %in% opts$motifs]

win <- data.table(start = seq(-1000, 1000, by = 25)) %>%
  .[, end := start + 100] %>%
  setkey()

dt <- map2(files, names(files), ~{
  plan(multiprocess)
  d=future_map(.x, ~fread(.) %>%
    .[, .(chr, start = pos, end = pos, rate)] %>%
    setkey(chr, start, end) %>%
    foverlaps(fimo_select %>% setkey(chr, start, end), nomatch = 0L) %>%
    .[, c("start", "end") := .(mid - i.start, mid - i.start)] %>%
      .[, c("i.start", "i.end") := NULL] %>%
      setkey(start, end) %>%
      foverlaps(win %>% setkey(start, end), nomatch = 0L)) %>%
    rbindlist() %>%
    .[, .(rate = mean(rate), sample = .y), .(tf, start)]
}) %>%
  rbindlist()

sel <- c("CTCF", "Foxa2", "Gata1", "Gata4", "Pou5f1::Sox2", "TWIST1", "Sox17", "POU3F1" )
sel=opts$motifs

ggplot(dt[tf %in% sel], aes(start, rate, colour = sample)) +
  geom_line() +
  facet_wrap(~tf, nrow=2) +
  xlab("bp from motif") +
  ylab("mean accessibility in 100bp windows")





