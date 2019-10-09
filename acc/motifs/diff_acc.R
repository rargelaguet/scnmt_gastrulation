library(data.table)
library(purrr)
library(dplyr)

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$meta_data <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$parsed_data_dir <- paste0(io$data_dir, "/acc/parsed/motifdb_100bp/all")
io$out <- paste0(io$data_dir, "/acc/output/motif_dif_acc/motifdb_100bp")

if (!dir.exists(io$out)) dir.create(io$out, recursive = TRUE)

opts <- list()
opts$stage_lineage1 <- "E7.5_Endoderm"
opts$stage_lineage2 <- c("E7.5_Ectoderm", "E7.5_Mesoderm")#c("E6.5_PS")
opts$min_weight <- 100
opts$min_cells <- 10 # min cells per group


# fun
fread_gz = function(filename, ...){
  if (grepl(".gz", filename)) filename <- paste("zcat", filename)
  fread(filename, ...)
}

meta <- fread(io$meta_data) %>%
  .[stage=="E6.75", stage := "E6.5"] %>%
  .[pass_accQC == TRUE] %>%
  .[, stage_lineage := paste0(stage, "_", lineage)] %>%
  .[stage_lineage %in% c(opts$stage_lineage1, opts$stage_lineage2)] %>%
  .[, group := case_when(stage_lineage %in% opts$stage_lineage1 ~ "group1",
                         stage_lineage %in% opts$stage_lineage2 ~ "group2")]



cells <- meta[, sample]


acc_data <- dir(io$parsed_data_dir, pattern = ".tsv.gz", full = TRUE) %>%
  map(~{
    fread_gz(.) %>%
      .[sample %in% cells] %>%
      .[N >= opts$min_weight]
  }) %>%
  rbindlist() %>%
  setkey(sample) %>%
  merge(meta[, .(sample, group)] %>% setkey(sample)) %>%
  dcast(sample + anno ~ group, value.var = "rate")

do_t_test <- function(x, y){
  nx <- sum(!is.na(x)) 
  ny <- sum(!is.na(y)) 
  
  ret <- as.double(NA) %>%
    list(p = ., nx = nx, ny = ny, mean_x = ., mean_y = .)
  
  if (nx < opts$min_cells || ny < opts$min_cells) return(ret)
  
  t <- t.test(x, y)
  
  list(p = t$p.value, nx = nx, ny = ny, mean_x = t$estimate[1], mean_y = t$estimate[2])
}

st_lin <- paste0(opts$stage_lineage1, "_vs_", paste(opts$stage_lineage2, collapse = "_"))

t_test <- acc_data[, do_t_test(group1, group2), anno] %>%
  .[, padj := p.adjust(p, method = "fdr")] %>%
  .[, type := st_lin]



out_file <- paste0(io$out, "/", st_lin, ".csv")
fwrite(t_test, out_file)

