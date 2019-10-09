library(pheatmap)
library(purrr)
library(data.table)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$meta_data <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$parsed_data_dir <- paste0(io$data_dir, "/acc/parsed/motifdb_100bp/filt_by_peaks")

opts <- list()
opts$stage <- c("E5.5", "E6.5", "E7.5") #"all"#"E7.5" #"all"#"E7.5" #"all"
opts$lineage <- c("EPI", "PE", "PS", "Ectoderm", "Mesoderm", "Endoderm")#"all"
opts$KO_3b <- "not"
opts$anno_regex <- ""
opts$min_weight <- 100
opts$min_cells <- 10 # min cells per group
opts$mean_by_cell_type <- FALSE
opts$var_cutoff <- 20 # top x% to be kept

# fun
fread_gz = function(filename, ...){
  if (grepl(".gz", filename)) filename <- paste("zcat", filename)
  fread(filename, ...)
}

meta <- fread(io$meta_data) %>%
  .[stage=="E6.75", stage := "E6.5"] %>%
  .[pass_accQC == TRUE]

if (opts$stage[1] != "all") meta <- meta[stage %in% opts$stage]
if (opts$lineage[1] != "all") meta <- meta[lineage %in% opts$lineage]
if (opts$KO_3b[1] != "all") meta <- meta[KO_3b %in% opts$KO_3b]

cells <- meta[, sample]


acc_data <- dir(io$parsed_data_dir, pattern = ".tsv.gz", full = TRUE) %>%
  .[grep(opts$anno_regex, .)] %>%
  map(~{
    dt=fread_gz(.) %>%
      .[sample %in% cells] %>%
      .[N >= opts$min_weight]
  }) %>%
  rbindlist() 

if (opts$mean_by_cell_type) {
  acc_data <- merge(acc_data, meta[, .(sample, type = paste0(stage, "_", lineage))], by = "sample") %>%
    .[, .(rate = mean(rate), sum_n = sum(N), .N), .(type, anno)] %>%
    .[N >= opts$min_cells] %>%
    setnames("type", "sample")
}

top_var <- acc_data[, .(var = var(rate), .N), anno] %>% 
  .[order(-rank(var))] %>% 
  .[!is.na(var)] %>% 
  .[N/length(cells) >= 0.9] %>%
  .[1:(.N * opts$var_cutoff/100), anno]

acc_mat <- dcast(acc_data[anno %in% top_var], anno~sample, value.var = "rate") %>%
  setDF() %>%
  tibble::column_to_rownames("anno") %>%
  as.matrix()

# trim names of annotations
rownames(acc_mat) <- gsub("HOCOMOCOv10-|JASPAR_2014-|JASPAR_CORE-|jaspar2016-|jolma2013-|UniPROBE-", "", rownames(acc_mat)) %>%
  strsplit(split = "_|-|\\.") %>%
  map_chr(1)


acc_mat <- acc_mat[unique(rownames(acc_mat)), ]


anno_col <- meta[, .(sample, stage, lineage)] %>%
  .[sample %in% colnames(acc_mat)] %>%
  .[order(stage, lineage)] %>%
  setDF() %>%
  tibble::column_to_rownames("sample")

acc_mat <- acc_mat[, rownames(anno_col)]



if (opts$mean_by_cell_type){
  pheatmap(acc_mat)
} else {
  pheatmap(acc_mat, annotation_col = anno_col, show_colnames = FALSE)
}


#


