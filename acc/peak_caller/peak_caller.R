##################################################################
## Script to define open peaks based on chromatin accessibility ##
##################################################################

# Strategy: 

library(data.table)
library(purrr)
library(furrr)

options(future.globals.maxSize = 1000 * 1024 ^ 2)

#########
## I/O ##
#########

io <- list()
io$basedir <- "/Users/ricard/data/gastrulation"
# io$datadir <- paste0(io$basedir,"/acc/gpc_level/pseudobulk")
io$datadir <- paste0(io$basedir,"/met/cpg_level/pseudobulk")
# io$acc_pseudobulk <- paste0(io$data_dir, "/acc/raw/gpc_level/pseudobulk/E7.5_all.tsv.gz")
io$outdir <- paste0(io$basedir,"/met/cpg_level/pseudobulk")

#############
## Options ##
#############

opts <- list()

# Window settings
# opts$window <- 100  # window size
# opts$step <- 50     # step size
opts$window <- 1000
opts$step <- 250

# Minimum
opts$min_coverage <- 100

# merge adjacent windows if they have similar rates
opts$merge_overlapping_windows <- TRUE 

# maximum rate difference between adjacent sites to merge sites
opts$max_diff <- 10 

# Select top sites (in percentage)
opts$cutoff <- 0.05

# pick the most ("top") or least ("bottom") accessible sites
opts$top_or_bottom <- "top" 

# Define which stage and lineages to look at 
opts$stage_lineage <- c(
  
  # E4.5
  "E4.5_Epiblast"
  
  # # E5.5
  # "E5.5_Epiblast",
  # 
  # # E6.5
  # "E6.5_Epiblast",
  # "E6.5_Primitive_Streak",
  # "E6.5_Mesoderm",
  # 
  # # E7.5
  # "E7.5_Epiblast",
  # "E7.5_Primitive_Streak",
  # "E7.5_Ectoderm",
  # "E7.5_Endoderm",
  # "E7.5_Mesoderm"
)

###############
## Load data ##
###############

data <- fread(sprintf("%s/%s.tsv.gz",io$datadir,opts$stage_lineage)) %>% 
  .[complete.cases(.)]

############################
## Define running windows ##
############################

windows <- data[, .(start = seq(min(pos), max(pos), by = opts$step)), chr] %>%
  .[, end := start + opts$window] %>%
  setkey(chr, start, end)

data <- data %>% setnames("pos", "start") %>%
  .[, end := start] %>%
  setkey(chr, start, end) 

##############################################
## Summarise accessibility over the windows ##
##############################################

# split by chr to process as chunks
data <- split(data, by = "chr")
windows <- split(windows, by = "chr")

plan(multiprocess)
data <- future_map2(data, windows, foverlaps, nomatch=0L) %>% rbindlist %>%
  .[, .(met_cpgs = sum(met_cpgs), nonmet_cpgs = sum(nonmet_cpgs)), .(chr, start, end)] %>%
  .[, total_counts := met_cpgs+nonmet_cpgs]
  .[, rate := round(100*met_cpgs/total_counts)]

#################################
## Combine overlapping windows ##
#################################

if (opts$merge_overlapping_windows){

  data <- data[order(chr, start)] %>%
    .[, dif_r := rate - shift(rate, type = "lag")] %>%
    .[is.na(dif_r), dif_r := 100] %>%
    .[, dif_s := shift(start, type = "lead") - start] %>%
    .[, not_overlap := TRUE] %>%
    .[dif_s < opts$window, not_overlap := FALSE]  %>% 
    .[, group := as.numeric(abs(dif_r) > opts$max_diff | not_overlap) %>% cumsum] %>%
    .[, .(start = min(start), end = max(end), met_cpgs = sum(met_cpgs), total_counts = sum(total_counts)), .(chr, group)] %>%
    .[, rate := round(100*met_cpgs/total_counts)]
}


##################
## Filter sites ##
##################

anno <- data[total_counts>=opts$min_coverage] 

######################
## Select top sites ##
######################

if (opts$top_or_bottom == "bottom") {
  anno <- anno[order(rate)]
} else {
  anno <- anno[order(-rank(rate))]
}
anno <- anno[1:(.N*opts$cutoff)]

##########
## Save ##
##########

io$outfile <- sprintf("%s/%s_window%s_step%s_cutoff%s_mincov%s.tsv.gz",io$outdir,paste0(opts$stage_lineage,collapse="_"), opts$window, opts$step, opts$cutoff,opts$min_coverage)
if (opts$merge_overlapping_windows) io$outfile <- gsub(".tsv.gz", "_merged.tsv.gz", io$outfile)

# re-format and save

# anno_name <- basename(out_file) %>% gsub(".bed", "", .) %>% paste0("peaks_", .)
# anno <- anno[, c("strand", "id", "anno") := .("*", paste0(anno_name, "_", .I), anno_name)]
# anno <- anno[, .(chr, start, end, strand, id, anno)]

dir.create(io$outdir, recursive=TRUE, showWarnings = F)
fwrite(anno, io$outfile, sep="\t")


