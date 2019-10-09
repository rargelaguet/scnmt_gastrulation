suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(furrr))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
options(future.globals.maxSize = 1000 * 1024 ^ 2)
options(datatable.showProgress = FALSE)
options(datatable.verbose = FALSE)

# take unbiased windows accross the whole genome and test for sig differences in
# accessibility rates between lineages

## in/out ##
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$sample_meta <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$acc_data <- paste0(io$data_dir, "/acc/parsed/unbiased/e4.5/win100_step50/parsed_stageE4.5") # pre-generated single-cell accessibility data at short overlapping windows. data is split into chunks for memory efficiency (or for parallel computing)
io$acc_anno <- paste0(io$data_dir, "/acc/parsed/unbiased/e4.5/win100_step50/annotation.bed") # annotation used to generate accessibility data (e.g. overlapping windows)
io$out_dir <- paste0(io$acc_data, "/t_tests")

## options ##
opts <- list()
opts$KO_3b <- "not"
opts$min_weight <- 5
opts$min_sample_size_per_group <- 10L # for t-tests
#opts$weighted_t_test <- FALSE # should probably use weighted....
opts$parallel <- TRUE
opts$permutation <- TRUE
opts$n_perm <- 100L

if (length(args > 1)) opts$parallel <- args[2]

opts$stage_lineage1 <-  "E4.5_EPI"  # only 1 lineage
opts$stage_lineage2 <- "E4.5_PE"
# opts$stage_lineage1 <- args[1]
# lineages <- c("E4.5_EPI", "E4.5_PE")
# opts$stage_lineage1 <- lineages[grep(toupper(opts$stage_lineage1), toupper(lineages))]
# opts$stage_lineage2 <- lineages[-grep(toupper(opts$stage_lineage1), toupper(lineages))]


print(paste("testing", opts$stage_lineage1, "versus", paste(opts$stage_lineage2, collapse = " and ")))


# TESTING #
opts$testing <- FALSE

# functions
#source("/bi/home/clarks/fread_gz.R")
fread_gz <- function(filename, ...){
  if (grepl("gz", filename)) return(fread(paste("zcat", filename), ...))
  fread(filename, ...)
}

t_test_dt <- function(x, y, ...){
  #na_return <- list(p=as.numeric(NA), mean_x=as.numeric(NA), mean_y=as.numeric(NA))
  na_return <- list(p=as.numeric(NA), mean_dif=as.numeric(NA))
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x)<opts$min_sample_size_per_group||length(y)<opts$min_sample_size_per_group) return(na_return)
  if( var(c(x, y)) == 0) return(na_return)
  dt <- t.test(x, y, ...)
  list(p = dt$p.value, mean_dif = dt$estimate[1] - dt$estimate[2])
}

find_dif_acc <- function(x, lineage1, other_lineages){
  x[, t_test_dt(get(lineage1), mget(other_lineages) %>% unlist), id] 
} 



## load metadata and find files
stage_lineage <- c(opts$stage_lineage1, opts$stage_lineage2)

meta <- fread(io$sample_meta) %>% 
  .[, lineage := paste(stage, lineage, sep = "_")] %>%
  .[lineage %in% stage_lineage & 
      KO_3b %in% opts$KO_3b & 
      pass_accQC == TRUE ]

lineages <- meta[, unique(lineage)]
cells <- meta[, sample]

chunks <- dir(io$acc_data, full = TRUE, pattern = ".tsv")

# now perform t-tests in chunks

# parallel processing
if (opts$parallel) {
  plan(multiprocess)
} else {
  plan(sequential)
}

if (opts$testing) {
  chunks <- chunks[1:2]
}

dir.create(io$out_dir, recursive = TRUE)

chunks <- rev(chunks)

future_map(chunks, ~{
# dif_acc <- map(chunks, ~{
  meta <- meta[, .(sample, lineage = as.factor(lineage))]
  acc <- fread_gz(.) %>%
    setnames(c("sample", "id", "anno", "rate", "weight")) %>%
    setkey(sample) %>%
    .[sample %in% cells] %>%
    merge(meta, by = "sample") %>%
    setkey(weight) %>%
    .[weight >= opts$min_weight]
  
  if (opts$testing) {
    acc=acc[id %in% sample(id, 100)]
  }
  
  t_test <- dcast(acc, sample + id ~ lineage, value.var = "rate") %>%
    find_dif_acc(opts$stage_lineage1, opts$stage_lineage2) %>%
    .[complete.cases(.)]
  
  if (!opts$permutation) return (t_test)
  
  perms <- map(1:opts$n_perm , ~acc[, lineage := sample(lineage), id] %>%
                 dcast(sample + id ~ lineage, value.var = "rate") %>%
                 find_dif_acc(opts$stage_lineage1, opts$stage_lineage2) %>%
                 .[complete.cases(.), .(id, p)] %>%
                 setkey(id)
    ) %>%
    map2(paste0("perm", c(1:opts$n_perm)), ~setnames(.x, "p", .y)) %>%
    map(~.[, 2]) %>%
    do.call(cbind, .) %>%
    cbind(t_test)
  

# compute proportion of permuted p-values smaller than observed p-value 

    perms[, p_permuted := apply(.SD, 1, function(x) sum(x[1] >= x[2:length(x)])) / opts$n_perm, .SDcol = c("p", paste0("perm", 1:opts$n_perm))] %>%
      .[, .(id, p, mean_dif, p_permuted)]

out_file <- paste0(io$out_dir, "/", basename(.))
fwrite(perms, out_file, sep = "\t")
out_file
}) 

dif_acc <- dir(io$out_dir, full = TRUE) %>%
  map(fread) %>% 
  rbindlist() %>%
  setkey(id)


anno <- fread(io$acc_anno) %>%
  setkey(id)

dif_acc <- setkey(dif_acc, id) %>% anno[.]


# save the result
dir.create(io$out_dir, recursive = TRUE)


lin_name <- paste0(opts$stage_lineage1, "_vs_" , paste(opts$stage_lineage2, collapse = "_"))
if (opts$permutation) perm_name <- paste0("_", opts$n_perm, "perms") else perm_name <- "_"

out_file <- paste0(io$out_dir, "/", lin_name, perm_name, "_t-tests.csv")

if (opts$testing) out_file <- sub("_t-tests.csv", "_t-tests_TESTING.csv", out_file)

fwrite(dif_acc, out_file)



# # old code
# dif_acc <- future_map(chunks, ~{
#   meta <- meta[, .(sample, lineage = as.factor(lineage))]
#   acc <- fread_gz(.) %>%
#     setnames(c("sample", "id", "anno", "rate", "weight")) %>%
#     setkey(sample) %>%
#     .[sample %in% cells] %>%
#     merge(meta, by = "sample") %>%
#     setkey(weight) %>%
#     .[weight >= opts$min_weight]
#   
#   if (opts$testing) {
#     acc=acc[id %in% sample(id, 100)]
#   }
#   
#   
#   if (opts$permutation){
#     perms <- c(paste0("perm", 1:opts$n_perm))
#     # generate permutations by shuffling lineage
#     acc <- acc[, c(perms) := map(1:opts$n_perm, ~sample(lineage)), id]    
#   } else {
#     perms <- as.character(NULL)
#   }
#   # make a list of dcast dt's - one for observed, and the rest perumations
#   acc <- map(c("lineage", perms), ~dcast(acc, sample + id ~ get(.), value.var = "rate")) %>%
#     setNames(c("observed", perms))
#   
#   # t-tests on each locus, and for each dataset (observed and permuations)
#   map(acc, find_dif_acc, opts$stage_lineage1, opts$stage_lineage2) %>%
#     map(~.[complete.cases(.)]) %>%
#     # make into 1 data.table
#     map2(names(.), ~setnames(.x, c("p", "mean_dif"), paste0(c("p_", "mean_dif_"), .y))) %>%
#     map(setkey, id) %>%
#     purrr::reduce(merge) %>%
#     # compute proportion of permuted p-values smaller than observed p-value   
#     .[, p_permuted := apply(.SD, 1, function(x) sum(x[1] >= x[2:length(x)])) / opts$n_perm, .SDcol = c("p_observed", paste0("p_perm", 1:opts$n_perm))] %>%
#     .[, .(id, p_observed, mean_dif_observed, p_permuted)]
#   
# }) %>%
#   rbindlist() %>%
#   setkey(id)

