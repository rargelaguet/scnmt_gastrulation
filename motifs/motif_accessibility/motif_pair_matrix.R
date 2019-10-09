library(data.table)
library(purrr)
library(pheatmap)

# find pairwise combinations of motifs and produce matrix with score for each pair


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$motif_bed <- paste0(io$data_dir,"/features/motifs/020519/200bp")
io$dif_acc <- paste0(io$data_dir, "/acc/differential/feature_level")

io$out_dir <- "/bi/home/clarks/gastrulation/plots//motifs/motif_pair_matrix"

dir.create(io$out_dir, recursive = TRUE)
# io$out <- 

opts<-list()
opts$motfis <- c("SP8", "POU5F1", "POU3F1", "TFAP2A", "SOX2", "SOX3", 
                 "FOXA1", "HNF1B",  "FOXA2", "GATA1", "SOX17",
                 "TWIST1", "GATA4",   "HAND1", "MEIS1", "LEF1")
opts$max_dist <- 1000
opts$min_dist <- 50
opts$extend_bp <- 0
opts$dif_acc_files <- c("E7.5Mesoderm_vs_E7.5EndodermEctoderm_H3K27ac_distal_E7.5_Mes_intersect12.txt.gz",
                         "E7.5Ectoderm_vs_E7.5MesodermEndoderm_H3K27ac_distal_E7.5_Ect_intersect12.txt.gz",
                         "E7.5Endoderm_vs_E7.5MesodermEctoderm_H3K27ac_distal_E7.5_End_intersect12.txt.gz",
                        "E7.5Epiblast_vs_E7.5Ectoderm_H3K27ac_distal_E7.5_Ect_intersect12.txt.gz")


# fun
fread_gz <- function(filename, ...){
  f <- file(filename)
  type <- summary(f)$class
  close.connection(f)
  if (type == "gzfile") {
    filename <- paste("zcat", filename)
    return(fread(cmd = filename, ...))
  }
  fread(filename, ...)
}

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")

find_anno <- function(id){
  sp <- strsplit(id, "_")
  map_chr(sp, ~.[c(1:4)] %>% paste(collapse = "_"))
}

################################################################################

# find lineage-associated sites

dif_acc <- paste0(io$dif_acc, "/", opts$dif_acc_files) %>%
  .[file.exists(.)] %>%
  set_names(basename(.) %>% gsub(".txt.gz", "", .)) %>%
  map(fread_gz) %>%
  map(~.[sig == TRUE, id])



# find pairs

pairs <- expand.grid(opts$motfis, opts$motfis) %>%
  setDT() %>%
  .[Var1 != Var2]

motifs <- dir(io$motif_bed, full = TRUE, pattern = ".bed$") %>%
  map(fread) %>%
  rbindlist() %>%
  .[, tf := toupper(tf)] %>%
  .[, mid := round(mid, -1) %>% as.integer]  %>%  # round motif mid point to avoid multi motifs at the same location
  .[, anno := strsplit(id, "_") %>% map_chr(~paste(.[c(1:5)], collapse = "_"))] 

motifs_dif_acc <- map(dif_acc, ~motifs[id %in% .]) %>%
  map2(names(.), ~.x[, anno := .y]) %>%
  rbindlist()

motifs <- rbind(motifs, motifs_dif_acc, use.names = TRUE)

id_chr <- motifs[, .(id, chr)] %>%
  setkey() %>%
  unique()

total_loci <- motifs[, .(N = length(unique(id))), anno ]

co_motifs <- map(1:nrow(pairs), ~{
  
  motif_pair <- unlist(pairs[.])
  dt <- motifs[tf %in% motif_pair, .(tf, id, mid, anno)] %>%
    unique(by = c("tf", "id", "mid", "anno")) %>%
    split(by = "tf")
  
  if (length(dt)!= 2) return(NULL)
  
  by_chance <- map(dt, ~.[, .(n=.N), anno] ) %>%
    purrr::reduce(merge, by = "anno") %>%
    merge(total_loci, by = "anno") %>%
    .[, .(freq = N * (n.x / N) * (n.y / N)), .(anno)]
  
  
  
  #     map2(names(.), ~setnames(.x, "mid", paste0("mid_", .y))) %>%
  paired <- purrr::reduce(dt, merge, by = c("id", "anno"), allow.cartesian = TRUE) %>%
    .[, dist := mid.x - mid.y] %>%
    .[abs(dist) <= opts$max_dist & abs(dist) >= opts$min_dist] 
  
  name <- paste(motif_pair[1], motif_pair[2], sep = "_")
  
  raw <- paired[, length(unique(id)), anno]
  
  enrichment <- merge(raw, by_chance, by = "anno") %>%
    .[, .(anno = anno, enrichment = V1 / freq)] %>%
    setnames("enrichment", name)
  
  raw <- setnames(raw, "V1", name)
  
  list(raw = raw, enrichment = enrichment)
}) %>%
  purrr::compact() %>%
  purrr::transpose()

# process into matrix

co_motif_mat <- map(co_motifs, ~{
  compact(.) %>%
    purrr::reduce(merge, by = "anno") %>%
    split(by = "anno") %>%
    map(melt, id.var = "anno") %>%
    map(tidyr::separate, sep = "_", variable, c("tf1", "tf2")) %>%
    map(dcast, tf1~tf2, value = value) %>%
    map(setDF) %>%
    map(tibble::column_to_rownames, "tf1") %>%
    map(as.matrix)
})

################################################################################

# plot heatmaps



pheatmap(co_motif_mat$raw$H3K27ac_distal_E7.5_union_intersect12, cluster_cols = FALSE, cluster_rows = FALSE, main = "Raw counts at K27ac union sites")
pheatmap(co_motif_mat$enrichment$H3K27ac_distal_E7.5_union_intersect12, cluster_cols = FALSE, cluster_rows = FALSE, main = "Enrichment at K27ac union sites")


walk2(co_motif_mat$raw, names(co_motif_mat$raw), ~{
  title = paste("Raw counts at", .y)
  filename = paste0(io$out_dir, "/raw/", .y, ".pdf")
  dir.create(dirname(filename), recursive = TRUE)
  pheatmap(.x, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE, 
           main = title,
           filename = filename)
})

walk2(co_motif_mat$enrichment, names(co_motif_mat$enrichment), ~{
  title = paste("Enrichment at", .y)
  filename = paste0(io$out_dir, "/enrich/", .y, ".pdf")
  dir.create(dirname(filename))
  pheatmap(.x, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE, 
           main = title,
           filename = filename)
})


