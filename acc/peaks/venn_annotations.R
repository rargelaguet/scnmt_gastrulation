library(data.table)
library(purrr)
library(VennDiagram)
library(grDevices)

fread_gz <- function(file, ...) {
  paste("zcat", file) %>%
    fread(cmd = ., ...)
}

find_lineage <- function(string){
  pos <- regexpr("End|Ect|Mes", string)
  substr(string, pos, pos + 2)
}

io <- list()
opts <- list()
io$data_dir = "/bi/scratch/Stephen_Clark/gastrulation_data"
io$histone_files <- paste0(io$data_dir, "/features/ChIP")
io$peak_files <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/sjc/peaks/metacc"

opts$histone_regex <- "specific"
opts$peaks_regex <- "E7.5_"



histone_regex <- paste(opts$histone_regex, collapse = "|")
peaks_regex <- paste(opts$peaks_regex, collapse = "|")

his_files <- dir(io$histone_files, ".bed", full = TRUE) %>%
  .[grep(histone_regex, .)] %>%
  setNames(find_lineage(.))

histones <- map2(his_files, names(his_files), ~{
  fread(.x, colClasses = list("factor" = 1)) %>%
    setnames(c("chr", "start", "end")) %>%
    .[, anno := .y] %>%
    .[, id := paste0(anno, "_", .I)]
}) 

peak_files <- dir(io$peak_files, pattern = ".bed", full = TRUE) %>%
  .[grep(peaks_regex, .)] %>%
  setNames(find_lineage(.))

peaks <- map(peak_files, fread)

overlaps <- map(c("Ect", "End", "Mes"), ~{
  list(peaks[[.]], histones[[.]]) %>%
    map(setkey, chr, start, end) %>%
    purrr::reduce(foverlaps, nomatch = 0L)
})





venn_numbers <- map2(overlaps, pairwise, ~{
  anno1 <- histones[[.y[1]]][, .(anno = unique(anno), .N)]
  anno2 <- histones[[.y[2]]][, .(anno = unique(anno), .N)]
  
  list(anno1 = anno1[, anno],
       anno2 = anno2[, anno],
       cross.area = .x[, min(length(unique(id)), length(unique(i.id)))],
       area1 = anno1[, N],
       area2 = anno2[, N]
  )
  
})


venns <- map(venn_numbers, ~{
  venn.plot <- draw.pairwise.venn(.$area1, .$area2, .$cross.area, 
                                  category = c(.$anno1, .$anno2), 
                                  fill = c("blue", "red"), 
                                  alpha = c(0.5, 0.5),
                                  cat.pos = c(0, 0),
                                  cat.dist = c(-0.05, -0.05))
  grid.draw(venn.plot)
  grid.newpage()
})
venns