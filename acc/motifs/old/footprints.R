library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(doParallel)

# some TF motif sites appear to show footprinting. This script loads all motif
# annotations, selects those which (accross sites) show evidence of a footprint
# then uses these as an annotation to quantify accessibility at the motif site
# and flanking site. It does this for every cell at every position.

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$base_dir <- "/bi/group/reik/Stephen/github/gastrulation"

io$acc_motif_data <- paste0(io$base_dir, "/acc/motifs/profiles/out/accessibility_at_motifs.csv")
io$anno_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb/filt/motif_only")
io$plots_out <- paste0(io$base_dir, "/acc/motifs/out/footprints/plots")

io$sample_meta <- paste0(io$data_dir, "/sample_metadata.txt")
io$raw_data_dir <- paste0(io$data_dir, "/acc/raw")
io$data_out_dir <- paste0(io$data_dir, "/acc/parsed/motifs/footprints")

opts <- list()
opts$flank <- 60
opts$min_footprint <- 0.05

opts$anno_flank <- 60
opts$stages <- c("E4.5", "E5.5", "E6.5", "E7.5")
opts$KO_3b <- "not"

dir.create(io$plots_out, recursive=TRUE)

acc_data <- fread(io$acc_motif_data) 
center <- acc_data[dist==0, .(center=mean(rate)), .(anno)]
flank <- acc_data[dist!=0&(dist<=opts$flank&dist>=-opts$flank), .(flank=mean(rate)), .(anno)]

footprint <- merge(center, flank, by="anno") %>%
  .[flank-center>opts$min_footprint, anno]

plots <- map(footprint, ~{
  ggplot(acc_data[anno==.], aes(dist, rate, colour=stage_lineage, fill=stage_lineage)) + 
    geom_line() +
    geom_ribbon(aes(ymin=rate-sd, ymax=rate+sd), alpha=0.2, colour=NA) +
    ggtitle(.)
})


walk2(footprint, plots, ~{
  file_name <- paste0(io$plots_out, "/", .x, ".pdf")
  save_plot(file_name, .y)
})


# now quantify accessibility footprints at individual loci in each cell

motif <- dir(io$anno_dir, pattern=paste(footprint, collapse="|"), full=TRUE) %>%
  map(fread, colClasses=list(character=1), select=1:6) %>%
  rbindlist()

flank <- copy(motif)[, mid := start+(end-start)/2] %>%
  .[, c("start", "end", "mid") := .(mid-opts$anno_flank, mid+opts$anno_flank, NULL)]



meta <- fread(io$sample_meta) %>%
  .[pass_accQC==TRUE & stage %in% opts$stages & KO_3b %in% opts$KO_3b]

cells <- meta[, id_acc]

files <- paste0(io$raw_data_dir, "/", cells, ".tsv.gz") %>%
  .[file.exists(.)] %>%
  paste("zcat", .)

dir.create(paste0(io$data_out_dir,"/tmp"), recursive=TRUE)

anno <- list(motif=motif, flank=flank) %>%
  map(setkey, chr, start, end)

registerDoParallel(cores=opts$cores)
foreach(file=files) %dopar% {
  cell <- basename(file) %>%
    gsub(".tsv.gz", "", .)
  tmp_file <- paste0(io$data_out_dir, "/tmp/", cell, ".tsv")
  dt <- fread(file, select=c(1:2, 5)) %>%
    setnames("pos", "start") %>%
    .[, end := start] %>%
    setkey(chr, start, end)
  dt <- map(anno, ~foverlaps(dt, ., nomatch=0L) %>% 
              .[, .(rate=round(mean(rate*100)), .N), .(id, anno)])
  dt <- map2(dt, names(dt), ~setnames(.x, c("rate", "N"), c(paste0("rate_", .y), paste0("N_", .y))))
  walk(dt, setkey, id, anno)
  dt <- reduce(dt, merge)
  dt <- dt[, sample := cell]
  fwrite(dt, tmp_file, sep="\t")
}

# merge files and save

files <- dir(paste0(io$data_out_dir,"/tmp"), full=TRUE)

out_file <- paste0(io$data_out_dir, "/footprints.tsv")

# this only works without colnames....
# file.create(out_file)
# walk(files, ~file.append(out_file, .))

data <- map(files, fread) %>%
  rbindlist()


out_file <- paste0(io$data_out_dir, "/footprints.tsv")
fwrite(data, out_file)
system(paste0("gzip -f ", out_file))
file.remove(files)


