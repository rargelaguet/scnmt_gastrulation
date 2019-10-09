library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(stringr)
library(doParallel)
source("/bi/group/reik/Stephen/gastrulation/fread_gz.R")

# this script generates accessibility profiles at motif sites. Profile data is 
# saved to disk for further processing and filtering

## I/O ##
io <- list()

io$basedir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$outdir <- "/bi/group/reik/Stephen/gastrulation/acc/motifs/profiles/out"
io$out_file <- paste0(io$outdir, "/profiles_at_motifs.csv")
io$data.indir <- paste(io$basedir,"acc/raw/",sep="/")


io$data.indir <- paste(io$basedir,"acc/raw",sep="/")
io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$features.indir <- paste0(io$basedir,"/features/sjc/motifs/motifdb")

## Options ##
opts <- list()
opts$window_size <- 1000 # in bp
opts$tile <- 30          # window size to calculate accessibility rates 
opts$numb_cells_per_group = 20 # take a random sample of this many cells for each group


opts$annotations = list.files(io$features.indir, pattern=".bed") %>% 
  str_replace(".bed", "")
names(opts$annotations) = opts$annotations
opts$annotations[1:length(opts$annotations)] = "center"

# Define stages and lineages
opts$stage_lineage <- c("E4.5", "E5.5","E6.5", "E7.5")

# Pseudobulk cells in the same stage_lineage?
opts$pseudobulk <- F

# Define which cells to use
opts$cells <- fread(io$in.sample_metadata) %>% 
  # .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[,stage_lineage:=stage] %>%
  .[pass_accQC==T & stage_lineage%in%opts$stage_lineage & KO_3b=="not"] %>%
  split(by="stage") %>%
  map(~.[sample(min(opts$numb_cells_per_group, .N)), id_acc]) %>%  
  # subset to xx cells per stage
  unlist(use.names=FALSE)


#load sample metadata
sample_metadata <- fread(io$in.sample_metadata) %>% .[id_acc%in%opts$cells] %>% 
  .[,c("id_acc","stage","lineage")] %>% setnames("id_acc","sample") %>%
  # .[,stage_lineage:=paste(stage,lineage,sep="_")]
  .[,stage_lineage:=stage]

#Load genomic contexts and define windows
anno_list <- list()
for (anno in names(opts$annotations)) {
  tmp <- fread(sprintf("%s/%s.bed",io$features.indir,anno))[,c(1,2,3,4,5,6)]
  colnames(tmp) <- c("chr","start","end","strand","id","anno")
  tmp <- tmp[complete.cases(tmp)]
  tmp[end<start, c("start", "end") := list(end, start)]
  
  # Define central position for the window approach
  if (opts$annotations[anno] == "start") {
    tmp <- rbind(tmp[strand=="+",.(chr,start,strand,id,anno)] %>% .[,center:=start] %>% .[,c("start"):=NULL], 
                 tmp[strand=="-",.(chr,end,strand,id,anno)] %>% .[,center:=end] %>% .[,c("end"):=NULL]) 
  }
  if (opts$annotations[anno] == "center") {
    stopifnot(all(tmp[,end] > tmp[,start]))
    tmp <- tmp[,.(chr,start,end,strand,id,anno)][,center:=round(end+start)/2][,c("start","end"):=NULL]
  }
  if (opts$annotations[anno] == "end") {
    tmp <- rbind(tmp[strand=="+",.(chr,end,strand,id,anno)][,center:=end][,c("end"):=NULL], 
                 tmp[strand=="-",.(chr,start,strand,id,anno)][,center:=start][,c("start"):=NULL])
  }
  anno_list[[anno]] <- tmp %>% .[, c("start","end") := list(center-opts$window_size,center+opts$window_size)]
}









files <- dir(io$data.indir, pattern=paste(opts$cells, collapse="|"), full=TRUE)
cells <- basename(files) %>% sub(".tsv.gz", "", .)
acc_raw <- map2(files, cells,  ~fread_gz(.x, select=c(1:2, 5)) %>% .[,sample := .y]) %>%
  rbindlist() %>%
  setnames("pos", "start") %>%
  .[, end := start] %>%
  setkey(chr, start, end)

# generate profiles and save to disk

sample_metadata <- unique(sample_metadata)
setkey(sample_metadata, sample)

dir.create(io$outdir, recursive=TRUE)


acc_dt <- map(anno_list, ~{
  anno_name <- .[,unique(anno)]
  .[, chr := as.character(chr)]
  setkey(acc_raw, chr, start, end)
  setkey(., chr, start, end)
  dt <- foverlaps(acc_raw, ., nomatch=0L) %>%
    .[,dist:=ifelse(strand %in% c("+","*"),i.start-center,center-i.start)] %>%
    .[, dist:=opts$tile*round(dist/opts$tile)] %>%
    .[, .(rate=mean(rate), .N), .(dist,anno,sample)] %>%
    setkey(sample) %>%
    merge(sample_metadata[, .(sample, stage_lineage)] %>% setkey(sample)) %>%
    .[, .(rate=mean(rate), sd=sd(rate)), .(dist, anno, stage_lineage)]
  plot <- ggplot(dt, aes(dist, rate, colour=stage_lineage, fill=stage_lineage)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=rate-sd, ymax=rate+sd), alpha=0.2, colour=NA) +
    ylim(0.1, 0.8) +
    ggtitle(anno_name)
  file_name <- paste0(io$outdir, "/", anno_name, ".pdf")
  save_plot(file_name, plot, base_height=6)
  print(plot)
  return(dt)        
}) 
acc_dt <- rbindlist(acc_dt)
fwrite(acc_dt, io$out_file)




