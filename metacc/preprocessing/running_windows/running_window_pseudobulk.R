library(data.table)
library(purrr)

# Script to summarise methylation data using an unsupervised running window approach

#########
## I/O ##
#########

io <- list()
io$data_dir <- "/Users/ricard/data/gastrulation"
io$sample_meta <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$pseudobulk.dir <- paste0(io$data_dir, "/met/raw/pseudobulk")
io$out.dir <- "/Users/ricard/data/gastrulation/met/parsed/pseudobulk"; dir.create(io$out.dir, showWarnings=F)

#############
## Options ##
#############

opts <- list()
opts$stage_lineage <- c(
  "E4.5_EPI","E4.5_PE",
  "E5.5_EPI","E5.5_PE",
  "E6.5_EPI","E6.5_PS","E6.5_VE",
  "E7.5_Mesoderm","E7.5_Ectoderm","E7.5_Endoderm"
)
opts$chr <- c(1:19,"X")
opts$window.size <- 5000
opts$step.size <- 500

##########################
## Load pseudobulk data ##
##########################

data <- lapply(opts$stage_lineage, function(n) {
  fread(sprintf("zcat < %s/%s.tsv.gz",io$pseudobulk.dir,n), showProgress=F, stringsAsFactors=T, header=T, sep="\t")
}) %>% rbindlist
colnames(data) <- c("chr","pos","group","Nmet","N","rate")
data %>% setnames("group","sample")

############################
## Define running windows ##
############################

# Prepare windows
# NOTE: BE CAREFUL, ALL GROUPS MUST BE PROCESSED AT THE SAME TIME TO OBTAIN THE SAME EXACT WINDOWS
windows <- data[, .(start=seq(min(pos), max(pos), by=opts$step.size)), chr] %>%
  .[, end:=start+opts$window.size] %>%
  .[,id:=as.factor(paste0("window_",chr,":",start,"-",end))] %>%
  .[,anno:=as.factor(sprintf("window%s_step%s",opts$window.size,opts$step.size))] %>%
  setkey(chr,start,end)

#############
## Overlap ##
#############

# Prepare methylation data
data <- data %>% setnames("pos","start") %>%
  .[,end:=start] %>% setkey(chr,start,end) 

# Do the overlap (by chromosome to save memory)
ov_list <- list()
for (i in opts$chr) {
  print(i)
  ov_list[[i]] <- foverlaps(data[chr==i], windows[chr==i], nomatch=0, mult="all") %>%
    .[,c("i.start","i.end"):=NULL] %>%
    .[,.(rate=round(mean(rate)), N=.N), by=c("sample","id","anno","chr","start","end")]
}
ov <- rbindlist(ov_list)
rm(ov_list)

##########
## Save ##
##########

# select and re-order relevant columns
ov %>% setcolorder(c("sample","chr","start","end","id","anno","rate"))

# Write
fwrite(ov,
  file=sprintf("%s/window%d_step%d.tsv",io$out.dir, opts$window.size, opts$step.size),
  sep="\t", quote=F, col.names=F, row.names=F, na="NA"
)

# Compress
system(sprintf("gzip -f %s/*.tsv",io$out.dir))
