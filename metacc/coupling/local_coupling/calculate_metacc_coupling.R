here::here("metacc/coupling/calculate_metacc_coupling.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--tss',  type="character",              help='TSS annotation file')
p$add_argument('--up',  type="integer", default=3000,              help='Window upstream')
p$add_argument('--down',  type="integer", default=3000,              help='Window downstream')
p$add_argument('--window',  type="integer", default=3000,              help='Window size')
p$add_argument('--tile',  type="integer", default=100,              help='Tile size')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outfile',  type="character",              help='Output file')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$tss <- "/Users/argelagr/data/mm10_regulation/genes/TSS_protein_coding.bed"
# args$up <- 3000
# args$down <- 3000
# args$window <- 150
# args$tile <- 50
# args$test <- FALSE
# args$outfile  <- file.path(io$basedir,"results/metacc/coupling/precomputed_metacc_tss_coupling.txt.gz")
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_metQC==TRUE & pass_accQC==TRUE]

# [TESTING MODE] Subset cells
if (args$test) {
  sample_metadata <- sample_metadata %>% .[,head(.SD,n=3),by="sample"]
}

print(table(sample_metadata$stage))

#########################
## Load gene metadata  ##
#########################

tss <- fread(args$tss, colClasses=c("factor","integer","integer","factor","factor","factor")) %>%
  setnames(c("chr","start","end","id","score","strand")) %>%
  # .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,score:=NULL] %>%
  .[,tss:=start] %>%
  .[,c("start","end"):=.(start-args$up,end+args$down)] %>%
  setkey(chr,start,end)

# Define genomic window around TSS
tmp <- seq(from=0-args$up, to=0+args$down-args$window, by=args$tile)
foo <- data.table(window_center=tmp+(args$window/2), rel_start=tmp, rel_end=tmp+args$window) %>%
  .[,(c("window_center","rel_start","rel_end")):=lapply(.SD, as.integer),.SDcols=(c("window_center","rel_start","rel_end"))] %>%
  setkey(rel_start,rel_end)

################################
## Calculate Met/Acc coupling ##
################################

metacc_coupling.dt <- sample_metadata$cell %>% map(function(i) {

  print(i)
  
  # Chromatin accessibility
  acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,sample_metadata[cell==i,id_acc]), showProgress = F, header = T,
                  select = c("chr"="character", "pos"="integer", "rate"="integer")) %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    setnames("pos","bp") %>% .[,c("start","end"):=list(bp,bp)] %>%
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(args$tile*round(dist/args$tile))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%
    
    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,.(acc_rate=mean(rate)), by=c("id","window_center")] %>%
    .[,cell:=i]
  
  if (nrow(acc_dt)==0) {
    warning(sprintf("No overlaps found between %s CpG methylation and TSS annotation...",i))
  }
  
  # DNA methylation
  met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,sample_metadata[cell==i,id_met]), showProgress = F, header = T,
        select = c("chr"="character", "pos"="integer", "rate"="integer")) %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    setnames("pos","bp") %>% .[,c("start","end"):=list(bp,bp)] %>%
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(args$tile*round(dist/args$tile))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%

    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,.(met_rate=mean(rate)), by=c("id","window_center")] %>%
    .[,cell:=i]
    
  if (nrow(met_dt)==0) {
    warning(sprintf("No overlaps found between %s CpG methylation and TSS annotation...",i))
  }
    
    # Calculate Met vs Acc correlation coefficient per window
   merge(met_dt, acc_dt, by = c("window_center","id","cell")) %>%
    .[,.(r = cor(x=met_rate, y=acc_rate, method="pearson") %>% round(2)), by = c("cell","window_center")] %>%
     return
}) %>% rbindlist

# print(length(unique(metacc_coupling.dt$cell)))

##########
## Save ##
##########

fwrite(metacc_coupling.dt, args$outfile)
