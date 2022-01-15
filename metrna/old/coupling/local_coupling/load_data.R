suppressPackageStartupMessages(library(SingleCellExperiment))

# Define functions
mean_sd <- function(x) data.frame(y=mean(x), ymin=mean(x)-sd(x)/2, ymax=mean(x)+sd(x)/2) 

######################
## Define settings  ##
######################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/settings.R")
  io$tss <- "/Users/ricard/data/hg38_regulation/promoters/TSS_mRNA.bed"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/human_embryo_multiomics/settings.R")
  io$tss <- "/hps/nobackup2/research/stegle/users/ricard/hg38_regulation/promoters/TSS_mRNA.bed"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/metrna/coupling")

# Define window options
opts$up <- 3000
opts$down <- 3000
opts$window <- 250
opts$slide <- 50

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[pass_metQC==TRUE & pass_rnaQC==TRUE]

# [TESTING MODE] Subset cells to reduce memory burden  ##
# sample_metadata <- sample_metadata[,head(.SD,n=1),by="lineage"]

#########################
## Load gene metadata  ##
#########################

tss <- fread(io$tss, colClasses=c("factor","integer","integer","factor","factor","factor")) %>%
  setnames(c("chr","start","end","id","score","strand")) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,score:=NULL] %>%
  .[,tss:=start] %>%
  .[,c("start","end"):=.(start-opts$up,end+opts$down)] %>%
  setkey(chr,start,end)

# Split between CGI or non-CGI promoters
# non_cgi_promoters <- fread("/Users/ricard/data/gastrulation/features/filt/prom_2000_0_noncgi.bed",stringsAsFactors=TRUE)[,V5]
# cgi_promoters <- fread("/Users/ricard/data/gastrulation/features/filt/prom_2000_0_cgi.bed",stringsAsFactors=TRUE)[,V5]
# tss <- tss[,type:=ifelse(id%in%non_cgi_promoters,'non_cgi','cgi')]
# tss <- tss[type=="non_cgi"]

####################
## Load RNA data  ##
####################

# Load SingleCellExperiment object
sce <- readRDS(io$rna.sce)[,sample_metadata$id_rna]


##########
## Run  ##
##########

# Define genomic window around TSS
tmp <- seq(from=0-opts$up, to=0+opts$down-opts$window, by=opts$slide)
foo <- data.table(window_center=tmp+(opts$window/2), rel_start=tmp, rel_end=tmp+opts$window) %>%
  .[,(c("window_center","rel_start","rel_end")):=lapply(.SD, as.integer),.SDcols=(c("window_center","rel_start","rel_end"))] %>%
  setkey(rel_start,rel_end)

dt <- list()
for (i in sample_metadata$id_met) {
  print(i)
  
  # Define RNA expression data.table
  id.rna <- sample_metadata[id_met==i,id_rna]
  rna_dt <- logcounts(sce[,id.rna]) %>% t %>% as.data.table(keep.rownames="id_rna") %>% 
    melt(id.vars="id_rna", value.name="expr", variable.name="id")
  
  # Define DNA methylation data.table
  dt[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), showProgress = F, header = T,
        select = c("chr"="factor", "pos"="integer", "rate"="integer")) %>%
    .[,c("start","end"):=list(pos,pos)] %>%
    setnames("pos","bp") %>% 
    setkey("chr","start","end") %>%
    
    # Overlap with TSS annotations
    foverlaps(.,tss, nomatch=0) %>% 
    .[,dist:=ifelse(strand %in% c("+","*"),bp-tss,tss-bp)] %>% 
    .[,dist:=as.integer(opts$slide*round(dist/opts$slide))] %>%
    .[, c("chr","i.start","i.end","bp","tss","strand","start","end") := NULL] %>%

    # Overlap data with windows 
    .[,c("rel_start","rel_end"):=dist] %>% setkey(rel_start,rel_end) %>%
    foverlaps(foo) %>% 
    .[,c("i.rel_start","i.rel_end","rel_start","rel_end"):=NULL] %>%
    .[,.(rate=mean(rate)), by=c("id","window_center")] %>%
    
    # Merge with RNA expression
    merge(rna_dt, by="id") %>%
    
    # Calculate correlation coefficient per window
    .[,.(r = cor(x=rate, y=expr, method="pearson")), by = c("id_rna","window_center")]
}  
  
dt <- rbindlist(dt)

##########
## Save ##
##########

fwrite(dt, paste0(io$outdir,"/data.txt.gz"))
