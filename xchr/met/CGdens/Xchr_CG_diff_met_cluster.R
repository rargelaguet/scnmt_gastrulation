# Load packages and source data
library(data.table)
library(purrr)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(argparse)

source("/bi/group/reik/Carine/gastrulation/met/enhancers/utils.R")


## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-s', '--stage',  type="character",  nargs='+',  help='stage (E3.5, E4.5,...)')
p$add_argument('-c1',  '--chr1',  type="character",  nargs='+',  help='chromosome (chr1, chr2, chrX,...')
p$add_argument('-c2',  '--chr2',  type="character",  nargs='+',  help='chromosome nr (1, 2, X,...')
p$add_argument('-b', '--bins',  type="integer",  nargs='+',  help='number of bins (10000, 100000,...)')
args <- p$parse_args(commandArgs(TRUE))


## Define I/O ##
io <- list()
io$basedir <- "/bi/scratch/for_Christel_from_Carine/gastrulation_data"
io$sample_metadata <- "/bi/group/reik/Carine/sample_metadata_scNMT_sex.txt"
io$met.data.indir <- paste(io$basedir,"met/raw",sep="/")
io$outdir <- "/bi/group/reik/Carine/Results/CGbin"

## Define options ##
opts <- list()

# Define stages and lineages (only embryonic tissues! Xchr dynamics are different in extraembryonic)
opts$stage_lineage <- c("E4.5_EPI","E5.5_EPI","E6.5_EPI","E6.5_PS","E6.75_EPI","E6.75_PS","E7.5_Ectoderm", "E7.5_Mesoderm", "E7.5_Endoderm")
opts$stage <- args$stage

# Select chromosome to analyse
opts$chr1 <- args$chr1
opts$chr2 <- args$chr2
# Define bins
opts$nr_bins <- args$bins

# Filtering criteria
#opts$min.weight <- 3
#opts$min.coverage <- 0.3
#opts$fraction.sites <- 0.75
opts$min.n <- 5

# Define cells to use
opts$cells <- fread(io$sample_metadata) %>% 
  .[KO_3b == "not"] %>%
  .[,stage_lineage:=paste(stage,lineage,sep="_")] %>%
  .[pass_metQC==T & pass_sex==T & stage%in%opts$stage & stage_lineage%in%opts$stage_lineage,id_met]


## Load sample metadata ##
sample_metadata <- fread(io$sample_metadata) %>% 
  .[id_met %in% opts$cells] %>% 
  .[,stage_lineage:=paste(stage,lineage,sep="_")] 


## Load genome sequence, define bins, calculate CG density per bin ##

# Load genome sequence mouse
Mmusculus <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
chr <-Mmusculus[[opts$chr1]]

# Define bins
length_chr <- length(chr)  # can analyze for different chromosomes, chr3 similar length to X
nr_bins <- opts$nr_bins
step <- length_chr / nr_bins

start <- seq(from=1, to=(length_chr-step), by=step)

# Calculate CG density per bin
CGdens <- list()

for (i in 1:length(start)){
  tmp <- chr[start[i]:(start[i]+step)]
  CGdens$CG_density[i] <- dinucleotideFrequency(tmp)["CG"] / length(tmp)
  CGdens$bin[i] <- i
}
CGdens <- as.data.table(CGdens)


## Load and bin methylation data ##
result <- list()
met_list_bin <- list()

for (cell in opts$cells) {
  
  foo <- fread(sprintf("zcat < %s/%s.tsv.gz",io$met.data.indir,cell), showProgress=F) %>%
    .[,c("chr","pos","rate")] %>% .[,sample:=cell] %>% .[chr %in% opts$chr2]  # can analyze different chromosomes
  
  for(i in 1:length(start)){
    bar <- foo[pos >= start[i] & pos <= (start[i]+step)]
    if(length(bar$rate) != 0){
      result$rate[i] <- mean(bar$rate)
    }
    else{
      result$rate[i] <- 0
    }
    result$bin[i] <- i
    result$sample[i] <- cell
  }
  
  met_list_bin[[cell]] <- result
}
met_bin <- rbindlist(met_list_bin) %>% merge(sample_metadata[,c("sample", "stage", "stage_lineage", "trueSex")], by="sample")


## Parse data ##

# sum met data by sex
met_rate_bin <- met_bin[,.(meanRate=mean(rate)), by=c("bin", "stage", "trueSex")] %>% .[stage %in% opts$stage]

met_bin_F <- met_rate_bin[trueSex=="female"] %>% .[,meanRate_F:=meanRate] %>% .[,meanRate:=NULL] %>% .[,trueSex:=NULL]
met_bin_M <- met_rate_bin[trueSex=="male"] %>% .[,meanRate_M:=meanRate] %>% .[,meanRate:=NULL] %>% .[,trueSex:=NULL]

met_bin_all <- merge(met_bin_F, met_bin_M, by=c("bin", "stage")) %>% 
  .[!meanRate_M==0] %>% .[!meanRate_F==0] %>% 
  .[,ratioFM:=(meanRate_F/meanRate_M)] %>% 
  .[,log_ratio:=log2(ratioFM)] %>% 
  .[,diff:=meanRate_F-meanRate_M]

# Merge met data and CG dens
data <- merge(met_bin_all, CGdens, by="bin") %>% setkey(., CG_density)


## Save results ##
write.csv(data, file=(paste0(io$outdir, "/data_CGdens_bins_", opts$stage, "_chr", opts$chr2, "_", opts$nr_bins,".csv")))
