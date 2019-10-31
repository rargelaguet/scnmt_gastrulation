
##########################################################################
## Script to do volcano plots from the differential expression analysis ##
##########################################################################

library(data.table)
library(purrr)
library(ggplot2)


#####################
## Define settings ##
#####################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
	io$input.dir <- "/Users/ricard/data/gastrulation/rna/differential"
	io$outdir <- "/Users/ricard/data/gastrulation/rna/differential/pdf"
	source("/Users/ricard/gastrulation/rna/differential/utils.R")
}

## Options ##
opts <- list()
opts$comparisons <- c(
  "E4.5Epiblast_vs_E7.5Ectoderm",
  "E5.5Epiblast_vs_E7.5Ectoderm",
  "E6.5Epiblast_vs_E7.5Ectoderm",
  "E7.5Epiblast_vs_E7.5Ectoderm"
)


###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  fread(cmd=sprintf("zcat < %s/%s.txt.gz",io$input.dir,i)) %>% .[,comparison:=i] 
) %>% rbindlist

###################
## Volcano plots ##
###################

for (i in unique(diff.results$comparison)) {
  
    p <- gg_volcano_plot(diff.results[comparison==i], top_genes = 0)
    
    pdf(sprintf("%s/%s.pdf",io$outdir,i), width=6, height=4)
    print(p)
    dev.off()
}

