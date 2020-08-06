#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/settings.R")

# Define I/O
io$outdir <- paste0(io$basedir,"/met/results/stats")

# Update sample metadata
sample_metadata <- sample_metadata %>% .[!is.na(id_met)]

#####################
## Calculate stats ##
#####################

stats <- data.table(expand.grid(sample_metadata$id_met,opts$chr)) %>% 
  setnames(c("id_met","chr")) %>%
  .[,c("nreads","coverage","mean"):=as.numeric(NA)]

for (i in sample_metadata$id_met) {
  if (file.exists(sprintf("%s/%s.tsv.gz",io$met_data_raw,i))) {
    print(i)

    # Load sample methylation data
    data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i)) %>%
      .[,chr:=paste0("chr",chr)]

    # Compute methylation statistics per chromosome
    for (j in opts$chr) {
      data_j <- data[chr==j]
      stats[id_met==i & chr==j, c("nreads","coverage","mean"):=list(sum(data_j$met_reads+data_j$nonmet_reads), nrow(data_j),round(mean(data_j$rate)*100,2))]  
    }

  } else {
    print(sprintf("Sample %s not found for methylation",i))
  }
}
# stats <- stats[complete.cases(stats)]
fwrite(stats, paste0(io$outdir,"/stats_per_chromosome.txt.gz"), sep="\t")
