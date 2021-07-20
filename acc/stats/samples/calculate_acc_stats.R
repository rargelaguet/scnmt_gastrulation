#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

io$outdir <- paste0(io$basedir,"/acc/results/stats")

# Define cells
opts$cells <- list.files(io$acc_data_raw, pattern = "*.tsv.gz") %>% 
# opts$cells <- list.files(io$acc_data_raw, pattern = "E3.5_*.tsv.gz") %>% 
  stringr::str_replace_all(".tsv.gz","") %>% head(n=3)

####################################################
## Load accessibility data and calculate statistics ##
####################################################

stats <- data.table(id_acc=opts$cells) %>% 
  .[,c("nreads","nGC","rate"):=as.numeric(NA)]

for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i))) {
    print(i)

    data <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), sep="\t", verbose=F, showProgress=F)

    # Compute genome-wide statistics
    stats[id_acc==i, c("nreads","nGC","rate"):=list(sum(data$acc_reads+data$nonacc_reads), nrow(data), round(100*mean(data$rate),2))]
    # stats[id_acc==i, c("nreads","nGC","rate"):=list(sum(data$acc_reads+data$nonacc_reads), nrow(data), round(100*mean(data$rate),2))]

  } else {
    print(sprintf("Sample %s not found",i))
  }
}

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/sample_stats.txt")
fwrite(stats, outfile, sep="\t", na = "NA")
