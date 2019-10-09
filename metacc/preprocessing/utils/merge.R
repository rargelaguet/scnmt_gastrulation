
##########################################
## Script to merge bismark output files ##
##########################################

# The output of bismark comes as two files per sample: R1 and R2
# This script merges them and generates a single file per sample.
# If the same site is present in the two files, then the reads are simply summed

library(data.table)
library(purrr)
library(doParallel)

indir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/E5.5_scBS/2683_2684"
outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/raw/E5.5_scBS/2683_2684/merged"; dir.create(outdir, showWarnings = F)
# cores <- 1

filenames <- list.files(indir,pattern=".*(.tsv.gz)$")
samples <- unique(substr(filenames,1,nchar(filenames)-10))
# samples <- unique(sapply(strsplit(filenames,split="_"), function(x) paste(x[c(1,2,3)],collapse="_") ) )

# registerDoParallel(cores=cores)
# invisible(foreach(i=1:length(samples)) %dopar% {
for (i in 1:length(samples)) {
  print(sprintf("%s (%d/%d)", samples[i], i, length(samples)))
  fname.in <- sprintf("%s/%s",indir,filenames[grep(paste0(samples[i],"_"),filenames)])
  
  if (length(fname.in) == 2) {
    dat1 <- fread(sprintf("zcat %s",fname.in[1]), sep=" ", header=F, select=c(1,2,3,4), stringsAsFactors=T, verbose=F, showProgress = F)
    dat2 <- fread(sprintf("zcat %s",fname.in[2]), sep=" ", header=F, select=c(1,2,3,4), stringsAsFactors=T, verbose=F, showProgress = F)
    colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
    dat <- rbind(dat1,dat2) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)
    
    # In some cases, samples have four files (WHY?)
  } else if (length(fname.in) == 4) {
    dat1 <- fread(sprintf("zcat %s",fname.in[1]), sep=" ", header=F, select=c(1,2,5,6), stringsAsFactors=T, verbose=F, showProgress = F)
    dat2 <- fread(sprintf("zcat %s",fname.in[2]), sep=" ", header=F, select=c(1,2,5,6), stringsAsFactors=T, verbose=F, showProgress = F)
    dat3 <- fread(sprintf("zcat %s",fname.in[3]), sep=" ", header=F, select=c(1,2,5,6), stringsAsFactors=T, verbose=F, showProgress = F)
    dat4 <- fread(sprintf("zcat %s",fname.in[4]), sep=" ", header=F, select=c(1,2,5,6), stringsAsFactors=T, verbose=F, showProgress = F)
    colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat3) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat4) <- c("chr","pos","met_reads","nonmet_reads")
    dat <- rbind(dat1,dat2,dat3,dat4) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)
    
  } else {
    stop("error")
  }
  fwrite(dat, file=sprintf("%s/%s.tsv",outdir,samples[i]), sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
}

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
# system(sprintf("pigz -p %d -f %s/*.tsv", cores, outdir))