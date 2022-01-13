# options(warn=-1)
suppressMessages(library(stringr))
# suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('--context', type="character",  help='cg/CG or gc/GC')
p$add_argument('--anno',    type="character",    help='Genomic annotation')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

args$context <- toupper(args$context)
stopifnot(args$context %in% c("CG","GC"))


## START TEST ##
# args <- list()
# args$context <- "CG"
# args$anno <- "prom_2000_2000"
## END TEST ##

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

if (args$context == "GC") {
  io$in.folder <- paste0(io$acc_data_raw,"/pseudobulk")
  io$out.folder <- paste0(io$acc_data_parsed,"/pseudobulk")
} else if (args$context=="CG") {
  io$in.folder <- paste0(io$met_data_raw,"/pseudobulk")
  io$out.folder <- paste0(io$met_data_parsed,"/pseudobulk")
}

####################
## Define options ##
####################

# Define genomic contexts
if (is.null(args$anno)) {
  args$anno <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}

# Define samples
samples <- list.files(io$in.folder, pattern="*.tsv.gz") %>% gsub(".tsv.gz","",.)

###########
## Print ##
###########

cat("\nProcessing samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotation: %s\n", paste(args$anno, collapse=" ")))
cat(sprintf("- Context:?: %s\n",args$context))
cat("\n")

######################
## Load annotations ##
######################

anno_list <- list()
for (i in 1:length(args$anno)) {
  anno_list[[i]] <- fread(sprintf("%s/%s.bed.gz",io$features.dir,args$anno[i]), header=F, select=c(1,2,3,4,5)) %>% 
    setnames(c("chr","start","end","strand","id")) %>%
    setkey(chr,start,end)
}
names(anno_list) <- args$anno


###########################
## Load and process data ##
########################### 

for (i in 1:length(samples)) {
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",samples[i],args$anno) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for %s...\n",samples[i],paste(args$anno,collapse=", "))) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",samples[i]))  
    
    filename <- sprintf("%s/%s.tsv.gz",io$in.folder,samples[i])
    print(filename)
    if (file.exists(filename)) {
      dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,samples[i]), sep="\t", verbose=F, showProgress=F) %>%
        .[!is.na(pos)] %>% setnames(c("chr", "pos", "met_reads", "nonmet_reads", "rate"))
      if (nrow(dat_sample)>1) {
        dat_sample <- dat_sample %>%
          # .[,c("chr","pos","rate")] %>%
          # .[,rate:=rate*100] %>%
          .[,c("start","end") := list(pos,pos)] %>% .[,pos:=NULL] %>%
          .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
          setkey(chr,start,end)

        # stopifnot(all(dat_sample$rate %in% c(0,100)))
        
        # Overlap data with annotations
        for (anno in args$anno) {
          # fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
          fname.out <- sprintf("%s/tmp/%s_%s.gz",io$out.folder,samples[i],anno)
          if (file.exists(fname.out)) {
            cat(sprintf("Annotation for %s with %s already found, loading...\n",samples[i],anno))
          } else {
            cat(sprintf("Annotating %s with %s annotations...\n",samples[i],anno))
            
            # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
            ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
              .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
              .[,c("sample","anno") := list(samples[i],anno)] %>%
              # .[,.(rate=round(mean(rate)), met_reads=sum(met_reads), total_reads=sum(met_reads+nonmet_reads)), keyby=.(sample,id,anno)] %>%
              # .[,.(rate=sum(met_reads)/sum(met_reads+nonmet_reads), met_reads=sum(rate==100), total_reads=.N), by = c("sample","id","anno")] %>%
              .[,.(
                # binomial_rate = round(100*sum(met_reads)/sum(met_reads+nonmet_reads)), 
                rate = round(mean(rate)), 
                met_reads = sum(met_reads), 
                total_reads = sum(met_reads+nonmet_reads)
              ), by = c("sample","id","anno")] %>%
              .[,c("sample","id","anno","met_reads","total_reads","rate")]
            
            # Store and save results
            fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE)
          }
        }
        rm(dat_sample)
      }
    }
  }
}


# Concatenate everything and save it
for (i in args$anno) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    # files <- list.files(paste(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample' | gzip > %s", io$out.folder, i, outfile))
    cat("\n")
  }
}
