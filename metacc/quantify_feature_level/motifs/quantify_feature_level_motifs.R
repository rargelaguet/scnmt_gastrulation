suppressMessages(library(stringr))
suppressMessages(library(argparse))

#####################
## Parse arguments ##
#####################

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('--context',  type="character",              required=TRUE,  help='cg/CG or gc/GC')
p$add_argument('--annos',    type="character",  nargs="+",  required=TRUE,  help='Genomic annotation')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$context <- "GC"
# args$annos <- c("prom_2000_2000_jaspar2020_ARGFX_359")
## END TEST ##

# Define what context to look at: CG (MET) or GC (ACC)
args$context <- toupper(args$context)
stopifnot(args$context %in% c("CG","GC"))

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

#########
## I/O ##
#########

if (args$context == "GC") {
  io$in.folder <- io$acc_data_raw
  io$out.folder <- paste0(io$acc_data_parsed,"/motifs");  dir.create(io$out.folder)
} else if (args$context=="CG") {
  io$in.folder <- io$met_data_raw
  io$out.folder <- paste0(io$met_data_parsed,"/motifs");  dir.create(io$out.folder)
}

#############
## Options ##
#############

# Annotations to analyse
if (is.null(args$annos)) {
  args$annos <- list.files(io$motifs.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}

# Define cells
if (args$context=="CG") {
  samples_keep <- fread(io$metadata) %>% .[!is.na(id_met),id_met]
} else if (args$context=="GC") {
  samples_keep <- fread(io$metadata) %>% .[!is.na(id_acc),id_acc]
}

# ###########
# ## Print ##
# ###########

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$motifs.dir))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(args$annos, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",args$context))

######################
## Load annotations ##
######################

anno_list <- list()
for (i in 1:length(args$annos)) {
  anno_list[[i]] <- fread(sprintf("%s/%s.bed.gz",io$motifs.dir,args$annos[i]), header=T, select=c(1,2,3,5)) %>%
    setnames(c("chr","start","end","id")) %>%
    .[,chr:=as.factor(chr)] %>%
    setkey(chr,start,end)
}
names(anno_list) <- args$annos


#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  if (all(file.exists(sprintf("%s/tmp/%s/%s.gz",io$out.folder,args$annos,samples_keep[i])))) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",samples_keep[i])) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",samples_keep[i]))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$in.folder,samples_keep[i]), sep="\t", verbose=F, showProgress=F) %>%
      .[,c("chr","pos","rate")] %>%
      .[,c("start","end") := list(pos,pos)] %>% # Add 'start' and 'end' columns to do the overlap
      .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
      setkey(chr,start,end)
    
    # Sanity check
    stopifnot(all(dat_sample$rate %in% c(0,1)))
    
    # Overlap data with annotations
    for (anno in args$annos) {
      dir.create(sprintf("%s/tmp/%s",io$out.folder,anno))
      fname.out <- sprintf("%s/tmp/%s/%s.gz",io$out.folder,anno,samples_keep[i])
      if (file.exists(fname.out)) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",samples_keep[i],anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",samples_keep[i],anno))
        
        # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
          .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
          .[,c("sample","anno") := list(samples_keep[i],anno)] %>%
          # Compute number of methylated CpGs and the corresponding methylation rates
          .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), N=.N), keyby=c("sample","id","anno")] %>%
          # Reorder columns
          .[,c("sample","id","anno","Nmet","N","rate")]
        
        # Store and save results
        print(fname.out)
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE)
      }
    }
  }
}

warnings()

#####################################
## Concatenate everything and save ##
#####################################

for (i in args$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    system(sprintf("cat %s/tmp/%s/*.gz | zgrep -E -v 'sample|id_met|id_acc' | pigz > %s", io$out.folder, i, outfile))
  }
}

