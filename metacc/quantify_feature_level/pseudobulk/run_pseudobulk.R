
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$tmpdir <- "/Users/ricard/scnmt_gastrulation/metacc/quantify_feature_level/tmp"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  io$tmpdir <- "/homes/ricard/scnmt_gastrulation/metacc/quantify_feature_level/tmp"
} else {
  stop("Computer not recognised")
}

#############
## Options ##
#############

# CpG or GpC
opts$contexts <- c(
  "CG",
  "GC"
)

# Define genomic contexts
opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)# %>% head(n=1)


#########
## Run ##
#########

for (i in opts$contexts) {
  for (j in opts$annos) {

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      # lsf <- sprintf("bsub -M 12000 -n 1 -q research-rh74 -o %s/%s_%s.txt", io$tmpdir,i,j)
      lsf <- sprintf("bsub -M 24000 -n 1 -q research-rh74 -o %s/%s_%s.txt", io$tmpdir,i,j)
    }

    cmd <- sprintf("%s Rscript quantify_feature_level_pseudobulk.R --context %s --anno %s", lsf,i,j)
    system(cmd)
  }
}
