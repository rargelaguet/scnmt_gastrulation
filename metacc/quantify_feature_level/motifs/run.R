
##############
## Settings ##
##############

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation/metacc/quantify_feature_level/motifs/quantify_feature_level_motifs.R"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation/metacc/quantify_feature_level/motifs/quantify_feature_level_motifs.R"
} else {
  stop("Computer not recognised")
}

# I/O
io$features.dir <- sprintf("%s/features/motifs",io$basedir)
io$outdir <- sprintf("%s/motifs",io$met_data_parsed);  dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp");  dir.create(io$tmpdir, showWarnings = F)

#############
## Options ##
#############

opts$annos <- c(
  "prom_2000_2000"
)

# opts$context <- c("CG","GC")
opts$context <- "CG"

#########
## Run ##
#########

for (i in opts$annos %>% head(n=1)) {
  opts$anno_motif <- list.files(io$features.dir, pattern=sprintf("%s.*.bed.gz",i)) %>% gsub(".bed.gz","",.) 
  for (j in opts$anno_motif %>% head(n=3)) {
    for (k in opts$context) {
      if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 4096 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir,i)
      } else {
        lsf <- ""
      }
      cmd <- sprintf("%s Rscript %s --anno %s --context %s", lsf, io$script, j, k)
      print(cmd)
      # system(cmd)
    }
  }
}

