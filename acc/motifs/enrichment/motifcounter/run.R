
#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/hipsci_scMT/utils.R")
  io$script <- "/Users/ricard/hipsci_scMT/met/differential/feature_level/motif_enrichment/motif_enrichment.R"
  lsf <- ""
} else if (grepl("yoda",Sys.info()['nodename'])) {
  source("/homes/ricard/hipsci_scMT/utils.R")
  io$script <- "/homes/ricard/hipsci_scMT/met/differential/feature_level/motif_enrichment/motif_enrichment.R"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/met/differential/feature_level/motif_enrichment/v2")
io$tmpdir <- paste0(io$outdir,"/tmp")

dir.create(io$outdir, showWarnings = F)
dir.create(io$tmpdir, showWarnings = F)

#############
## Options ##
#############

opts$anno <- c(
  # "CGI",
  # "H3K27ac",
  # "H3K27ac_500_500",
  "distal_H3K27ac",
  # "distal_H3K27ac_500_500",
  # "H3K36me3",
  # "H3K4me1",
  # "H3K4me3",
  # "H3K4me3_500_500",
  # "H3K27me3",
  # "H3K27me3_500_500",
  # "prom_2000_2000",
  "prom_2000_2000_noncgi",
  "prom_2000_2000_cgi"
  # "genebody"
)


# Conditions to test enrichment
opts$conditions <- c("UD","D3")

# Minimum differential DNA methylation for significant hits
opts$min.diff <- 0.25

# order of the markov model
opts$order <- c(0,1)
# opts$order <- c(0)

# Testing mode (subsetting the number of motifs)
opts$test <- FALSE

for (i in opts$anno) {
  for (j in opts$order) {
    for (k in opts$conditions) {

      if (grepl("yoda",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M 2048 -n 1 -q standard -o %s/%s_%s_%s.txt", io$tmpdir,i,j,k)
      }

      cmd <- sprintf("%s Rscript %s --condition %s --anno %s --min.diff %.2f --order %d", lsf, io$script, k, i, opts$min.diff, j)
      if (opts$test) cmd <- paste0(cmd," --test")
        
      print(cmd)
      system(cmd)
    }

  }
}
