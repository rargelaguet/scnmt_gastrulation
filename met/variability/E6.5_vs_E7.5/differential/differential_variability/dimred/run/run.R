
#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run/fit.R"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run/fit.R"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir, "/met/results/variability/differential/differential_variability/dimred")
io$tmpdir <- paste0(io$basedir, "/met/results/variability/differential/differential_variability/dimred/tmp")

#############
## Options ##
#############

# Define genomic contexts
opts$anno <- list(
  "H3K27ac_distal_E7.5_union_intersect12" = "H3K27ac_distal_E7.5_union_intersect12",
  "H3K27ac_distal_E7.5_union_intersect12_500" = "H3K27ac_distal_E7.5_union_intersect12_500",
  "prom_2000_2000" = c("prom_2000_2000")
)
# opts$anno <- list("prom_2000_2000" = c("prom_2000_2000"))

# Number of highly variable genes
opts$number.hvg <- seq(25,500,by=25)

#########
## Run ##
#########

# Differential variability hits
for (i in names(opts$anno)) {
  for (j in opts$number.hvg) {
    outprefix <- sprintf("%s/diffvar_%s_%d", io$outdir, i, j)

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 3000 -n 1 -q research-rh74 -o %s/diffvar_%s_%d.txt", io$tmpdir,i,j)
    }
    # lsf <- ""

    cmd <- sprintf("%s Rscript %s --anno %s --hvg %d --outprefix %s", 
                   lsf, io$script, paste(opts$anno[[i]], collapse=" "), j, outprefix)
    system(cmd)
  }
}
