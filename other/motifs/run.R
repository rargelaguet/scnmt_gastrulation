
#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation/motifs/motif_matcher.R"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation/motifs/motif_matcher.R"
} else {
  stop("Computer not recognised")
}
io$outdir <- sprintf("%s/features/motifs",io$basedir);  dir.create(io$outdir, showWarnings = F)
io$tmpdir <- paste0(io$outdir,"/tmp");  dir.create(io$tmpdir, showWarnings = F)

#############
## Options ##
#############

opts$anno <- c(
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_union_intersect12",
  # "prom_2000_2000"
  "multiome_peaks"
)

if (is.null(opts$annos)) { opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.) }

opts$motif_set <- c("jaspar2020")   # c("jaspar2020", "cisbp")
opts$width <- 7
opts$cutOff <- 5e-05
opts$extend <- 50

#########
## Run ##
#########

for (i in opts$anno) {
  for (j in opts$motif_set) {

    if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 4096 -n 1 -q research-rh74 -o %s/%s_%s.txt", io$tmpdir,i,j)
    } else {
      lsf <- ""
    }

    cmd <- sprintf("%s Rscript %s --anno %s --motif_set %s --extend %s --width %d --cutOff %.10f --outprefix %s", lsf, io$script, i, j, opts$extend, opts$width, opts$cutOff, sprintf("%s/%s_%s",io$outdir,i,j))
      
    print(cmd)
    system(cmd)
  }
}

