##############################################
## Script to compute RNA expression markers ##
##############################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/gastrulation/rna/find_markers/find_markers.R"
  io$outdir <- "/Users/ricard/data/gastrulation/rna/find_markers"; dir.create(io$outdir, showWarnings = F)
  # io$tmpdir <- "/Users/ricard/data/gastrulation_norsync_stuff/rna/test/tmp"; dir.create(io$tmpdir)  
} else {
  stop()
  # io$script <- "/homes/ricard/gastrulation/rna/find_markers/find_markers.R"
  # io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_norsync_stuff/rna/test"; dir.create(io$outdir)
  # io$tmpdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_norsync_stuff/rna/test/tmp"; dir.create(io$tmpdir)
}


## Options ##
opts <- list()

opts$groups <- list(
  
  # E4.5 markers
  # E5.5 markers

  # E6.5 markers
  # "E6.5Epiblast" = list(c("E6.5_Epiblast"), c("E6.5_Visceral_endoderm","E6.5_Primitive_Streak","E6.5_Mesoderm")),
  # "E6.5Visceral_endoderm" = list(c("E6.5_Visceral_endoderm"), c("E6.5_Epiblast","E6.5_Primitive_Streak","E6.5_Mesoderm")),
  "E6.5Mesoderm" = list(c("E6.5_Mesoderm"), c("E6.5_Epiblast","E6.5_Primitive_Streak","E6.5_Visceral_endoderm"))
  # "E6.5Primitive_Streak" = list(c("E6.5_Primitive_Streak"), c("E6.5_Epiblast","E6.5_Mesoderm","E6.5_Visceral_endoderm")),

  # E7.5 markers
  # "E7.5_Epiblast" = list(c("E7.5_Epiblast"), c("E7.5_Primitive_Streak","E7.5_Mesoderm","E7.5_Endoderm","E7.5_Ectoderm")),
  # "E7.5_Ectoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Primitive_Streak","E7.5_Mesoderm","E7.5_Endoderm","E7.5_Epiblast")),
  # "E7.5_Primitive_Streak" = list(c("E7.5_Primitive_Streak"), c("E7.5_Ectoderm","E7.5_Mesoderm","E7.5_Endoderm","E7.5_Epiblast")),
  # "E7.5_Mesoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm","E7.5_Primitive_Streak","E7.5_Endoderm","E7.5_Epiblast")),
  # "E7.5_Endoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm","E7.5_Primitive_Streak","E7.5_Mesoderm","E7.5_Epiblast"))
)


for (group in names(opts$groups)) {
  stage_lineage1 <- opts$groups[[group]][[1]]
  stage_lineage2 <- opts$groups[[group]][[2]]
  outfile <- sprintf("%s/%s.txt", io$outdir, group)
  # lsf <- sprintf("bsub -M 8192 -n 1 -q research -o %s/%s.txt", io$tmpdir, group)
  lsf <- ""
  cmd <- sprintf("%s Rscript %s --stage_lineage1 %s --stage_lineage2 %s --outfile %s", 
                 lsf, io$script, paste(stage_lineage1, collapse=" "), paste(stage_lineage2, collapse=" "), outfile)
  system(cmd)
}
