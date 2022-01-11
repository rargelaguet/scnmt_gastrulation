#################################################################
## Script to compute (in parallel) differential RNA expression ##
#################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/gastrulation/rna/differential/differential_expr.R"
  io$outdir <- "/Users/ricard/data/gastrulation/rna/differential"; dir.create(io$outdir)
} else {
  io$script <- "/homes/ricard/gastrulation/rna/differential/differential_expr.R"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/rna/differential"; dir.create(io$outdir, showWarnings=F)
  io$tmpdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/rna/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
}



## Options ##
opts <- list()

opts$groups <- list(
  
  # Stage transition: E4.5 to E5.5
  # "E4.5Epiblast_vs_E5.5Epiblast" = list(c("E4.5_Epiblast"), c("E5.5_Epiblast")),
  
  # Stage transition: E5.5 to E6.5
  # "E5.5Epiblast_vs_E6.5Epiblast" = list(c("E5.5_Epiblast"), c("E6.5_Epiblast")),
  # "E5.5Visceral_endoderm_vs_E6.5Visceral_endoderm" = list(c("E5.5_Visceral_endoderm"), c("E6.5_Visceral_endoderm")),
  
  # Stage transition: E6.5 to E7.5
  # "E6.5Epiblast_vs_E7.5Epiblast" = list(c("E6.5_Epiblast"), c("E7.5_Epiblast")),
  # "E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast"), c("E7.5_Ectoderm")),
  # "E6.5Primitive_Streak_vs_E7.5Endoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Endoderm")),
  # "E6.5Primitive_Streak_vs_E7.5Mesoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Mesoderm")),
  # "E6.5Primitive_Streak_vs_E7.5Primitive_Streak" = list(c("E6.5_Primitive_Streak"), c("E7.5_Primitive_Streak")),
  
  # Mixed E6.5 and E7.5
  "E6.5E7.5Primitive_Streak_vs_E6.5E7.5Mesoderm" = list(c("E6.5_Primitive_Streak","E7.5_Primitive_Streak"), c("E6.5_Mesoderm","E7.5_Mesoderm")),
  "E6.5E7.5Epiblast_vs_E6.5E7.5Primitive_Streak" = list(c("E6.5_Epiblast","E7.5_Epiblast"), c("E6.5_Primitive_Streak","E7.5_Primitive_Streak")),
  "E6.5E7.5Primitive_Streak_vs_E7.5Endoderm" = list(c("E6.5_Primitive_Streak","E7.5_Primitive_Streak"), c("E7.5_Endoderm")),
  
  # E4.5 Lineage comparisons
  "E4.5Epiblast_vs_E4.5Primitive_endoderm" = list(c("E4.5_Epiblast"), c("E4.5_Primitive_endoderm")),
  
  # E5.5 Lineage comparisons
  "E5.5Epiblast_vs_E5.5Visceral_endoderm" = list(c("E5.5_Epiblast"), c("E5.5_Visceral_endoderm")),
  
  # E6.5 sublineage comparisons
  "E6.5Epiblast_vs_E6.5Visceral_endoderm" = list(c("E6.5_Epiblast"), c("E6.5_Visceral_endoderm")),
  "E6.5Epiblast_vs_E6.5Primitive_Streak" = list(c("E6.5_Epiblast"), c("E6.5_Primitive_Streak")),
  "E6.5Primitive_Streak_vs_E6.5Mesoderm" = list(c("E6.5_Primitive_Streak"), c("E6.5_Mesoderm")),
  
  # E7.5 lineage comparison
  "E7.5Ectoderm_vs_E7.5MesodermEndoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Mesoderm","E7.5_Endoderm")),
  "E7.5Ectoderm_vs_E7.5Endoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Endoderm")),
  "E7.5Ectoderm_vs_E7.5Mesoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Mesoderm")),
  
  "E7.5Mesoderm_vs_E7.5EndodermEctoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm","E7.5_Endoderm")),
  "E7.5Mesoderm_vs_E7.5Endoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Endoderm")),
  "E7.5Mesoderm_vs_E7.5Ectoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm")),

  "E7.5Endoderm_vs_E7.5MesodermEctoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm","E7.5_Mesoderm")),
  "E7.5Endoderm_vs_E7.5Mesoderm" = list(c("E7.5_Endoderm"), c("E7.5_Mesoderm")),
  "E7.5Endoderm_vs_E7.5Ectoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm"))
)

opts$groups <- list(
  "E5.5E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast","E5.5_Epiblast"), c("E7.5_Ectoderm"))
)

for (group in names(opts$groups)) {
  stage_lineage1 <- opts$groups[[group]][[1]]
  stage_lineage2 <- opts$groups[[group]][[2]]
  outfile <- sprintf("%s/%s.txt", io$outdir, group)
  # lsf <- sprintf("bsub -M 8192 -n 1 -q standard -o %s/%s.txt", io$tmpdir, group)
  lsf <- ""
  cmd <- sprintf("%s Rscript %s --stage_lineage1 %s --stage_lineage2 %s --outfile %s", 
                 lsf, io$script, paste(stage_lineage1, collapse=" "), paste(stage_lineage2, collapse=" "), outfile)
  system(cmd)
}
