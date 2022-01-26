here::i_am("rna/differential/differential_expr.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$script <- here::here("rna/differential/differential_expr.R")
io$outdir <- file.path(io$basedir,"results/rna/differential"); dir.create(io$outdir, showWarnings = F)

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
  "E3.5_vs_E4.5" = list(c("E3.5_ICM"), c("E4.5_Epiblast","E4.5_Primitive_endoderm")),
  "E4.5Epiblast_vs_E4.5Primitive_endoderm" = list(c("E4.5_Epiblast"), c("E4.5_Primitive_endoderm")),
  "E5.5Epiblast_vs_E5.5Visceral_endoderm" = list(c("E5.5_Epiblast"), c("E5.5_Visceral_endoderm"))
)

for (i in names(opts$groups)) {
  groupA <- opts$groups[[i]][[1]]
  groupB <- opts$groups[[i]][[2]]
  outfile <- sprintf("%s/%s.txt", io$outdir, i)
  
  if (grepl("BI",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
    lsf <- sprintf("sbatch -n 1 --mem 8G --wrap")
  }
  
  cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --group_label stage_lineage3 --outfile %s'", lsf, io$script, paste(groupA, collapse=" "), paste(groupB, collapse=" "), outfile)
  print(cmd)
  system(cmd)
}
