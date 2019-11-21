###################################################################################
## Script to compute (in parallel) differential methylation at the feature level ##
###################################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/gastrulation/met/results/differential/diffmet_supervised.R"
  io$outdir <- "/Users/ricard/data/gastrulation/met/results/differential/test"
} else {
  io$script <- "/homes/ricard/gastrulation/met/results/differential/diffmet_supervised.R"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/results/differential/test"
  io$tmpdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/results/differential/tmp"
}
dir.create(io$outdir, showWarnings=F)


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
  "E7.5Endoderm_vs_E7.5Ectoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm")),

  # Stage transition: E6.5 to E7.5
  "E6.5Epiblast_vs_E7.5Epiblast" = list(c("E6.5_Epiblast"), c("E7.5_Epiblast")),
  "E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast"), c("E7.5_Ectoderm")),
  "E6.5Primitive_Streak_vs_E7.5Endoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Endoderm")),
  "E6.5Primitive_Streak_vs_E7.5Mesoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Mesoderm"))
)

opts$groups <- list(
  # "E5.5E6.5E7.5Epiblast_vs_E7.5Ectoderm" = list(c("E5.5_Epiblast","E6.5_Epiblast","E7.5_Epiblast"), c("E7.5_Ectoderm")),
  # "E6.5E7.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast","E7.5_Epiblast"), c("E7.5_Ectoderm")),
  "E5.5E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast","E5.5_Epiblast"), c("E7.5_Ectoderm")),
  "E5.5Epiblast_vs_E7.5Ectoderm" = list(c("E5.5_Epiblast"), c("E7.5_Ectoderm")),
  "E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast"), c("E7.5_Ectoderm"))
)

# Minimum number of cells per group
opts$min.cells <- 10

# Genomic contexts
opts$anno <- c(
  # "CGI",
  # "H3K27ac_distal_E7.5_Ect_500",
  # "H3K27ac_distal_E7.5_End_500",
  # "H3K27ac_distal_E7.5_Mes_500",
  # "H3K27ac_distal_E7.5_union_500",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K27ac_distal_E7.5_union_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_union_intersect12",
  # "H3K27ac_promoter_E7.5_Mes_500",
  # "H3K27ac_promoter_E7.5_Ect_500",
  # "H3K27ac_promoter_E7.5_End_500",
  # "H3K27ac_promoter_E7.5_union_500",
  # "H3K27ac_E7.5_Ect_500",
  # "H3K27ac_E7.5_End_500",
  # "H3K27ac_E7.5_Mes_500",
  # "H3K27ac_E7.5_union_500",
  "H3K4me3_E7.5_Ect",
  "H3K4me3_E7.5_End",
  "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_union",
  "genebody",
  "prom_2000_2000"
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi"
  # "window2000_step1000"
)

# opts$anno <- c("H3K27ac_distal_E7.5_Ect_intersect12")

for (group in names(opts$groups)) {
  stage_lineage1 <- opts$groups[[group]][[1]]
  stage_lineage2 <- opts$groups[[group]][[2]]
  for (anno in opts$anno) {
    outfile <- sprintf("%s/%s_%s.txt", io$outdir, group, anno)
    lsf <- sprintf("bsub -M 2048 -n 1 -q standard -o %s/%s_%s.txt", io$tmpdir, group, anno)
    # lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --stage_lineage1 %s --stage_lineage2 %s --min.cells %d --outfile %s", 
                   lsf, io$script, anno, paste(stage_lineage1, collapse=" "), paste(stage_lineage2, collapse=" "), opts$min.cells, outfile)
    system(cmd)
  }
}
