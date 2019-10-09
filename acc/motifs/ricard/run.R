###################################################################################
## Script to compute (in parallel) differential accessibility at the feature level ##
###################################################################################

library(data.table)
library(purrr)

## I/O ##
io <- list()
io$data <- "/hps/nobackup/stegle/users/ricard/gastrulation/acc/parsed/fimo_motifs"
io$script <- "/homes/ricard/gastrulation/acc/motifs/ricard/diffacc_supervised.R"
io$outdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_norsync_stuff/acc/test/diff_motifs"; dir.create(io$outdir, showWarnings=F)
io$tmpdir <- "/hps/nobackup/stegle/users/ricard/gastrulation_norsync_stuff/acc/test/diff_featurelevel/tmp"; dir.create(io$tmpdir)

## Options ##
opts <- list()

opts$groups <- list(
  "E7.5Ect_vs_E7.5MesEnd" = list(c("E7.5_Ectoderm"), c("E7.5_Mesoderm","E7.5_Endoderm")),
  "E7.5Mes_vs_E7.5EctEnd" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm","E7.5_Endoderm")),
  "E7.5End_vs_E7.5EctMes" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm","E7.5_Mesoderm"))
)

# Minimum number of cells per group
opts$min.cells <- c(
  "E7.5Ect_vs_E7.5MesEnd" = 5,
  "E7.5Mes_vs_E7.5EctEnd" = 5,
  "E7.5End_vs_E7.5EctMes" = 5
)

# Genomic contexts
opts$anno <- list.files(io$data, pattern = "\\.gz$") %>%
  stringr::str_replace_all(".tsv.gz","")# %>% head(n=1)


for (group in names(opts$groups)) {
  print(group)
  stage_lineage1 <- opts$groups[[group]][[1]]
  stage_lineage2 <- opts$groups[[group]][[2]]
  for (anno in opts$anno) {
    print(anno)
    outfile <- sprintf("%s/%s_%s.txt", io$outdir, group, anno)
    lsf <- sprintf("bsub -M 4096 -n 1 -q research -o %s/%s_%s.txt", io$tmpdir, group, anno)
    # lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --stage_lineage1 %s --stage_lineage2 %s --min.cells %d --outfile %s", 
                   lsf, io$script, anno, paste(stage_lineage1, collapse=" "), paste(stage_lineage2, collapse=" "), opts$min.cells[group], outfile)
    system(cmd)
  }
}
