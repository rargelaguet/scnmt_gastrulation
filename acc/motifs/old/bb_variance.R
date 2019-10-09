
## I/O ##
io <- list()
io$outdir <- "/homes/ricard/gastrulation/met/variability/betabinomial_model/out"
io$script <- "/homes/ricard/gastrulation/met/variability/betabinomial_model/beta_binomial.R"
io$tmpdir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met/test/tmp"

# io$outdir <- "/Users/ricard/gastrulation/met/variability/betabinomial_model/out"
# io$script <- "/Users/ricard/gastrulation/met/variability/betabinomial_model/beta_binomial.R"
# io$tmpdir <- "/Users/ricard/tmp"

## Options ##
opts <- list()

opts$stage_lineage <- list(
  "E4.5" = c("E4.5_EPI","E4.5_PE"),
  "E5.5" = c("E5.5_EPI"),
  "E6.5" = c("E6.5_EPI","E6.5_PS","E6.75_EPI","E6.75_PS"),
  "E7.5" = c("E7.5_Ectoderm","E7.5_Endoderm","E7.5_Mesoderm")
)

opts$anno <- c(
  "prom_2000_2000_noncgi", "prom_2000_2000_cgi",
  # "exons", "introns",
  "LINE", "LTR",
  "E3.5_Distal_H3K27ac", "E6.5_Distal_H3K27ac","Wei_Distal_K27ac_intersect",
  # "E3.5_H3K27ac", "E6.5_H3K27ac", "Wei_K27ac_intersect",
  # "Wei_Ect_K27ac", "Wei_Mes_K27ac", "Wei_End_K27ac",
  "CGI"
)


for (i in names(opts$stage_lineage)) {
  for (anno in opts$anno) {
    outfile <- sprintf("%s/%s_%s.txt", io$outdir, i, anno)
    lsf <- sprintf("bsub -M 8192 -n 1 -q research -o %s/%s_%s.txt", io$tmpdir, i, anno)
    # lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --stage_lineage %s --outfile %s", 
                   lsf, io$script, anno, paste(opts$stage_lineage[[i]], collapse=" "), outfile)
    system(cmd)
  }
}
