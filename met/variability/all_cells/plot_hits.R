if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

dt <- fread("/Users/ricard/data/gastrulation/metrna/overdispersion/metrna_variability_basics.csv")

###############
## Plot hits ##
###############

io$script <- "/Users/ricard/scnmt_gastrulation/metrna/dynamics_individual_examples/boxplots_argparse.R"

args <- list()
args$met.anno <- "prom_2000_2000"
# args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
args$stage_lineage <- c("E4.5_Epiblast","E4.5_Primitive_endoderm", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Primitive_Streak", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Mesoderm")
args$outdir <- "/Users/ricard/data/scnmt_gastrulation/metrna/plot_individual_examples"

# to.plot <- dt[color=="red"]
# to.plot <- dt[color=="blue"]
# to.plot <- dt[gamma>0.3 & bio.disp<(-2)]
# to.plot <- dt[epsilon<(-1)]
to.plot <- dt[symbol%in%c("Dppa2","Dppa5a","Spp1","Zfp42","Id3","Mesp1","Krt8","Phlda2","Mixl1")]
for (i in 1:nrow(to.plot)) {
  # args$gene <- args$met.id  <- i
  args$met.id <- to.plot[i,ens_id]
  args$gene <- to.plot[i,symbol]

  cmd <- sprintf("Rscript %s --gene %s --met.id %s --met.anno %s --stage_lineage %s --outdir %s",
                 io$script, args$gene, args$met.id, args$met.anno, paste(args$stage_lineage,collapse=" "), args$outdir)
  system(cmd)
}