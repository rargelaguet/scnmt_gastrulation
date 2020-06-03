
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$dir.lib <- "/Users/ricard/bbreg/lib/"
  io$script <- "/Users/ricard/scnmt_gastrulation/metaccrna/dynamics_individual_examples/boxplots_argparse.R"
} else {
  stop("Computer not recognised")
}

R.utils::sourceDirectory(io$dir.lib, modifiedOnly = FALSE)

io$dir.bayes <- paste0(io$scmet,"/betabinomial/bayes/stan")
# io$out.dir <- paste0(io$scmet,"/betabinomial/bayes/pdf")


opts$anno <- c(
  # "H3K27ac_distal_E7.5_union_intersect12_500",
  # "H3K27ac_distal_E7.5_union_intersect12",
  "prom_2000_2000"
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_Ect_intersect12"
)

###############
## Load data ##
###############

# Load HVFs
hvf_res <- list()
for (i in 1:length(opts$anno)) {
  stan <- readRDS(sprintf("%s/allcells_%s_vb.rds", io$dir.bayes, opts$anno[i]))
  hvf_res[[opts$anno[i]]] <- detect_hvf(stan, delta_e = NULL, delta_g = 0.25, efdr = 0.1)
  # hvf_res[[opts$anno[i]]] <- detect_hvf(stan, delta_e = 0.9, delta_g = NULL, efdr = 0.1)
}
rm(stan)

# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>% 
  .[,c("ens_id","symbol")]

################
## Parse data ##
################

dt <- hvf_res$prom_2000_2000$summary %>% as.data.table %>%
  setnames("feature_name","ens_id") %>%
  merge(gene_metadata, by="ens_id") %>%
  .[prob>0.95 & epsilon>0.5] %>%
  setorder(-prob,-epsilon)

hits <- dt$ens_id

##########
## Plot ##
##########

args <- list()
args$met.anno <- args$acc.anno <- "prom_2000_2000"
args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
args$outdir <- "/Users/ricard/data/scnmt_gastrulation/metaccrna/plot_individual_examples"

for (i in 1:nrow(dt)) {
  # args$gene <- args$met.id <- args$acc.id <- i
  args$met.id <- args$acc.id <- dt[i,ens_id]
  args$gene <- dt[i,symbol]
  
  cmd <- sprintf("Rscript %s --gene %s --met.id %s --met.anno %s --acc.id %s --acc.anno %s --stage_lineage %s --outdir %s", 
    io$script, args$gene, args$met.id, args$met.anno, args$acc.id, args$acc.anno, paste(args$stage_lineage,collapse=" "), args$outdir)
  system(cmd)
}

