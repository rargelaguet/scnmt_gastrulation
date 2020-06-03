
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$dir.lib <- "/Users/ricard/bbreg/lib/"
} else {
  stop("Computer not recognised")
}

R.utils::sourceDirectory(io$dir.lib, modifiedOnly = FALSE)

io$dir.bayes <- paste0(io$scmet,"/betabinomial/bayes/stan")
io$outfile <- paste0(io$scmet,"/betabinomial/bayes/txt/hvf.txt.gz")


opts$anno <- c(
  "H3K27ac_distal_E7.5_union_intersect12_500",
  "H3K27ac_distal_E7.5_union_intersect12",
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

dt <- opts$anno %>% map(function(i) {
  stan <- readRDS(sprintf("%s/allcells_%s_vb.rds", io$dir.bayes, i))
  hvf <- detect_hvf(stan, delta_e = NULL, delta_g = 0.25, efdr = 0.1)$summary %>% 
    as.data.table %>% .[,anno:=i]
  return(hvf)
}) %>% rbindlist %>% .[,c("feature_idx"):=NULL]

##########
## Save ##
##########

fwrite(dt, io$outfile)
