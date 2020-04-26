
#########
## I/O ##
#########

source("/Users/ricard/scnmt_gastrulation/settings.R")
io$script <- "/Users/ricard/scnmt_gastrulation/met/plot_examples/boxplots.R"
io$diff.dir <- "/Users/ricard/data/betabinomial/argelaguet2019/met/variability/differential/bayes/txt"
io$outdir <- "/Users/ricard/data/betabinomial/argelaguet2019/met/variability/differential/differential_variability/pdf"

#############
## Options ##
#############

opts$anno <- c(
  # "H3K27ac_distal_E7.5_union_intersect12_500"
  # "H3K27ac_distal_E7.5_union_intersect12",
  "prom_2000_2000"
)


# Top number of hits to plot
opts$nhits <- 25

# Threshold on the tail probability
opts$tail_prob_threshold <- 0.90

# Minimum log odds ratio
opts$LOR_threshold <- log(2.5)

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s_epsilon.txt.gz",io$diff.dir,.))
) %>% rbindlist

dt[,abs_epsilon_diff:=abs(res_disp_A-res_disp_B)]

# Call differential features
# diff.dt %>%
#   .[,sig:=abs(res_disp_LOR)>=opts$min.LOR & res_disp_diff_tail_prob>=opts$tail_prob] %>%
#   .[sig==F,res_disp_diff_test:="Nodiff"] %>%
#   .[sig==T & res_disp_LOR>0,res_disp_diff_test:=paste0(opts$groups[1],"+")] %>%
#   .[sig==T & res_disp_LOR<0,res_disp_diff_test:=paste0(opts$groups[2],"+")]

table(dt$res_disp_diff_test)

####################################################
## Plot hits that are more variable in E6.5 cells ##
####################################################

# Genomic features 
diff <- dt[res_disp_diff_test=="E6.5+"] %>%
  .[res_disp_diff_tail_prob>=opts$tail_prob_threshold & abs(res_disp_LOR)>=opts$LOR_threshold] %>% 
  setorder(-abs_epsilon_diff) %>%
  head(n=opts$nhits)
io$outdir2 <- paste0(io$outdir,"/E6.5+"); dir.create(io$outdir2, showWarnings = F)

# Boxplots of the DNA methylation rate per cell type
for (i in opts$anno) {
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script, i, paste(diff[anno==i,id],collapse=" "), io$outdir2)
  system(cmd)
}


####################################################
## Plot hits that are more variable in E7.5 cells ##
####################################################

# Genomic features 
diff <- dt[res_disp_diff_test=="E7.5+"] %>%
  .[res_disp_diff_tail_prob>=opts$tail_prob_threshold & abs(res_disp_LOR)>=opts$LOR_threshold] %>% 
  setorder(-abs_epsilon_diff) %>%
  head(n=opts$nhits)
io$outdir2 <- paste0(io$outdir,"/E7.5+"); dir.create(io$outdir2, showWarnings = F)

# Boxplots of the DNA methylation rate per cell type
for (i in opts$anno) {
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script, i, paste(diff[anno==i,id],collapse=" "), io$outdir2)
  system(cmd)
}
