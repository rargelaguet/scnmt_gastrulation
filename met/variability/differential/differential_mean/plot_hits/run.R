
#########
## I/O ##
#########

source("/Users/ricard/scnmt_gastrulation/settings.R")
io$script <- "/Users/ricard/scnmt_gastrulation/met/plot_examples/boxplots.R"
io$diff.dir <- "/Users/ricard/data/betabinomial/argelaguet2019/met/variability/differential/bayes/txt"
io$outdir <- "/Users/ricard/data/betabinomial/argelaguet2019/met/variability/differential/differential_mean/pdf"

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
opts$mean_LOR_threshold <- log(2.5)

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s_mu.txt.gz",io$diff.dir,.))
) %>% rbindlist

dt[,abs_mean_diff:=abs(mean_A-mean_B)]

table(dt$mean_diff_test)

##########################################################
## Plot hits that are more methylated in E6.5 cells ##
##########################################################

# Genomic features 
diff <- dt[mean_diff_test=="E6.5+"] %>%
  .[mean_diff_tail_prob>=opts$tail_prob_threshold & abs(mean_LOR)>=opts$mean_LOR_threshold] %>% 
  setorder(-abs_mean_diff) %>%
  head(n=opts$nhits)

io$outdir2 <- paste0(io$outdir,"/E6.5+"); dir.create(io$outdir2, showWarnings = F)

# Boxplots of the DNA methylation rate per cell type
for (i in opts$anno) {
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script, i, paste(diff[anno==i,id],collapse=" "), io$outdir2)
  system(cmd)
}


##########################################################
## Plot hits that are mores methylated in E7.5 cells ##
##########################################################

# Genomic features 
diff <- dt[mean_diff_test=="E7.5+"] %>%
  .[mean_diff_tail_prob>=opts$tail_prob_threshold & abs(mean_LOR)>=opts$mean_LOR_threshold] %>% 
  setorder(-abs_mean_diff) %>%
  head(n=opts$nhits)
io$outdir2 <- paste0(io$outdir,"/E7.5+"); dir.create(io$outdir2, showWarnings = F)


# Boxplots of the DNA methylation rate per cell type
for (i in opts$anno) {
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script, i, paste(diff[anno==i,id],collapse=" "), io$outdir2)
  system(cmd)
}
