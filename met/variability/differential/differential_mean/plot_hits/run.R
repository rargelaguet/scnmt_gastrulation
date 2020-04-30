
#########
## I/O ##
#########

source("/Users/ricard/scnmt_gastrulation/settings.R")
io$script <- "/Users/ricard/scnmt_gastrulation/met/plot_examples/boxplots.R"
io$diff.dir <- paste0(io$scmet,"/differential/bayes/txt")
io$outdir <- paste0(io$scmet,"/differential/differential_mean/pdf")

#############
## Options ##
#############

opts$anno <- c(
  "H3K27ac_distal_E7.5_union_intersect12_500",
  "H3K27ac_distal_E7.5_union_intersect12",
  "H3K4me3_E7.5_union",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12"
)

# Top number of hits to plot
opts$nhits <- 25

# Threshold on the tail probability
opts$tail_prob_threshold <- 0.90

#############################################
## Load results from differential analysis ##
#############################################

diff.dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s_mu.txt.gz",io$diff.dir,.))
) %>% rbindlist %>% .[,abs_mean_diff:=abs(mean_A-mean_B)]

######################################################
## Plot hits that are more methylated in E6.5 cells ##
######################################################

for (i in opts$anno) {

  # Select hits  
  hits <- diff.dt %>%
    .[anno==i & mean_diff_tail_prob>=opts$tail_prob_threshold & mean_LOR>0] %>% 
    setorder(-abs_mean_diff) %>%
    head(n=opts$nhits) %>% .$id
  
  # Plot
  io$outdir2 <- sprintf("%s/E6.5+/%s",io$outdir,i); dir.create(io$outdir2, showWarnings = F, recursive = T)
  cmd <- sprintf("Rscript %s --anno %s --id %s --min_cpg 3 --outdir %s", io$script, i, paste(hits,collapse=" "), io$outdir2)
  system(cmd)
}


######################################################
## Plot hits that are more methylated in E7.5 cells ##
######################################################

for (i in opts$anno) {

  # Select hits  
  hits <- diff.dt %>%
    .[anno==i & mean_diff_tail_prob>=opts$tail_prob_threshold & mean_LOR<0] %>% 
    setorder(-abs_mean_diff) %>%
    head(n=opts$nhits) %>% .$id
  
  # Plot
  io$outdir2 <- sprintf("%s/E7.5+/%s",io$outdir,i); dir.create(io$outdir2, showWarnings = F, recursive = T)
  cmd <- sprintf("Rscript %s --anno %s --id %s --min_cpg 3 --outdir %s", io$script, i, paste(hits,collapse=" "), io$outdir2)
  system(cmd)
}
