io$diff.dir <- paste0(io$scmet,"/differential/bayes/txt")

# Define genomic contexts
# opts$anno <- c(
#   # "H3K27ac_distal_E7.5_union_intersect12_500",
#   # "H3K27ac_distal_E7.5_union_intersect12",
#   # "H3K4me3_E7.5_union",
#   # "prom_2000_2000"
#   "H3K27ac_distal_E7.5_Mes_intersect12_500" = "Mesoderm enhancers",
#   # "H3K27ac_distal_E7.5_Mes_intersect12",
#   "H3K27ac_distal_E7.5_End_intersect12_500" = "Endoderm enhancers",
#   # "H3K27ac_distal_E7.5_End_intersect12",
#   "H3K27ac_distal_E7.5_Ect_intersect12_500" = "Ectoderm enhancers"
#   # "H3K27ac_distal_E7.5_Ect_intersect12"
# )

# Mean
diff.mu.dt <- names(opts$anno) %>% 
  map(~ fread(sprintf("%s/%s_mu.txt.gz", io$diff.dir, .))) %>%
  rbindlist %>% .[,mean_diff:=mean_B-mean_A] %>%
  .[,c("id","anno","mean_diff","mean_diff_tail_prob")] %>%
  setnames(c("id","anno","mean_diff","mean_prob"))

# Residual overdispersion
diff.resdisp.dt <- names(opts$anno) %>%
  map(~ fread(sprintf("%s/%s_epsilon.txt.gz", io$diff.dir, .))) %>%
  rbindlist %>% .[,c("id","anno","res_disp_change","res_disp_diff_tail_prob")] %>%
  setnames(c("id","anno","res_disp_diff","res_disp_prob"))

# Overdispersion
diff.disp.dt <- names(opts$anno) %>% 
  map(~ fread(sprintf("%s/%s_gamma.txt.gz", io$diff.dir, .))) %>%
  rbindlist %>% .[,disp_change:=disp_B-disp_A] %>%
  .[,c("id","anno","disp_change","disp_diff_tail_prob")] %>%
  setnames(c("id","anno","disp_diff","disp_prob"))

# Merge 
diff.dt <- merge(diff.mu.dt, diff.disp.dt, by=c("id","anno")) %>%
  merge(diff.resdisp.dt, by=c("id","anno")) %>%
  .[,anno:=stringr::str_replace_all(anno,opts$anno)]

# Calculate quadrants
diff.dt %>%
  .[mean_diff>0 & disp_diff>0,quadrant:="2"] %>%
  .[mean_diff>0 & disp_diff<0,quadrant:="4"] %>%
  .[mean_diff<0 & disp_diff<0,quadrant:="3"] %>%
  .[mean_diff<0 & disp_diff>0,quadrant:="1"]