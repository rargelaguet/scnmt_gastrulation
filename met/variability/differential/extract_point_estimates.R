#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$dir.lib <- "/Users/ricard/bbreg/lib/"
  io$dir.bayes <- paste0(io$scmet,"/differential/bayes")
  io$out.dir <- paste0(io$scmet,"/differential/bayes/txt")
  source("/Users/ricard/scnmt_gastrulation/met/variability/differential/utils.R")
} else {
  stop("Computer not recognised")
}

io$pdf.dir <- paste0(io$out.dir,"/pdf")

# Source project lib directory
R.utils::sourceDirectory(io$dir.lib, modifiedOnly = FALSE)

if (!dir.exists(io$out.dir)) { dir.create(io$out.dir, showWarnings = FALSE) }
if (!dir.exists(io$pdf.dir)) { dir.create(io$pdf.dir, showWarnings = FALSE) }

#############
## Options ##
#############

# Define genomic contexts
opts$anno <- c(
  "H3K27ac_distal_E7.5_union_intersect12_500",
  "H3K27ac_distal_E7.5_union_intersect12",
  "prom_2000_2000",
  "H3K4me3_E7.5_union"
)

# Define groups
opts$groups <- c("E6.5", "E7.5")


#######################################################
## Load fitted objects and run differential analysis ##
#######################################################

dt_bayes <- diff_analysis <- list()
for (i in opts$anno) {
  
  # Load stan models
  dt_bayes[[i]] <- opts$groups %>% 
    map(~ readRDS(sprintf("%s/stan/%s_%s_vb.rds", io$dir.bayes,i,.)))
  
  # Run differential analysis  
  diff_analysis[[i]] <- differential_test(
    obj_A = dt_bayes[[i]][[1]], 
    obj_B = dt_bayes[[i]][[2]],
    group_label_A = opts$groups[1], 
    group_label_B = opts$groups[2]
  )
}


###########################
## Store point estimates ##
###########################

for (i in opts$anno) {
  
  # Save differential mean (mu)
  diff_mu <- diff_analysis[[i]]$mean_summary %>% as.data.table %>% 
    setnames("feature_name","id") %>%
    .[, evidence_thresh := diff_analysis[[i]]$diff_mean_thresh$evidence_thresh] %>%
    .[, efdr := diff_analysis[[i]]$diff_mean_thresh$efdr] %>%
    .[, efnr := diff_analysis[[i]]$diff_mean_thresh$efnr] %>%
    .[,anno := i]  %>%
    setorder(-mean_diff_tail_prob)
  fwrite(diff_mu, sprintf("%s/%s_mu.txt.gz",io$out.dir,i), sep="\t")
  
  # Save differential dispersion (gamma)
  diff_gamma <- diff_analysis[[i]]$disp_summary %>% as.data.table %>% 
    setnames("feature_name","id") %>%
    .[, evidence_thresh := diff_analysis[[i]]$diff_disp_thresh$evidence_thresh] %>%
    .[, efdr := diff_analysis[[i]]$diff_disp_thresh$efdr] %>%
    .[, efnr := diff_analysis[[i]]$diff_disp_thresh$efnr] %>%
    .[,anno := i] %>%
    setorder(-disp_diff_tail_prob)
  fwrite(diff_gamma, sprintf("%s/%s_gamma.txt.gz",io$out.dir,i), sep="\t")
  
  # Save differential residual dispersion (epsilon)
  diff_epsilon <- diff_analysis[[i]]$res_disp_summary %>% as.data.table %>% 
    setnames("feature_name","id") %>%
    .[, evidence_thresh := diff_analysis[[i]]$diff_res_disp_thresh$evidence_thresh] %>%
    .[, efdr := diff_analysis[[i]]$diff_res_disp_thresh$efdr] %>%
    .[, efnr := diff_analysis[[i]]$diff_res_disp_thresh$efnr] %>%
    .[,anno := i] %>%
    setorder(-res_disp_diff_tail_prob)
  fwrite(diff_epsilon, sprintf("%s/%s_epsilon.txt.gz",io$out.dir,i), sep="\t")
  
  # Save differential estimates using a simple gaussian model 
  # TO DO: CHECK THAT EXCITATORY IS GROUP A
  # diff_gaussian <- names(dt_bayes[[i]]) %>% 
  #   map(function(j) dt_bayes[[i]][[j]][["Y"]] %>% 
  #       .[,rate:=met_reads/total_reads] %>% 
  #       .[,.(mean=sum(met_reads)/sum(total_reads), var=var(rate)), by="Feature"] %>%
  #       .[,c("anno","group"):=list(i,j)]
  # ) %>% rbindlist %>% setnames("Feature","id") %>%
  #   .[group=="Excitatory",group:="A"] %>% .[group=="Inhibitory",group:="B"] %>%
  #   data.table::dcast(id+anno~group, value.var=c("mean","var"))
  #         
  # fwrite(diff_gaussian, sprintf("%s/%s_gaussian.txt.gz",io$out.dir,i))
}
