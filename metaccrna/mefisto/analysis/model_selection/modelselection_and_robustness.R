library(MOFA2)


#####################
## Define settings ##
#####################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$model.dir <- "/Users/ricard/data/gastrulation_norsync_stuff/mofa/hdf5/robustness"
  io$outdir <- "/Users/ricard/data/gastrulation_norsync_stuff/mofa/pdf/robustness"
} else {
  io$model.dir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_norsync_stuff/mofa/hdf5/robustness"
  io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_norsync_stuff/mofa/pdf/robustness"
}

opts <- list()
opts$ntrials <- 5

#############################
## Load precomputed models ##
#############################

MOFAlist <- list()
for (i in 1:opts$ntrials) {
  outfile <- sprintf("%s/test_%d.hdf5",io$model.dir,i)
  MOFAlist[[i]] <- load_model(outfile)
}

###########################
## Subset active factors ##
###########################

# for (i in 1:opts$ntrials) {
#   r2 <- MOFAlist[[i]]@cache$variance_explained$r2_per_factor
#   factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
#   MOFAlist[[i]] <- subset_factors(MOFAlist[[i]], which(apply(factors,1,sum) >= 1))
# }

##################################
## Assess robustness of factors ##
##################################

pdf(paste0(io$outdir,"/robustness_factors.pdf"), height=7, width=8.5)
compare_factors(MOFAlist, show_rownames=F, show_colnames=F)
dev.off()


#####################
## Model selection ##
#####################

pdf(paste0(io$outdir,"/model_selection_elbo.pdf"), height=7, width=8.5)
compare_elbo(MOFAlist)# + 
#   coord_cartesian(ylim = c(28736299-1000, 28736299+1000)) 
dev.off()

