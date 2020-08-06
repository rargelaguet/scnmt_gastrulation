#####################
## Define settings ##
#####################

io <- list(); opts <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$standard.mnn.script <- "/Users/ricard/scnmt_gastrulation/rna/iterative_mapping/standard_mnn.R"
  io$iterative.mnn.script <- "/Users/ricard/scnmt_gastrulation/rna/iterative_mapping/iterative_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  io$standard.mnn.script <- "/homes/ricard/scnmt_gastrulation/rna/iterative_mapping/standard_mnn.R"
  io$iterative.mnn.script <- "/homes/ricard/scnmt_gastrulation/rna/iterative_mapping/iterative_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation/rna/results/iterative_mapping/tmp"
} 

# Atlas stages
opts$atlas_stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75"
  # "E8.0"
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)

# Test mode (subset cells)?
opts$test <- FALSE


#########
## Run ##
#########

# LSF
if (grepl("ricard",Sys.info()['nodename'])) {
  lsf <- ""
} else {
  lsf <- sprintf("bsub -M 50000 -n 1 -q research-rh74 -o %s/mapping_mnn.txt", io$tmpdir)
}

# Run standard MNN
# cmd <- sprintf("%s Rscript %s --atlas_stages %s", lsf, io$standard.mnn.script, paste(opts$atlas_stages, collapse=" "))
# if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
# system(cmd)

# Run tree-guided MNN
cmd <- sprintf("%s Rscript %s --atlas_stages %s", lsf, io$iterative.mnn.script, paste(opts$atlas_stages, collapse=" "))
if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
system(cmd)
