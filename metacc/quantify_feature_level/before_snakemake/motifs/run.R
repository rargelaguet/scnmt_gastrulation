
##############
## Settings ##
##############

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation/metacc/quantify_feature_level/motifs/quantify_feature_level_motifs.R"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation/metacc/quantify_feature_level/motifs/quantify_feature_level_motifs.R"
} else {
  stop("Computer not recognised")
}

#########
## I/O ##
#########

#############
## Options ##
#############

opts$annos <- c(
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_union_intersect12",
  # "prom_2000_2000"
  "multiome_peaks"
)

# opts$context <- c("CG","GC")
opts$context <- "GC"

opts$motif.annotation <- "jaspar2020"




#########
## Run ##
#########

for (i in opts$annos) {
  opts$anno_motif <- list.files(io$motifs.dir, pattern=sprintf("%s_%s_.*.bed.gz",i,opts$motif.annotation)) %>% gsub(".bed.gz","",.)# %>% head(n=2)
  for (j in opts$anno_motif) {
    for (k in opts$context) {

      if (k=="CG") {
        opts$memory <- 4*1e3
        io$tmpdir <- paste0(io$met_data_motifs,"/tmp")
      } else if (k=="GC") {
        opts$memory <- 25*1e3
        io$tmpdir <- paste0(io$acc_data_motifs,"/tmp")
      } else {
        stop()
      }

      if (grepl("ebi",Sys.info()['nodename'])) {
        lsf <- sprintf("bsub -M %d -n 1 -q research-rh74 -o %s/%s_%s.txt", opts$memory,io$tmpdir,j,k)
        # lsf <- sprintf("bsub -M %d -n 1 -q research-rh74", opts$memory)
      } else {
        lsf <- ""
      }
      cmd <- sprintf("%s Rscript %s --anno %s --context %s", lsf, io$script, j, k)
      print(cmd)
      system(cmd)
    }
  }
  
  # for (j in opts$anno_motif %>% head(n=3)) {
  #   for (k in opts$context) {
  #     if (grepl("ebi",Sys.info()['nodename'])) {
  #       lsf <- sprintf("bsub -M 4096 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir,i)
  #     } else {
  #       lsf <- ""
  #     }
  #     cmd <- sprintf("%s Rscript %s --anno %s --context %s", lsf, io$script, j, k)
  #     print(cmd)
  #     system(cmd)
  #   }
  # }
  
}

