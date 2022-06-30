here::i_am("rna/differential/differential_expr.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
io$script <- here::here("rna/differential/differential_expr.R")
io$outdir <- file.path(io$basedir,"results/rna/differential"); dir.create(io$outdir, showWarnings = F)

opts$group_label <- "stage_lineage3"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & !is.na(celltype)] %>%
  .[,stage_lineage:=paste(stage,celltype,sep="_")] %>%
  .[,stage_lineage2:=paste(stage,celltype2,sep="_")] %>%
  .[,stage_lineage3:=paste(stage,celltype3,sep="_")]

table(sample_metadata$stage_lineage3)

#######################################
## Run selected pairwise comparisons ##
#######################################

# opts$groups <- list(
#   
#   # Stage transition: E4.5 to E5.5
#   # "E4.5Epiblast_vs_E5.5Epiblast" = list(c("E4.5_Epiblast"), c("E5.5_Epiblast")),
#   
#   # Stage transition: E5.5 to E6.5
#   # "E5.5Epiblast_vs_E6.5Epiblast" = list(c("E5.5_Epiblast"), c("E6.5_Epiblast")),
#   # "E5.5Visceral_endoderm_vs_E6.5Visceral_endoderm" = list(c("E5.5_Visceral_endoderm"), c("E6.5_Visceral_endoderm")),
#   
#   # Stage transition: E6.5 to E7.5
#   # "E6.5Epiblast_vs_E7.5Epiblast" = list(c("E6.5_Epiblast"), c("E7.5_Epiblast")),
#   # "E6.5Epiblast_vs_E7.5Ectoderm" = list(c("E6.5_Epiblast"), c("E7.5_Ectoderm")),
#   # "E6.5Primitive_Streak_vs_E7.5Endoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Endoderm")),
#   # "E6.5Primitive_Streak_vs_E7.5Mesoderm" = list(c("E6.5_Primitive_Streak"), c("E7.5_Mesoderm")),
#   # "E6.5Primitive_Streak_vs_E7.5Primitive_Streak" = list(c("E6.5_Primitive_Streak"), c("E7.5_Primitive_Streak")),
#   
#   # Mixed E6.5 and E7.5
#   "E6.5E7.5Primitive_Streak_vs_E6.5E7.5Mesoderm" = list(c("E6.5_Primitive_Streak","E7.5_Primitive_Streak"), c("E6.5_Mesoderm","E7.5_Mesoderm")),
#   "E6.5E7.5Epiblast_vs_E6.5E7.5Primitive_Streak" = list(c("E6.5_Epiblast","E7.5_Epiblast"), c("E6.5_Primitive_Streak","E7.5_Primitive_Streak")),
#   "E6.5E7.5Primitive_Streak_vs_E7.5Endoderm" = list(c("E6.5_Primitive_Streak","E7.5_Primitive_Streak"), c("E7.5_Endoderm")),
#   
#   # E4.5 Lineage comparisons
#   "E4.5Epiblast_vs_E4.5Primitive_endoderm" = list(c("E4.5_Epiblast"), c("E4.5_Primitive_endoderm")),
#   
#   # E5.5 Lineage comparisons
#   "E5.5Epiblast_vs_E5.5Visceral_endoderm" = list(c("E5.5_Epiblast"), c("E5.5_Visceral_endoderm")),
#   
#   # E6.5 sublineage comparisons
#   "E6.5Epiblast_vs_E6.5Visceral_endoderm" = list(c("E6.5_Epiblast"), c("E6.5_Visceral_endoderm")),
#   "E6.5Epiblast_vs_E6.5Primitive_Streak" = list(c("E6.5_Epiblast"), c("E6.5_Primitive_Streak")),
#   "E6.5Primitive_Streak_vs_E6.5Mesoderm" = list(c("E6.5_Primitive_Streak"), c("E6.5_Mesoderm")),
#   
#   # E7.5 lineage comparison
#   "E7.5Ectoderm_vs_E7.5MesodermEndoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Mesoderm","E7.5_Endoderm")),
#   "E7.5Ectoderm_vs_E7.5Endoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Endoderm")),
#   "E7.5Ectoderm_vs_E7.5Mesoderm" = list(c("E7.5_Ectoderm"), c("E7.5_Mesoderm")),
#   
#   "E7.5Mesoderm_vs_E7.5EndodermEctoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm","E7.5_Endoderm")),
#   "E7.5Mesoderm_vs_E7.5Endoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Endoderm")),
#   "E7.5Mesoderm_vs_E7.5Ectoderm" = list(c("E7.5_Mesoderm"), c("E7.5_Ectoderm")),
# 
#   "E7.5Endoderm_vs_E7.5MesodermEctoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm","E7.5_Mesoderm")),
#   "E7.5Endoderm_vs_E7.5Mesoderm" = list(c("E7.5_Endoderm"), c("E7.5_Mesoderm")),
#   "E7.5Endoderm_vs_E7.5Ectoderm" = list(c("E7.5_Endoderm"), c("E7.5_Ectoderm"))
# )
# 
# opts$groups <- list(
#   # "E3.5_vs_E4.5" = list(c("E3.5_ICM"), c("E4.5_Epiblast","E4.5_Primitive_endoderm")),
#   "E3.5ICM_vs_E4.5Epiblast" = list(c("E3.5_ICM"), c("E4.5_Epiblast")),
#   "E3.5ICM_vs_E4.5Primitive_endoderm" = list(c("E3.5_ICM"), c("E4.5_Primitive_endoderm")),
#   "E4.5E5.5Epiblast_vs_E6.5E7.5Epiblast" = list(c("E4.5_Epiblast","E5.5_Epiblast"), c("E6.5_Epiblast","E7.5_Epiblast"))
#   # "E4.5Epiblast_vs_E4.5Primitive_endoderm" = list(c("E4.5_Epiblast"), c("E4.5_Primitive_endoderm")),
#   # "E5.5Epiblast_vs_E5.5Visceral_endoderm" = list(c("E5.5_Epiblast"), c("E5.5_Visceral_endoderm"))
# )
# 
# for (i in names(opts$groups)) {
#   groupA <- opts$groups[[i]][[1]]
#   groupB <- opts$groups[[i]][[2]]
#   name_groupA <- strsplit(i,"_vs_")[[1]][1]
#   name_groupB <- strsplit(i,"_vs_")[[1]][2]
#   outfile <- sprintf("%s/%s.txt.gz", io$outdir, i)
#   
#   if (grepl("BI",Sys.info()['nodename'])) {
#     lsf <- ""
#   } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
#     lsf <- sprintf("sbatch -n 1 --mem 7G --wrap")
#   }
#   
#   cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --name_groupA %s --name_groupB %s --group_label %s --outfile %s'", lsf, io$script, paste(groupA, collapse=" "), paste(groupB, collapse=" "), name_groupA, name_groupB, opts$group_label, outfile)
#   print(cmd)
#   system(cmd)
# }


##################################
## Run all pairwise comparisons ##
##################################

# print(table(sample_metadata[[opts$group_label]]))

opts$groups <- names(which(table(sample_metadata[[opts$group_label]])>=25))# %>% head(n=4)

for (i in 1:length(opts$groups)) {
  groupA <- opts$groups[[i]]
  for (j in i:length(opts$groups)) {
    if (i!=j) {
      groupA <- name_groupA <- opts$groups[[i]]
      groupB <- name_groupB <- opts$groups[[j]]
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)
      if (!file.exists(outfile)) {
        
        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 6G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --name_groupA %s --name_groupB %s --group_label %s --outfile %s'", lsf, io$script, paste(groupA, collapse=" "), paste(groupB, collapse=" "), name_groupA, name_groupB, opts$group_label, outfile)

        # Run
        print(cmd)
        system(cmd)
      }
    }
  }
}
