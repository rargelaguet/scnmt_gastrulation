source("/Users/ricard/scnmt_gastrulation/settings.R")
io$diff.dir <- paste0(io$scmet,"/differential/bayes/txt")
io$outdir <- paste0(io$scmet,"/differential/metacc_profiles"); dir.create(io$outdir, showWarnings = F)

# Define which stage and lineages to plot
opts$stage_lineage <- c(
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)



# Define genomic contexts
opts$annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12_500"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12_500"="Endoderm enhancers"
  # "H3K27ac_distal_E7.5_Ect_intersect12_500"="Ectoderm enhancers"
  # "H3K27ac_distal_E7.5_union_intersect12"="All enhancers"
)

opts$window_size <- 2000
opts$met.tile <- 100
opts$acc.tile <- 50

# Differential methylation parmeters 
# opts$tail_prob_threshold <- 0.90  # Threshold on the tail probability
# opts$min.change <- 0.10           # Minimum rate change


############################
# Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[stage_lineage%in%opts$stage_lineage] %>%
  .[lineage10x!="Visceral_endoderm"] %>%
  .[stage=="E7.5" & lineage10x_2%in%c("Epiblast","Ectoderm"),lineage10x_2:="Epiblast/Ectoderm"] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]

opts$met.cells <- sample_metadata[pass_metQC==T,id_met]
opts$acc.cells <- sample_metadata[pass_accQC==T,id_acc]

sample_metadata <- sample_metadata[,c("sample","pass_metQC","pass_accQC","id_met","id_acc","stage","lineage10x_2","stage_lineage")]