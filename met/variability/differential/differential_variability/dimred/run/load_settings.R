
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  setwd("/Users/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  setwd("/homes/ricard/scnmt_gastrulation/met/variability/differential/differential_variability/dimred/run")
} else {
  stop("Computer not recognised")
}

io$scmet.diff <- paste0(io$scmet,"/differential/bayes/txt")

####################
## Define options ##
####################

# Data filtering options
# opts$min.CpGs <- 1        # minimum number of CpG sites per feature and cell
# opts$min.cells <- 50     # minimum coverage (number of cells with at least min.CpG measurements)

# Threshold on the tail probability
opts$tail_prob_threshold <- 0.90

# Minimum log odds ratio
opts$LOR_threshold <- log(2.5)

opts$stage_lineage <- c(
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[stage=="E7.5"] %>%
  .[stage_lineage%in%opts$stage_lineage] %>%
  .[pass_metQC==T,c("id_met","sample","stage","lineage10x_2")]

############################
## User-defined functions ##
############################

plot_dimred <- function(to.plot, color.by) {
  ggplot(to.plot, aes(x=V1, y=V2)) +
    geom_point(aes_string(fill=color.by), colour="black", shape=21, stroke = 0.3) +
    labs(x="", y="") +
    theme_classic() +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "none",
      legend.title = element_blank()
    )
}
