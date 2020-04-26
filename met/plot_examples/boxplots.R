suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--id',    type="character",   nargs='+',   help='feature id(s)')
p$add_argument('--anno',  type="character",                help='genomic context')
p$add_argument('--outdir',  type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# args$id <- c(
#   "H3K27ac_distal_E7.5_union_intersect12_500_12616",
#   "H3K27ac_distal_E7.5_union_intersect12_500_11546",
#   "H3K27ac_distal_E7.5_union_intersect12_500_15743",
#   "H3K27ac_distal_E7.5_union_intersect12_500_16689",
#   "H3K27ac_distal_E7.5_union_intersect12_500_11470"
# )
# args$anno <- "H3K27ac_distal_E7.5_union_intersect12_500"

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  source("~/Research/Projects/epigenetic_heterogeneity/code/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

if (is.null(args$outdir)) args$outdir <- paste0(io$basedir,"/boxplots")

# Define filtering criteria
opts$min.cpg <- 3

opts$stage_lineage <- c(
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)
sample_metadata <- sample_metadata[stage_lineage%in%opts$stage_lineage]

###############################
## Load DNA methylation data ##
###############################

data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,args$anno), showProgress=F) %>%
  setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>%
  .[id%in%args$id]

# Filter by coverage
data <- data[Ntotal>=opts$min.cpg]

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata, by="id_met")

# Rate from 0 to 1
if (max(data$rate)>1) {
  data[,rate:=rate/100]
}

###############
## Box plots ##
###############

for (i in args$id) {

  p1 <- ggplot(data[id==i], aes(x=stage, y=rate, fill=stage)) +
    geom_jitter(size=1.5, alpha=0.75, width=0.25, shape=21) +
    geom_violin(alpha=0.5) +
    geom_boxplot(alpha=0.6, outlier.shape=NA, width=0.3) +
    # scale_fill_manual(values=opts$colors) +
    labs(x="", y="Methylation rate", title="") +
    theme_classic() +
    theme(
      axis.title.y = element_text(colour="black", size=rel(1.1)),
      # axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
      axis.text = element_text(colour="black", size=rel(1.0)),
      legend.position = "none"
    )

  p2 <- ggplot(data[id==i], aes(x=stage_lineage, y=rate, fill=lineage10x_2)) +
    geom_jitter(size=1.5, alpha=0.75, width=0.25, shape=21) +
    geom_violin(alpha=0.5, width=1.0) +
    geom_boxplot(alpha=0.6, outlier.shape=NA, width=0.5) +
    # facet_wrap(~`stage`, scales="free_x") +
    facet_grid(~stage, scales="free_x", space = "free_x") +
    scale_fill_manual(values=opts$colors_lineages) +
    labs(x="", y="Methylation rate", title="") +
    theme_classic() +
    theme(
      axis.title.y = element_text(colour="black", size=rel(1.1)),
      axis.text.x = element_text(colour="black", size=rel(1.0), angle=20, hjust=1),
      axis.text.y = element_text(colour="black", size=rel(1.0)),
      strip.text = element_text(colour="black", size=rel(1.2)),
      legend.position = "none"
    )

  p <- cowplot::plot_grid(plotlist=list(p1,p2), ncol=2, rel_widths = c(1/4,3/4))
  pdf(sprintf("%s/boxplot_%s.pdf",args$outdir,i), width=10, height=4, useDingbats = F)
  print(p)
  dev.off()
}
