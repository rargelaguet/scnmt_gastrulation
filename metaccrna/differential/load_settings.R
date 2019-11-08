#####################
## Define settings ##
#####################

## Define I/O ##

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  stop()
}
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$met.stats <- paste0(io$basedir,"/met/results/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/samples/sample_stats.txt")
io$rna.file <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$annos_dir  <- paste0(io$basedir, "/features/filt")
io$outdir <- paste0(io$basedir,"/metaccrna/mesendoderm_commitment/ectoderm")

io <- list()

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/gastrulation/met/results/differential"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/results/differential"
io$diff.rna <- "/Users/ricard/data/gastrulation/rna/results/differential"

## Define options ##
opts <- list()

# Define genomic contexts for methylation
opts$met.annos <- c(
  # "genebody"="Gene body",
  "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  # "H3K4me3_E7.5_End"="End- H3K4me3",
  # "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  # "genebody"="Gene body",
  "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  # "H3K4me3_E7.5_End"="End- H3K4me3",
  # "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

# Overlap differential sites with nearby genes
opts$overlapGenes <- FALSE

# window length for the overlap between genes and features
opts$gene_window <- 25000

# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2

opts$min.fdr <- 0.10
opts$min.acc.diff <- 5
opts$min.met.diff <- 5

# Lineage colors
opts$colors <- c(
  Ectoderm = "steelblue", 
  Endoderm = "#43CD80", 
  Mesoderm = "violetred"
)

######################
## Define functions ##
######################

gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value)) +
    geom_bar(aes(fill=assay), color="black", stat="identity", position="dodge", size=0.25) +
    scale_fill_manual(values=c("met"="#F37A71", "acc"="#00BFC4")) +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title=i, x="", y="Number of hits") +
    theme_bw() +
    theme(
      plot.title = element_text(size=11, face='bold', hjust=0.5),
      axis.text = element_text(size=rel(1.0), color='black'),
      axis.text.x = element_text(size=rel(1.0), angle=60, hjust=1, vjust=1, color="black"),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position="none"
    )
  
  return(p)
}

theme_pub <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size=rel(1.2), angle=60, hjust=1, vjust=1, color="black"),
    axis.text.y = element_text(size=rel(1.2), color="black"),
    axis.title.y = element_text(size=rel(1.2), color="black"),
    legend.position = "right"
    )
}