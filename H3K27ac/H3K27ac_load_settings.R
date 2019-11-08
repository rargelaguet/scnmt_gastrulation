#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation/"
  io$outdir <- "/Users/ricard/data/gastrulation/H3K27ac/pdf"
} else {
  io$basedir <- "/Users/stapelc/Documents/gastrulation_data"
  # io$gene.metadata <- "/Users/stapelc/Documents/GastrulaProject/data/ensembl/mouse/v93/BioMart/Mmusculus_genes_BioMart.93_GRCm38.p6.txt"
  io$outdir <- "/Users/stapelc/Documents/GastrulaProject/Results/revisions/H3K27ac"
}
io$sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$features.dir <- paste0(io$basedir,"/features/filt")
io$gene.metadata <- paste0(io$basedir,"/features/genes/Mmusculus_genes_BioMart.87.txt")
io$rna <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

# Folders with enhancer data
io$enh <- paste0(io$basedir, "/H3K27ac/20190508_ProbeReport_enhancers.txt")  
io$enh.esc <- paste0(io$basedir, "/H3K27ac/20190508_ProbeReport_Marked_in_esc.txt") 
io$enh.brain <- paste0(io$basedir, "/H3K27ac/20190508_ProbeReport_Marked_in_brain.txt") 
io$enh.gut <- paste0(io$basedir, "/H3K27ac/20190508_ProbeReport_Marked_in_gut.txt")
io$enh.heart <- paste0(io$basedir, "/H3K27ac/20190508_ProbeReport_Marked_in_heart.txt") 


## Define options ##
opts <- list()

# Define genomic contexts
opts$annos <-c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers"
)

# number of differentially marked enhancers selected
opts$nr.diff <- 250

## Options enhancer overlap
# window length for the overlap between genes and features (for coupling enhancers to genes?)
opts$gene_window <- 25000


## Options RNA-seq data
# Define stages and lineages
opts$stage <- c("E4.5","E5.5","E6.5","E7.5")
opts$lineage10x_1 <- c("Anterior_Primitive_Streak","Ectoderm","Embryonic_endoderm","Epiblast","Mature_mesoderm","Nascent_mesoderm","Notochord","Primitive_Streak")
#opts$lineage10x_2 <- c("Anterior_Primitive_Streak","Ectoderm","Embryonic_endoderm","Epiblast","ExE_ectoderm","ExE_endoderm","Mature_mesoderm","Nascent_mesoderm","NOIDEA","Notochord","PGC","Primitive_endoderm","Primitive_Streak")
opts$lineage10x_2 <- c("Ectoderm","Endoderm","Epiblast","Mesoderm","Primitive_Streak")

# Select cells for RNA-seq analysis
opts$cells <- fread(io$sample_metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_rnaQC==T & pass_metQC==T & pass_accQC==T & stage%in%opts$stage & lineage10x_2%in%opts$lineage10x_2,sample]

# Define colors for plotting
opts$colors <- c(
  "brain" = "steelblue", 
  "ecto" = "lightsteelblue", 
  "esc" = "grey70", 
  "gut"="#43CD80", 
  "heart"="#CD3278", 
  "other"="black"
)

######################
## Define functions ##
######################

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
