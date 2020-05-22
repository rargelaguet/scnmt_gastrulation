
####################################################################
## Plot dimensionality reduction of scNMT-seq mapped to the atlas ##
####################################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/gastrulation/rna/mapping2atlas/Mapping2gastrulationAtlas.R")

################
## Define I/O ##
################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
io$mapping.results <- "/Users/ricard/data/gastrulation/rna/mapping_10x/mapping.rds"
io$outdir <- "/Users/ricard/data/gastrulation/rna/mapping_10x"

####################
## Define options ##
####################

opts <- list()

# Plotting options
opts$dot_size = 0.6
opts$dot_alpha = 1

################################
## Load 10x atlas information ##
################################

# Load sample metadata
meta_atlas <- readRDS(paste0(io$path2atlas, "/sample_metadata_atlas.rds")) %>% as.data.table 

# Filter
meta_atlas <- meta_atlas %>%
  .[stage%in%c("E6.5","E6.75","E7.0","E7.25","E7.5")]

# Load pre-computed UMAP
umap <- readRDS(paste0(io$path2atlas,"/tsne_coordinates.rds"))

# Load mapping 
mapping <- readRDS(io$mapping.results)$mapping


###################################
## Plot dimensionality reduction ##
###################################

# Prepare data.frame to plot
plot_df = as.data.frame(umap)
plot_df$cell = meta_atlas$cell
plot_df$sample = meta_atlas$sample
plot_df$stage = meta_atlas$stage
plot_df$cluster = meta_atlas$cluster
plot_df$celltype = meta_atlas$celltype
plot_df = plot_df[sample(nrow(plot_df), nrow(plot_df), replace = FALSE), ]

################
## Plot atlas ##
################

p <- ggplot(data = plot_df, mapping = aes(x=V1, y=V2, colour=celltype)) +
  # geom_point(size=opts$dot_size, alpha=opts$dot_alpha) +
  ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="t-SNE Dimension 1", y="t-SNE Dimension 2") +
  scale_colour_manual(values = celltype_colours, name = "Stage") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_bw()

pdf(paste0(io$outdir,"/umap_rasterised.pdf"), width=13, height=7, useDingbats = F)
print(p)
dev.off()

##################
## Plot mapping ##
##################

index <- match(plot_df$cell, as.character(mapping$closest.cell) )
plot_df$mapped <- as.factor(!is.na(index))
  
p <- ggplot(data=plot_df, mapping=aes(x=V1, y=V2, colour=mapped)) +
  ggrastr::geom_point_rast(aes(size=mapped), alpha=opts$dot_alpha) +
  # geom_point(aes(size=mapped), alpha=opts$dot_alpha) +
  scale_size_manual(values = c("TRUE"=0.55, "FALSE"=0.15)) +
  labs(x="", y="") +
  scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_bw() +
  theme(
    legend.position = "none"
  )
# print(p)

pdf(paste0(io$outdir,"/umap_mapped_rasterised.pdf"), width=9, height=6)
print(p)
dev.off()