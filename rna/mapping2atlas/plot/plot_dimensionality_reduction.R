
####################################################################
## Plot dimensionality reduction of scNMT-seq mapped to the atlas ##
####################################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/gastrulation/rna/mapping2atlas/Mapping2gastrulationAtlas.R")

#####################
## I/O and options ##
#####################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
io$path2scNMT <- "/Users/ricard/data/gastrulation"
io$outdir <- "/Users/ricard/data/gastrulation/mapping_10x/out"

opts <- list()

# Plotting options
opts$dot_size = 0.6
opts$dot_alpha = 1

####################
## Load 10x atlas ##
####################

# sce_atlas  <- readRDS(paste0(io$path2atlas, "/SingleCellExperimentAtlas.rds"))
meta_atlas <- readRDS(paste0(io$path2atlas, "/sample_metadata_atlas.rds")) %>% as.data.table 

# Filter
meta_atlas <- meta_atlas %>%
  .[stage%in%c("E6.5","E6.75","E7.0","E7.25","E7.5")]

####################
## Load scNMT-seq ##
####################

sce_nmt  <- readRDS(paste0(io$path2scNMT, "/rna/parsed/SingleCellExperiment.rds"))
meta_nmt <- fread(file = paste0(io$path2scNMT, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)

# Filter
meta_nmt <- meta_nmt %>% .[pass_rnaQC==T & stage%in%c("E6.5","E7.5") ,]
sce_nmt <- sce_nmt[,meta_nmt$id_rna]

#########
## PCA ##
#########

# Load pre-computed PCA coordinates
pca <- readRDS(paste0(io$path2atlas, "/corrected_pcas.rds"))$all
pca <- pca[meta_atlas$cell,]

###########
## t-SNE ##
###########

# Run t-SNE
# set.seed(42)
# tsne = Rtsne::Rtsne(pca, check_duplicates = FALSE, pca = FALSE, perplexity = 120)
# saveRDS(tsne, paste0(io$path2atlas,"/tsne_coordinates.rds"))

# Load pre-computed t-SNE
# tsne <- readRDS(paste0(io$path2atlas,"/tsne_coordinates.rds"))

###########
## UMAP ##
###########

# library(umap)
# umap.defaults$n_neighbors <- 20
# umap.defaults$min_dist <- 0.7
# 
# set.seed(42)
# umap <- umap(pca, config = umap.defaults)$layout
# saveRDS(umap, paste0(io$path2atlas,"/umap_coordinates.rds"))

# Load pre-computed UMAP
umap <- readRDS(paste0(io$path2atlas,"/tsne_coordinates.rds"))

###################################
## Plot dimensionality reduction ##
###################################

# Prepare data.frame to plot
# plot_df = as.data.frame(tsne$Y)
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

# p <- ggplot(data = plot_df, mapping = aes(x = V1, y = V2, col = factor(stage))) +
#   geom_point(size = opts$dot_size, alpha = opts$dot_alpha) +
#   labs(x = "t-sne 1", y = "t-sne 2") +
#   scale_colour_manual(values = stage_colours, labels = stage_labels, name = "Stage") +
#   ggtitle("Stage")  +
#   guides(colour = guide_legend(override.aes = list(size=6)))

# pdf(paste0(path2scNMT,"tsne_clusters.pdf"))
# print(p)
# dev.off()

p <- ggplot(data = plot_df, mapping = aes(x=V1, y=V2, colour=celltype)) +
  # geom_point(size=opts$dot_size, alpha=opts$dot_alpha) +
  ggrastr::geom_point_rast(size=opts$dot_size, alpha=opts$dot_alpha) +
  labs(x="t-SNE Dimension 1", y="t-SNE Dimension 2") +
  scale_colour_manual(values = celltype_colours, name = "Stage") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_bw()# +
  # theme(
  #   legend.position = "none"
  # )

pdf(paste0(io$outdir,"/umap_rasterised.pdf"), width=13, height=7, useDingbats = F)
print(p)
dev.off()

##################
## Plot mapping ##
##################

mapping <- readRDS("/Users/ricard/data/gastrulation/mapping_10x/out/mapping.rds")$mapping

index <- match(plot_df$cell, as.character(mapping$closest.cell) )
plot_df$mapped <- as.factor(!is.na(index))
  
p <- ggplot(data=plot_df, mapping=aes(x=V1, y=V2, colour=mapped)) +
  # geom_point(size=opts$dot_size, alpha=opts$dot_alpha) +
  ggrastr::geom_point_rast(aes(size=mapped), alpha=opts$dot_alpha) +
  # geom_point(aes(size=mapped), alpha=opts$dot_alpha) +
  scale_size_manual(values = c("TRUE"=0.55, "FALSE"=0.15)) +
  labs(x="t-SNE Dimension 1", y="t-SNE Dimension 2") +
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