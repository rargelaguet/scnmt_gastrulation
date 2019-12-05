################################################################
## Plot dimensionality reduction of cells mapped to the atlas ##
################################################################

# This script requires the cell metadata from the atlas, which contains the precomputed UMAP coordinates

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/scnmt_eb/rna/mapping/plot/plot_settings.R")

#########
## I/O ##
#########

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x"
io$path2query <- "/Users/ricard/data/gastrulation"
io$mapping <- "/Users/ricard/data/gastrulation/rna/results/mapping_10x/mapping.rds"
io$outdir <- "/Users/ricard/data/gastrulation/rna/results/mapping_10x"

#############
## Options ##
#############

opts <- list()

# Dot size
opts$size.mapped <- 0.6
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 1.0
opts$alpha.nomapped <- 0.5

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt"))

# Extract precomputed dimensionality reduction coordinates
umap <- meta_atlas[,c("cell","umapX","umapY")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#####################
## Load query data ##
#####################

# Load query cell metadata
meta_query <- fread(paste0(io$path2query, "/sample_metadata.txt"))# %>%
  # .[!lineage10x_2 %in% c("PGC","Allantois","NOIDEA",NA)] # remove weird lineages 

# Load precomputed mapping
mapping <- readRDS(io$mapping)$mapping
mapping.dt <- data.table(
  # id_rna            = mapping$cell, 
  id_rna            = stringr::str_replace_all(mapping$cell,"map_",""), 
  celltype.mapped = mapping$celltype.mapped,
  stage.mapped    = mapping$stage.mapped,
  closest.cell    = as.character(mapping$closest.cell)
)

###################################
## Plot dimensionality reduction ##
###################################

# Prepare query data.frame to plot
plot_df_query = mapping.dt %>% merge(meta_query, by=c("id_rna"))

# Prepare atlas data.frame to plot
plot_df_atlas = umap %>% merge(meta_atlas, by=c("cell"))
# plot_df_atlas <- plot_df_atlas[celltype!="ExE ectoderm"]

## Day 2 ##
plot_df_atlas[,index:=match(plot_df_atlas$cell, plot_df_query[,closest.cell] )]
plot_df_atlas[,mapped:=as.factor(!is.na(index))]

p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
  scale_size_manual(values = c("TRUE"=opts$size.mapped, "FALSE"=opts$size.nomapped)) +
  scale_alpha_manual(values = c("TRUE"=opts$alpha.mapped, "FALSE"=opts$alpha.nomapped)) +
  # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


pdf(paste0(io$outdir,"/umap_mapped.pdf"), width=4, height=4)
print(p)
dev.off()
