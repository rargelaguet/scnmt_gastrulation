# here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$query_metadata <- file.path(io$basedir,"results/rna/mapping/trajectories/blood/sample_metadata_after_mapping.txt.gz")
io$atlas_metadata <- file.path(io$atlas.basedir,"results/trajectories/blood_precomputed/blood_trajectory.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/mapping/trajectories/blood/pdf")

dir.create(io$outdir, showWarnings = F)
dir.create(file.path(io$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(io$outdir,"per_class"), showWarnings = F)


# Options

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(io$query_metadata) %>%
  .[!is.na(closest.cell)]

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas_metadata)# %>%
  # setnames(c("DC1","DC2"),c("V1","V2"))

##############################
## Define plotting function ##
##############################

plot.dimred <- function(plot_df, query.label, atlas.label = "Atlas") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(query.label, atlas.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(query.label, atlas.label)
  
  # Define dot colours  
  colour.values <- c("red", "lightgrey")
  names(colour.values) <- c(query.label, atlas.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}


###############################
## Plot one sample at a time ##
###############################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- meta_atlas %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

###############################
## Plot one class at a time ##
###############################

classes.to.plot <- unique(sample_metadata$class)

for (i in classes.to.plot) {
  
  to.plot <- meta_atlas %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

########################################
## Plot both classes at the same time ##
########################################

# TO-DO....

# Completion token
file.create(file.path(io$outdir,"completed.txt"))

