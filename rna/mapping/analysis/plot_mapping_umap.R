here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$query_metadata <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
args$outdir <- file.path(io$basedir,"results_new/rna/mapping/pdf")
## END TEST ##

dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

#####################
## Define settings ##
#####################

# Options

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

opts$remove_ExE_cells <- TRUE
opts$subset_atlas_cells <- TRUE

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(args$query_metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(closest.cell)]
  .[pass_rnaQC==TRUE & !is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

if (opts$remove_ExE_cells) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

if (opts$remove_ExE_cells) {
  meta_atlas <- meta_atlas %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# Subset cells to speed up plotting
if (opts$subset_atlas_cells) {
  umap.dt <- umap.dt %>% .[sample.int(n=nrow(.), size=nrow(.)/4)]
}

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


####################
## Plot all cells ##
####################

# to.plot <- umap.dt %>% copy %>%
#   .[,index:=match(cell, sample_metadata[,closest.cell] )] %>% 
#   .[,mapped:=as.factor(!is.na(index))] %>% 
#   .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas cells","Query cells"))] %>%
#   setorder(mapped) 
# 
# p <- plot.dimred(to.plot, query.label = "Query cells", atlas.label = "Atlas cells")
# 
# pdf(sprintf("%s/umap_mapped_allcells.pdf",args$outdir), width=8, height=6.5)
# print(p)
# dev.off()

###############################
## Plot one sample at a time ##
###############################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

###############################
## Plot one class at a time ##
###############################

classes.to.plot <- unique(sample_metadata$class)

for (i in classes.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))



##########
## TEST ##
##########

samples.to.plot <- c(
  "E7.5_WT",
  "E7.5_TET_TKO",
  # "E7.5_TET_TKO_crispr",
  "E8.5_WT",
  "E8.5_WT_CD41+",
  "E8.5_WT_KDR+",
  "E8.5_WT_KDR+_CD41+",
  "E8.5_TET_TKO",
  "E8.5_TET_TKO_KDR+",
  "E8.5_TET_TKO_CD41+",
  "E8.5_TET_TKO_KDR+_CD41+"
)
# samples.to.plot <- unique(sample_metadata$sample)

to.plot <- samples.to.plot %>% map(function(i) {
  umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
    .[,sample:=i] %>%
    setorder(mapped)
}) %>% rbindlist %>% .[,sample:=factor(sample,levels=samples.to.plot)]

p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") + 
  facet_wrap(~sample, ncol=2) +
  theme(legend.position = "none") 

pdf(file.path(args$outdir,"per_class/umap_mapped_all_classes.pdf"), width=6, height=12)
print(p)
dev.off()