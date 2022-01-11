here::i_am("rna/mapping/analysis/compare_mappings.R")

# load default setings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata1',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--metadata2',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################


## START TEST ##
args$metadata1 <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$metadata2 <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping_all_samples.txt.gz")
args$outdir <- file.path(io$basedir,"results_new/rna/mapping/pdf/comparison")
## END TEST ##

dir.create(args$outdir, showWarnings = F)

#####################
## Load metadata 1 ##
#####################

sample_metadata1.dt <- fread(args$metadata1) %>%
  .[pass_rnaQC==TRUE & !is.na(celltype.mapped)] %>% 
  .[,sample_plate:=sprintf("%s_%s",sample,plate)] %>%
  .[,c("cell","sample","plate","celltype.mapped","celltype.score")] %>%
  setnames("celltype.mapped","celltype") %>%
  .[,class:="Per sample"]

#####################
## Load metadata 2 ##
#####################

sample_metadata2.dt <- fread(args$metadata2) %>%
  .[pass_rnaQC==TRUE & !is.na(celltype.mapped)] %>% 
  .[,sample_plate:=sprintf("%s_%s",sample,plate)] %>%
  .[,c("cell","sample","plate","celltype.mapped","celltype.score")] %>%
  setnames("celltype.mapped","celltype") %>%
  .[,class:="All samples"]

###########
## Merge ##
###########

# sample_metadata_merged.dt <- merge(sample_metadata1.dt, sample_metadata2.dt, by=c("cell"))
sample_metadata_cat.dt <- rbind(sample_metadata1.dt, sample_metadata2.dt)

to.plot <- table(sample_metadata_merged.dt$celltype.mapped.x, sample_metadata_merged.dt$celltype.mapped.y)

##########################################
## Plot celltype proportions per sample ##
##########################################

# Calculate cell type proportions per sample
celltype_proportions.dt <- sample_metadata_cat.dt %>%
  .[,.(N=.N),by=c("sample","celltype","class")]

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(celltype_proportions.dt$celltype)]
celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

samples.to.plot <- unique(celltype_proportions.dt$sample)

for (i in samples.to.plot) {
  
  to.plot <- celltype_proportions.dt[sample==i]
  
  p <- ggplot(to.plot, aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/celltype_proportions_%s_comparison.pdf",args$outdir,i), width=7, height=5)
  print(p)
  dev.off()
}

