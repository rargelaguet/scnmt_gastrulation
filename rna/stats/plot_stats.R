library(data.table)
library(purrr)
library(ggplot2)
library(RColorBrewer)

#####################
## Define settings ##
#####################

io <- list()
io$stats <- "/Users/ricard/gastrulation/rna/stats/out/rna_stats.txt"
io$metadata.file <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$outdir <- "/Users/ricard/gastrulation/rna/stats/out"

###############
## Load data ##
###############

# Load RNA statistics
stats <- fread(io$stats)

# Load sample metadata
sample_metadata <- fread(io$metadata.file) 

################
## Parse data ##
################

to.plot <- stats %>% merge(sample_metadata, by="id_rna") %>% 
  .[,stage:=factor(stage,levels=c("E4.5","E5.5","E6.5","E7.5"))] %>%
  .[,stage_plate:=paste(stage,plate,sep=" ")] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep=" ")] %>%
  .[,log_counts:=log(total_counts)]

p <- ggboxplot(to.plot, x = "stage_plate", y = "log_counts", fill="#3CB54E", outlier.shape=NA) +
  labs(x="", y="Library size (log)") +
  # coord_cartesian(ylim=c(0,3e6)) +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p

pdf(paste0(io$outdir,"/rna_stats.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()