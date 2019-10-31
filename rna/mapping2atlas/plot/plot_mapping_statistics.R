#######################################
## Script to plot mapping statistics ##
#######################################

library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
io$sample_metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$mapping.results <- "/Users/ricard/data/gastrulation/rna/mapping_10x/mapping.rds"
io$outdir <- "/Users/ricard/data/gastrulation/rna/mapping_10x/pdf"

###################
# Define options ##
###################

source("/Users/ricard/gastrulation_clean/rna/mapping2atlas/plot/utils.R")

#########################
# Load sample metadata ##
#########################

sample_metadata <- fread(io$sample_metadata) %>%
  .[!is.na(lineage10x_2) ]


to.plot <- sample_metadata[,.N, by=c("stage","lineage10x_2")] %>% 
  .[, lineage10x_2:=stringr::str_replace_all( lineage10x_2,"_"," ")] %>%
  .[, lineage10x_2:=factor( lineage10x_2,levels=names(colors_lineage10x_2))]
  
p <- ggplot(to.plot, aes(x= lineage10x_2, y=N)) +
  geom_bar(aes(fill= lineage10x_2), stat="identity", color="black") +
  scale_fill_manual(values=colors_lineage10x_2) +
  facet_wrap(~stage, nrow=1) +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.3)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

pdf(paste0(io$outdir,"/mapping_stats/mapping_stats.pdf"), width=11, height=5)
print(p)
dev.off()
