library(ggplot2)
library(RColorBrewer)
library(data.table)
library(purrr)
source("/Users/ricard/gastrulation/metaccrna/bubble_plots/load_data.R")

scale <- function(value, min.data, max.data, min.scaled, max.scaled) {
  stopifnot(is.numeric(value))
  stopifnot(value<=max.data & value>=min.data)
  return ((max.scaled - min.scaled) * (value - min.data) / (max.data - min.data)) + min.scaled
}

################
## Define I/O ##
################

io <- list()
io$outdir <- "/Users/ricard/data/gastrulation_norsync_stuff/metaccrna/dynamics"
io$sample.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"

####################
## Define options ##
####################

opts <- list()

# Define gene identity
rna.id <- "ENSMUSG00000021944"
met.id <- "ENSMUSG00000021944"
met.anno <- "prom_2000_2000"
acc.id <- "ENSMUSG00000021944"
acc.anno <- "prom_2000_2000"


# Define stages to plot
opts$stage_lineage <- c(
  
  # # # E4.5
  # "E4.5_Epiblast",
  # 
  # # E5.5
  # "E5.5_Epiblast",
  # 
  # # E6.5
  # "E6.5_Epiblast",
  # "E6.5_Primitive_Streak",
  # "E6.5_Mesoderm",
  
  # E7.5
  "E7.5_Epiblast",
  # "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Mesoderm"
)

# Define colors for the omics
opts$color <- c(
  "RNA expression"="#3CB54E",
  "Chromatin accessibility"="#00BFC4",
  "DNA methylation"="#F37A71"
)

# Define minimum coverage
opts$min.cpg <- 1
opts$min.gpc <- 5

# Define cells to use
tmp <- fread(io$sample.metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]
opts$met_cells <- tmp %>% .[pass_metQC==T & stage_lineage%in%opts$stage_lineage,id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage,id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T & stage_lineage%in%opts$stage_lineage,id_acc]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample.metadata,stringsAsFactors=F) %>%
  .[,c("sample","id_met","id_rna","id_acc","stage","lineage10x","lineage10x_2")] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ] %>%
  droplevels()

sample_metadata %>%
  # .[lineage10x_2=="Endoderm",lineage10x_2:=ifelse(lineage10x=="Notochord","Endoderm (notochord)","Endoderm (not notochord)")] %>%
  # .[lineage10x_2=="Mesoderm",lineage10x_2:=ifelse(lineage10x=="Nascent_mesoderm","Nascent mesoderm","Mature mesoderm")] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]

# Load the three omics
dt <- load_data(rna.id, met.id, met.anno, acc.id, acc.anno, opts$min.cpg, opts$min.gpc)

# Merge data with sample metadata
dt$acc <- merge(dt$acc, sample_metadata, by="id_acc") %>% droplevels()
dt$met <- merge(dt$met, sample_metadata, by="id_met") %>% droplevels()
dt$rna <- merge(dt$rna, sample_metadata, by="id_rna") %>% droplevels()

# bind in a single data table
dt_all <- do.call("rbind",list(
  dt$rna[,c("sample","gene","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="RNA expression"] %>% setnames("gene","id"),
  dt$met[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="DNA methylation"],
  dt$acc[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="Chromatin accessibility"]
))

dt_all[,stage_lineage:=gsub("_"," ",stage_lineage)]

dt_all[,assay:=factor(assay,levels=c("RNA expression","DNA methylation","Chromatin accessibility"))]
  
#######################
## Generate Boxplots ##
#######################

p <- ggplot(dt_all, aes(x=lineage10x_2, y=value)) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  geom_jitter(aes(color=assay), size=0.5) +
  geom_violin(aes(fill=assay), alpha=0.5, size=0.25) +
  geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
  scale_fill_manual(values=opts$color) +
  scale_color_manual(values=opts$color) +
  labs(x="", y="", title="") +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(colour="black", size=rel(1.1), vjust=1.5),
    # axis.text.x = element_text(size=rel(1.2), angle=30, hjust=1, vjust=1, color="black"),
    axis.text.x = element_text(size=rel(1.2), color="black"),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.3)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.y = element_line(colour="black"),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    # strip.text = element_text(size=rel(1.3), color="black"),
    strip.text = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position="none",
    legend.text=element_text(size=15),
    legend.title=element_blank(),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
print(p)

##########
## Save ##
##########

pdf(sprintf("%s/mofaplus/boxplot_rna%s_met%s_acc%s.pdf",io$outdir,rna.id,met.id,acc.id), useDingbats=F, width=6, height=5)
print(p)
dev.off()
