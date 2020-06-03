source("/Users/ricard/scnmt_gastrulation/settings.R")
source("/Users/ricard/scnmt_gastrulation/metaccrna/dynamics_individual_examples/load_data.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/metaccrna/plot_dynamics_individual_examples")

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
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
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
opts$min.cpg <- 3
opts$min.gpc <- 10

# Define cells to use
opts$met_cells <- sample_metadata %>% .[pass_metQC==T & stage_lineage%in%opts$stage_lineage,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage,id_rna]
opts$acc_cells <- sample_metadata %>% .[pass_accQC==T & stage_lineage%in%opts$stage_lineage,id_acc]

###############
## Load data ##
###############

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[,c("sample","id_met","id_rna","id_acc","stage","stage_lineage","lineage10x_2")] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ]

# Load the three omics
dt <- load_data(io, rna.id, met.id, met.anno, acc.id, acc.anno, opts$min.cpg, opts$min.gpc)

# Merge data with sample metadata
dt$acc <- merge(dt$acc, sample_metadata, by="id_acc")
dt$met <- merge(dt$met, sample_metadata, by="id_met")
dt$rna <- merge(dt$rna, sample_metadata, by="id_rna")

# bind in a single data table
dt_all <- do.call("rbind",list(
  dt$rna[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="RNA expression"],
  dt$met[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="DNA methylation"],
  dt$acc[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="Chromatin accessibility"]
)) %>% .[,stage_lineage:=gsub("_"," ",stage_lineage)] %>% 
  .[,assay:=factor(assay,levels=c("RNA expression","DNA methylation","Chromatin accessibility"))]
  
#######################
## Generate Boxplots ##
#######################

p <- ggplot(dt_all, aes(x=stage, y=value)) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  geom_jitter(aes(color=assay), size=0.5) +
  geom_violin(aes(fill=assay), alpha=0.5, size=0.25) +
  geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
  scale_fill_manual(values=opts$color) +
  scale_color_manual(values=opts$color) +
  labs(x="", y="", title="") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(colour="black", size=rel(1.1), vjust=1.5),
    axis.text.x = element_text(size=rel(1.2), color="black"),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.y = element_line(colour="black"),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="none"
  )
print(p)

##########
## Save ##
##########

# pdf(sprintf("%s/boxplot_rna%s_met%s_acc%s.pdf",io$outdir,rna.id,met.id,acc.id), useDingbats=F, width=6, height=5)
print(p)
# dev.off()
