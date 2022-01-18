here::i_am("rna/celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--stages',         type="character",       nargs="+",   help='Stages')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
# args$stages <- opts$stages
# args$celltype_label <- "celltype2"
# args$outdir <- file.path(io$basedir,sprintf("results/rna/celltype_proportions/%s",args$celltype_label))
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & stage%in%args$stages & !is.na(eval(as.name(args$celltype_label)))]

# if (!"stage"%in%colnames(sample_metadata)) {
#   sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]
# }
# sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]

table(sample_metadata$stage)

######################
## Stacked barplots ##
######################

# p <- ggplot(to.plot, aes(x=sample, y=celltype_proportion)) +
#   geom_bar(aes(fill=celltype), stat="identity", color="black") +
#   facet_wrap(~class, scales = "free_x", nrow=1) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title = element_blank(),
#     axis.text.y = element_blank(),
#     # axis.text.x = element_text(color="black", size=rel(0.75)),
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank(),
#     axis.line = element_blank()
#   )
# 
# pdf(file.path(args$outdir,"celltype_proportions_stacked_barplots.pdf"), width=7, height=5)
# print(p)
# dev.off()

#####################################
## Calculate cell type proportions ##
#####################################

celltype_proportions.dt <- sample_metadata %>%
  .[,N:=.N,by="stage"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("stage",args$celltype_label)] %>%
  setorder(stage)  %>% .[,stage:=factor(stage,levels=opts$stages)]

##########
## Plot ##
##########

# Define colours and cell type order
if (args$celltype_label=="celltype") {
  opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(celltype_proportions.dt$celltype)]
} else if (args$celltype_label=="celltype2") {
  opts$celltype.colors <- opts$celltype2.colors[names(opts$celltype2.colors) %in% unique(celltype_proportions.dt$celltype2)]
} else if (args$celltype_label=="celltype3") {
  opts$celltype.colors <- opts$celltype3.colors[names(opts$celltype3.colors) %in% unique(celltype_proportions.dt$celltype3)]
}
celltype_proportions.dt[[args$celltype_label]] <- factor(celltype_proportions.dt[[args$celltype_label]], levels=rev(names(opts$celltype.colors)))

p <- ggplot(celltype_proportions.dt, aes_string(x=args$celltype_label, y="N")) +
  geom_bar(aes_string(fill=args$celltype_label), stat="identity", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~stage, nrow=1, scales="fixed") +
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

pdf(file.path(args$outdir,"celltype_proportions_per_stage.pdf"), width=10, height=5)
print(p)
dev.off()
