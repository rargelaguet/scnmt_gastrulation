here::i_am("rna/plot_individual_genes/plot_individual_genes_boxplots.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################
# I/O
io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/individual_genes/fig"); dir.create(io$outdir, showWarnings = F)

# Define which stages to use
# opts$stages <- c("E3.5", "E4.5", "E5.5", "E6.5", "E7.5", "E8.5")
opts$stages <- c("E3.5", "E4.5", "E5.5", "E6.5")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==T & stage%in%opts$stages] %>%
  .[,stage:=factor(stage,levels=opts$stages)]

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$id_rna, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("id_rna") %>% DataFrame

########################################
## Plot one gene at a time, per stage ##
########################################

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
genes.to.plot <- grep("Tet", rownames(sce),value=T)

for (i in genes.to.plot) {
  
    to.plot <- data.table(
      id_rna = colnames(sce),
      expr = logcounts(sce)[i,]
    ) %>% merge(sample_metadata, by="id_rna")# %>% .[,N:=.N,by=c("sample","stage_lineage")] %>% .[N>=5]
    
    give_n <- function(x){ return(c(y = max(to.plot$expr), label = length(x))) }
    
    p <- ggplot(to.plot, aes(x=stage, y=expr, fill=stage)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
      geom_jitter(shape=21, size=1, alpha=0.5, width=0.05, stroke=0.15) +
      scale_fill_manual(values=opts$stage.colors, drop=F) +
      stat_summary(fun.data = give_n, geom = "text", size=2.5) +
      # facet_wrap(~stage, nrow=1, scales="free_x") +
      theme_classic() +
      labs(title=i, x="",y=sprintf("%s expression",i)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="black",size=rel(1.25)),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(1.0)),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(sprintf("%s/%s.pdf",io$outdir,i), width=10, height=5)
      print(p)
      dev.off()
}

#####################################################
## Plot multiple genes at the same time, per stage ##
#####################################################

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
# genes.to.plot <- grep("Tet", rownames(sce),value=T)
genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b")

to.plot <- logcounts(sce[genes.to.plot]) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars=("gene"), variable.name="id_rna", value.name = "expr") %>% 
  merge(sample_metadata, by="id_rna")
  
p <- ggplot(to.plot, aes(x=stage, y=expr, fill=stage)) +
  geom_violin(alpha=0.30) +
  geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
  # geom_jitter(shape=21, size=0.5, alpha=0.4, width=0.03, stroke=0.05) +
  ggrastr::geom_jitter_rast(shape=21, size=0.5, alpha=0.4, width=0.03, stroke=0.05) +
  scale_fill_manual(values=opts$stage.colors, drop=F) +
  facet_wrap(~gene, nrow=1, scales="free_x") +
  theme_classic() +
  labs(x="", y="Gene expression") +
  theme(
    strip.text = element_text(size=rel(1.25)),
    strip.background = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(1)),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )
  
# pdf(sprintf("%s/tet_genes.pdf",io$outdir), width=7.5, height=4)
pdf(sprintf("%s/dnmt_genes.pdf",io$outdir), width=7.5, height=3.5)
print(p)
dev.off()


############################################
## Plot one gene at a time, per cell type ##
############################################

# Define which cells to use
opts$stage_lineage <- c(
  "E3.5_ICM",
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  "E5.5_Epiblast",
  "E5.5_ExE_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_ExE_endoderm",
  "E6.5_ExE_ectoderm",
  "E6.5_Mesoderm"
  # "E7.5_Epiblast",
  # "E7.5_Primitive_Streak",
  # "E7.5_Ectoderm",
  # "E7.5_Endoderm",
  # "E7.5_Mesoderm"
  # "E7.5_ExE_endoderm"
)

sample_metadata_filt <- sample_metadata %>% 
  .[,stage_lineage3:=paste(stage,celltype3,sep="_")] %>%
  .[stage_lineage3%in%opts$stage_lineage]

sce_filt <- sce[,sample_metadata_filt$id_rna]

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
genes.to.plot <- c(grep("Tet|Dnmt", rownames(sce_filt),value=T),"Uhrf1")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    id_rna = colnames(sce_filt),
    expr = logcounts(sce_filt)[i,]
  ) %>% merge(sample_metadata_filt, by="id_rna")# %>% .[,N:=.N,by=c("sample","stage_lineage")] %>% .[N>=5]
  
  give_n <- function(x){ return(c(y = max(to.plot$expr), label = length(x))) }
  
  p <- ggplot(to.plot, aes(x=celltype3, y=expr, fill=celltype3)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
    geom_jitter(shape=21, size=1, alpha=0.5, width=0.05, stroke=0.15) +
    scale_fill_manual(values=opts$celltype3.colors, drop=F) +
    stat_summary(fun.data = give_n, geom = "text", size=2.5) +
    facet_grid(~stage, scales = "free_x", space='free') +
    theme_classic() +
    labs(title=i, x="",y=sprintf("%s expression",i)) +
    guides(x = guide_axis(angle = 90)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      plot.title = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(1.25)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_per_celltype.pdf",io$outdir,i), width=10, height=5)
  print(p)
  dev.off()
}
