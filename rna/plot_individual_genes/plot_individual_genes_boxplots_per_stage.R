here::i_am("rna/plot_individual_genes/plot_individual_genes_boxplots.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## I/O ##
io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/individual_genes"); dir.create(io$outdir, showWarnings = F)

## Define options ##

# Define which stages to use
opts$stages <- c("E3.5", "E4.5", "E5.5", "E6.5", "E7.5", "E8.5")

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

#############################
## Plot one gene at a time ##
#############################

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
genes.to.plot <- grep("Tet", rownames(sce),value=T)

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s.pdf",io$outdir,gene)
    
    if (!file.exists(outfile)) {
      
      to.plot <- data.table(
        id_rna = colnames(sce),
        expr = logcounts(sce)[gene,]
      ) %>% merge(sample_metadata, by="id_rna")# %>% .[,N:=.N,by=c("sample","stage_lineage")] %>% .[N>=5]
      
      # to.plot <- to.plot[stage_lineage%in%celltypes.tmp]
      
      p <- ggplot(to.plot, aes(x=stage, y=expr, fill=stage)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
        geom_jitter(shape=21, size=1, alpha=0.5, width=0.05) +
        scale_fill_manual(values=opts$stage.colors, drop=F) +
        stat_summary(fun.data = give.n, geom = "text", size=2.5) +
        # facet_wrap(~stage, nrow=1, scales="free_x") +
        theme_classic() +
        labs(title=gene, x="",y=sprintf("%s expression",gene)) +
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
      
        # pdf(outfile, width=11, height=10)
        pdf(outfile, width=10, height=5)
        # png(outfile, width = 1100, height = 1000)
        print(p)
        dev.off()
        
    } else {
      print(sprintf("%s already exists...",outfile))
    }

  } else {
    print(sprintf("%s not found",gene))
  }
}



##########################################
## Plot multiple genes at the same time ##
##########################################

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
pdf(sprintf("%s/dnmt_genes.pdf",io$outdir), width=7.5, height=4)
print(p)
dev.off()

