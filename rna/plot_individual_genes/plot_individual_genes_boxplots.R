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
  "E6.5_Mesoderm",
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm",
  "E7.5_ExE_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,celltype,sep="_")] %>%
  .[,stage_lineage2:=paste(stage,celltype2,sep="_")] %>%
  .[,stage_lineage3:=paste(stage,celltype3,sep="_")] %>%
  .[pass_rnaQC==T & stage_lineage3%in%opts$stage_lineage] %>%
  .[,stage_lineage3:=factor(stage_lineage3, levels=opts$stage_lineage)] 

stopifnot(unique(sample_metadata$stage_lineage3) %in% opts$stage_lineage)

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

##########
## Plot ##
##########

# genes.to.plot <- fread("/Users/argelagr/data/scnmt_gastrulation_argelaguet2019/results/rna/differential/marker_genes/marker_genes_up.txt.gz") %>% .[celltype=="E4.5_Epiblast",gene] %>% unique
# genes.to.plot <- foo[!is.na(logFC)] %>% setorder(logFC) %>% head(n=30) %>% .$gene
genes.to.plot <-  c("Pou5f1", "Zfp42", "Utf1", "Pim2", "Slc7a3", "Dppa5a", "Fgf5", "Dppa3", "Dppa4", "Tfap2c", "Nanog")
genes.to.plot <-  c("L1td1", "Dnmt3b", "Snrpn", "Gng3", "Pdzd4")
# genes.to.plot <- grep("Tet", rownames(sce),value=T)
genes.to.plot <-  c("Nanog","Pou5f1")

stopifnot(genes.to.plot %in% rownames(sce))

# genes.to.plot <- foo[logFC>6] %>% setorder(-logFC) %>% head(n=25) %>% .$gene
for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))

  to.plot <- data.table(
    id_rna = colnames(sce),
    expr = logcounts(sce)[gene,]
  ) %>% merge(sample_metadata, by="id_rna")# %>% .[,N:=.N,by=c("sample","stage_lineage")] %>% .[N>=5]
  
  p <- ggplot(to.plot, aes(x=celltype3, y=expr, fill=celltype3)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
    geom_jitter(shape=21, size=0.75, alpha=0.35, width=0.05, stroke=0.15) +
    scale_fill_manual(values=opts$celltype3.colors, drop=F) +
    # stat_summary(fun.data = give.n, geom = "text", size=2.5) +
    # facet_wrap(~stage, nrow=1, scales="free_x") +
    facet_grid(~stage, scales = "free_x", space='free') +
    theme_classic() +
    labs(title=gene, x="",y=sprintf("%s expression",gene)) +
    guides(x = guide_axis(angle = 90)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      plot.title = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(0.95)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_boxplots_single_cells_argelaguet2019.pdf",io$outdir,gene), width=8, height=3.5)
  print(p)
  dev.off()
}


#######################
## Plot marker genes ##
#######################

marker_genes.dt <- fread(file.path(io$basedir,"results/rna/differential/marker_genes/germ_layers/marker_genes_up.txt.gz")) %>%
  .[celltype%in%opts$stage_lineage & !grepl("Rik",gene)]
celltypes.to.plot <- unique(marker_genes.dt$celltype)
# celltypes.to.plot <- c("E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Mesoderm")

# i <- "E7.5_Mesoderm"; j <- "Mesp1"
for (i in celltypes.to.plot) {
  genes.to.plot <- marker_genes.dt[celltype==i,gene]
  genes.to.plot <- genes.to.plot[genes.to.plot %in% rownames(sce)]
  for (j in genes.to.plot) {
    
    to.plot <- data.table(
      id_rna = colnames(sce),
      expr = logcounts(sce)[j,],
      celltype3 = sce$celltype3,
      stage = sce$stage
    ) 
    
    p <- ggplot(to.plot, aes(x=celltype3, y=expr, fill=celltype3)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
      geom_jitter(shape=21, size=1, alpha=0.5, width=0.05) +
      scale_fill_manual(values=opts$celltype3.colors, drop=F) +
      stat_summary(fun.data = give.n, geom = "text", size=2.5) +
      # facet_wrap(~stage, nrow=1, scales="free_x") +
      facet_grid(~stage, scales = "free_x", space='free') +
      theme_classic() +
      labs(title=j, x="",y=sprintf("%s expression",j)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="black",size=rel(0.95)),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(1.0)),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
    pdf(sprintf("%s/marker_genes/%s_%s.pdf",io$outdir,i,j), width=8, height=3.5)
    print(p)
    dev.off()
  }
}
