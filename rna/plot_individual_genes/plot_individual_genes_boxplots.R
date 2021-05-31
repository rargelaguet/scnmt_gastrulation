#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/utils.R")
}

## I/O ##
io$outdir <- paste0(io$basedir,"/rna/results/individual_genes")

## Define options ##

# Define which cells to use
opts$stage_lineage <- c(
  "E3.5_ICM",
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Visceral_endoderm",
  "E6.5_Mesoderm",
  "E7.5_Epiblast",
  "E7.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage] %>%
  .[,stage_lineage:=factor(stage_lineage, levels=opts$stage_lineage)] 

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$id_rna, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Rename genes
rownames(sce) <- rowData(sce)$symbol

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("id_rna") %>% DataFrame

################
## Parse data ##
################

##########
## Plot ##
##########

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")

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
      
      p <- ggplot(to.plot, aes(x=lineage10x_2, y=expr, fill=lineage10x_2)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
        scale_fill_manual(values=opts$celltype2.colors, drop=F) +
        stat_summary(fun.data = give.n, geom = "text", size=2.5) +
        facet_wrap(~stage, nrow=1, scales="free_x") +
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

