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
io$outdir <- paste0(io$basedir,"/rna/results/individual_genes/umap")

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
## Load UMAP ##
###############

umap.dt <- fread(io$umap)

################
## Parse data ##
################

cells <- intersect(umap.dt$sample,sample_metadata$sample)

umap.dt <- umap.dt %>% .[sample%in%cells]
sample_metadata <- sample_metadata %>% .[sample%in%cells]

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

##########
## Plot ##
##########

genes.to.plot <- c("Foxa1","Foxa2","Hnf1a","Hnf1b","Foxd1","Foxd2")
# genes.to.plot <- c("Hoxc11")

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s_umap.pdf",io$outdir,gene)
    
    if (!file.exists(outfile)) {
      
      to.plot <- data.table(
        id_rna = colnames(sce),
        expr = logcounts(sce)[gene,]
      ) %>% merge(sample_metadata, by="id_rna")  %>%
        merge(umap.dt,by="sample")
      
      p <- ggplot(to.plot, aes(x=V1, y=V2, fill=expr)) +
        # geom_point(alpha=0.9, size=3.0, shape=21) +
        ggrastr::geom_point_rast(alpha=0.9, size=3.0, shape=21) +
        scale_fill_gradient(low = "white", high = "darkgreen") +
        # scale_fill_gradientn(colours = terrain.colors(10)) +
        theme_classic() +
        ggplot_theme_NoAxes() +
        theme(
          legend.position = "right"
        )
      
      pdf(outfile, width=4.5, height=6)
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

