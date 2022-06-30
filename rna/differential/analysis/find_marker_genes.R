#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")
io$differential_results <- file.path(io$basedir,"results/rna/differential")
io$outdir <- file.path(io$basedir,"results/rna/differential/marker_genes"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.log2FC <- 2.5
opts$fdr <- 0.01
opts$min.score <- 0.80
opts$group_label <- "stage_lineage3"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & !is.na(celltype)] %>%
  .[,stage_lineage:=paste(stage,celltype,sep="_")] %>%
  .[,stage_lineage2:=paste(stage,celltype2,sep="_")] %>%
  .[,stage_lineage3:=paste(stage,celltype3,sep="_")]

table(sample_metadata$stage_lineage3)

opts$celltypes <- names(which(table(sample_metadata[[opts$group_label]])>=25))

##################
## Load results ##
##################

# i <- "Gut"; j <- "NMP"
dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$differential_results,sprintf("%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2,4)) %>% 
      .[abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  # .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  # .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(dt$celltypeA)),as.character(unique(dt$celltypeB)))) %>% length

#########################
## Define marker genes ##
#########################

foo <- dt[,.(score=sum(direction=="up")), by=c("celltypeA","gene")] %>% setnames("celltypeA","celltype")
bar <- dt[,.(score=sum(direction=="down")), by=c("celltypeB","gene")] %>% setnames("celltypeB","celltype")

markers_genes.dt <- merge(foo,bar,by=c("celltype","gene"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  .[score>=opts$min.score] %>%
  setorder(celltype,-score)
rm(foo,bar)

table(markers_genes.dt$celltype)
length(unique(markers_genes.dt$gene))

# Save
fwrite(markers_genes.dt, file.path(io$outdir,"marker_genes_up.txt.gz"))


################################################
## Plot number of marker genes per cell types ##
################################################

to.plot <- markers_genes.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N") +
  # scale_fill_manual(values=opts$celltype3.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir,"barplot_number_marker_genes.pdf"), width = 8, height = 4)
print(p)
dev.off()

