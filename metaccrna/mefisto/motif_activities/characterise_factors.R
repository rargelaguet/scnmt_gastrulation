################
## Load model ##
################

# source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/analysis/load_model.R")

# mefisto <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/mefisto_model_motifs.rds")
mefisto <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/mefisto_model_motifs_v2.rds")

# mefisto <- subset_factors(mefisto, factors=3:get_dimensions(mefisto)[["K"]])

####################
## Subset factors ##
####################

factors <- names(which(mefisto@cache$variance_explained$r2_per_factor$single_group[,"RNA"]>1))
mefisto <- subset_factors(mefisto, factors)

##################
## Plot factors ##
##################

plot_factor_cor(mefisto)

plot_factor(mefisto, factors = 2, groups="all", color_by = "lineage10x", group_by = "stage", add_violin = T, dodge=T) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factor(mefisto, factors = 3, color_by = "stage", group_by = "stage", add_violin = T, add_boxplot = T, add_dots=F, dodge=T) +
  # scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(mefisto, factors = c(6,7), color_by = "stage_lineage")# +
  # scale_fill_manual(values=opts$celltype.colors)

plot_factors(mefisto, factors = c(1,2), color_by = colMeans(mefisto@data$RNA$single_group))
plot_factors(mefisto, factors = c(6,7), color_by = colMeans(mefisto@data$motif_acc$single_group,na.rm=T))

##################
## Plot weights ##
##################

plot_weights(mefisto, factor = 2, view="motif_met", nfeatures = 5, scale = F)
plot_weights(mefisto, factor = 4, view="motif_acc", nfeatures = 5, scale = F)
plot_top_weights(mefisto, factor = 1, view="motif_met", nfeatures = 10, scale = F)
plot_top_weights(mefisto, factor = 4, view="motif_acc", nfeatures = 15, scale = T)
plot_top_weights(mefisto, factor = 3, view="RNA", nfeatures = 25, scale = F)

#######################
## Cell cycle factor ##
#######################

options("ggrepel.max.overlaps")[[1]] <- 100

options(ggrepel.max.overlaps = Inf)
plot_weights(mefisto, factor = 4, view="RNA", nfeatures = 15, scale = F)

cell_cycle.dt <- fread(io$cell.cycle) 

to.plot <- cell_cycle.dt %>% 
  merge(sample_metadata[,c("id_rna","sample")],by="id_rna") %>%
  merge(get_factors(mefisto, factors=4, as.data.frame=T, scale = T) %>% as.data.table, by="sample") %>%
  .[phase2!="unknown"]

p <- ggboxplot(to.plot, x="phase2", y="value") +
  stat_compare_means() +
  labs(x="Cell cycle phase", y="Factor 4 values")

pdf(sprintf("%s/Factor4_by_cellcycle_phase.pdf",io$outdir), width=6, height=4)
print(p)
dev.off()

#############
## MEFISTO ##
#############

p <- plot_smoothness(mefisto)

pdf(sprintf("%s/smoothness.pdf",io$outdir), width=7, height=2)
print(p)
dev.off()

to.plot <- mefisto@samples_metadata[,c("sample","V1","V2")] %>% as.data.table %>%
  merge(get_factors(mefisto)[[1]] %>% as.data.table(keep.rownames = T) %>% setnames("rn","sample"), by="sample") %>%
  melt(id.vars=c("V1","V2","sample"), variable.name="factor") %>%
  .[,value:=value/max(abs(value)),by="factor"]

# ggscatter(to.plot, x="V1", y="V2", color="value") +
p <- ggscatter(to.plot[factor%in%paste0("Factor",c(1,3,4))], x="V1", y="V2", fill="value", shape=21, stroke=0.15) +
# p <- ggscatter(to.plot[factor%in%paste0("Factor",c(4))], x="V1", y="V2", fill="value", shape=21, stroke=0.15) +
  facet_wrap(~factor) +
  scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank()
  )

pdf(sprintf("%s/umap_by_factor.pdf",io$outdir), width=7, height=3.5)
print(p)
dev.off()


#############
## XXXXXXX ##
#############

# factors <- "all"
# factors <- c(1,3)
# factors <- c(4)

w.met <- get_weights(mefisto, views="motif_met", factor="all", as.data.frame=T) %>% as.data.table %>%
  .[,feature:=gsub("_met","",feature)]
w.acc <- get_weights(mefisto, views="motif_acc", factor="all", as.data.frame=T) %>% as.data.table %>%
  .[,feature:=gsub("_acc","",feature)]

# Scale loadings
w.met[,value:=value/max(abs(value)),by=c("factor")]
w.acc[,value:=value/max(abs(value)),by=c("factor")]

# Merge loadings
w.dt <- merge(
  w.met[,c("feature","factor","value")], 
  w.acc[,c("feature","factor","value")], 
  by = c("feature","factor")
) %>% .[,feature:=strsplit(feature,"_") %>% map_chr(c(1))]

# Scatterplots
for (i in unique(w.dt$factor)) {
  
  to.plot <- w.dt[factor==i]
  to.label <- w.dt %>% 
    .[factor==i] %>%
    .[,value:=abs(value.x)+abs(value.y)] %>% setorder(-value) %>% head(n=15)
  
  p <- ggpubr::ggscatter(to.plot, x="value.x", y="value.y", size=1.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) +
    scale_x_continuous(breaks=c(-1,0,1)) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    ggrepel::geom_text_repel(data=to.label, aes(x=value.x, y=value.y, label=feature), size=3,  max.overlaps=100) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    stat_cor(method = "pearson") +
    labs(x="Methylation weights", y="Accessibility weights")
  
  pdf(sprintf("%s/%s_met_vs_acc_weights.pdf",io$outdir,i), width=6.5, height=5)
  print(p)
  dev.off()
}

