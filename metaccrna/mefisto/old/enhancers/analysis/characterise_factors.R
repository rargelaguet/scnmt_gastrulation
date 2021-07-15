################
## Load model ##
################

source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/analysis/load_model.R")

#####################
## Plot covariates ##
#####################

covariates.dt <- get_covariates(mefisto, as.data.frame = T) %>% 
  as.data.table %>% dcast(sample~covariate) %>%
  merge(sample_metadata, by="sample")

p1 <- ggplot(covariates.dt, aes(x=V1, y=V2, color=lineage10x_2)) +
  geom_point(alpha=0.7, size=2.0) +
  scale_color_manual(values=opts$celltype2.colors) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(x="", y="") +
  theme_classic() +
  ggplot_theme_NoAxes()

pdf(sprintf("%s/MOFA_%s.pdf",io$outdir,algorithm), width=5, height=7, useDingbats = F)
print(p1)
dev.off()


#############################
## Plot variance explained ##
#############################

plot_variance_explained(mefisto, x="view", y="factor", max_r2 = 5)

##################
## Plot factors ##
##################

plot_factor_cor(mefisto)

plot_factor(mefisto, factors = 3, groups="all", color_by = "stage_lineage", group_by = "stage", add_violin = T, dodge=T) +
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

plot_factors(mefisto, factors = c(1,2), color_by = "stage")# +
  # scale_fill_manual(values=opts$celltype.colors)

plot_factors(mefisto, factors = c(1,2), color_by = colMeans(mefisto@data$RNA$single_group))


##########
## UMAP ##
##########

mefisto <- run_umap(mefisto, factors = "all")
# mefisto <- run_umap(mefisto, factors = 2:get_dimensions(mefisto)[["K"]])
plot_dimred(mefisto, method="UMAP", color_by = "stage", stroke=0.1)# +
  # scale_fill_manual(values=opts$celltype.colors)

##################
## Plot weights ##
##################

plot_weights(mefisto, factor = 9, nfeatures = 10, scale = F)

###############
## Plot data ##
###############

plot_data_heatmap(mefisto,
                  factor = 1,
                  features = 25,
                  view = 1,
                  denoise = FALSE,
                  legend = TRUE,
                  # min.value = 0, max.value = 6,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  scale="row",
                  annotation_samples = "Category",  annotation_colors = list("Category"=opts$colors), annotation_legend = F
)

plot_data_scatter(mefisto, factor=1, view=1, color_by = "embryo", features = 8, dot_size = 2)

#############
## MEFISTO ##
#############

plot_smoothness(mefisto)

to.plot <- mefisto@samples_metadata[,c("sample","V1","V2")] %>% as.data.table %>%
  merge(get_factors(mefisto)[[1]] %>% as.data.table(keep.rownames = T) %>% setnames("rn","sample"), by="sample") %>%
  melt(id.vars=c("V1","V2","sample"), variable.name="factor") %>%
  .[,value:=value/max(abs(value)),by="factor"]

ggscatter(to.plot, x="V1", y="V2", color="value") +
  facet_wrap(~factor) +
  scale_color_gradient2(low = "gray50", mid="gray90", high = "red")

compare_factors(list(mefisto,mefisto.old))
