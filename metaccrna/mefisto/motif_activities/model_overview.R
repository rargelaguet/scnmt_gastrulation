########################
## Plot data overview ##
########################

p <- plot_data_overview(mefisto)

pdf(sprintf("%s/data_overview.pdf",io$outdir), width=6, height=4)
print(p)
dev.off()

#####################
## Plot covariates ##
#####################

p <- ggplot(mefisto@samples_metadata, aes(x=V1, y=V2, fill=lineage10x_2)) +
  geom_point(alpha=0.9, size=2.0, stroke=0.1, shape=21) +
  scale_fill_manual(values=opts$celltype2.colors) +
  # guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(x="", y="") +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank()
  )

pdf(sprintf("%s/umap_by_celltype.pdf",io$outdir), width=6, height=4)
print(p)
dev.off()

mefisto@samples_metadata$foo <- mefisto@samples_metadata$pass_metQC==T | mefisto@samples_metadata$pass_accQC==T
mefisto@samples_metadata$foo[is.na(mefisto@samples_metadata$foo)] <- FALSE
mean(mefisto@samples_metadata$foo)

p <- ggplot(mefisto@samples_metadata, aes(x=V1, y=V2, fill=foo, size=foo, alpha=foo)) +
  geom_point(stroke=0.1, shape=21) +
  scale_size_manual(values = c("TRUE"=1.5, "FALSE"=1.5)) +
  scale_alpha_manual(values = c("TRUE"=1.0, "FALSE"=0.75)) +
  scale_fill_manual(values = c("TRUE"="red", "FALSE"="lightgrey")) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  labs(x="", y="") +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.title = element_blank()
  )
pdf(sprintf("%s/umap_by_logical.pdf",io$outdir), width=6, height=4)
print(p)
dev.off()

#############################
## Plot variance explained ##
#############################

# plot_variance_explained(mefisto, x="view", y="factor")
p <- plot_variance_explained(mefisto, x="view", y="factor", max_r2 = 9)

pdf(sprintf("%s/var_decomposition.pdf",io$outdir), width=6, height=4)
print(p)
dev.off()


##########
## GSEA ##
##########

# Load MsigDb gene set
io$msigFile <- "/Users/ricard/data/genesets/MSigDB/v6.0/mus_musculus/C5/bp_binary_matrix_ensembl.rds"
# io$msigFile <- "/Users/ricard/data/MSigDB/v6.0/mus_musculus/C2/binary_matrix_ensembl.rds"
gene.sets <- readRDS(io$msigFile)

mefisto.gsea <- mefisto
features_names(mefisto.gsea)[["RNA"]] <- toupper(features_names(mefisto.gsea)[["RNA"]])

enrichment.out <- run_enrichment(
  object = mefisto.gsea,
  view = "RNA",
  feature.sets = gene.sets,
  statistical.test = "parametric",
  alpha = 0.01
)
# View(enrichment.out$pval.adj)

p <- plot_enrichment(enrichment.out, factor=4, max.pathways = 15)

pdf(paste0(io$outdir,"/Factor4_enrichment.pdf"), width=10, height=4)
print(p)
dev.off()