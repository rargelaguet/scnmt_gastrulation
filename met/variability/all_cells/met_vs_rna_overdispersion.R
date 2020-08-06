
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}
io$met.overdispersion <- paste0(io$scmet,"/betabinomial/bayes/txt/hvf.txt.gz")
io$rna.overdispersion <- paste0(io$basedir,"/rna/results/variability/gene_variability.txt.gz")
io$outdir <- paste0(io$basedir,"/metrna/overdispersion")


opts$anno <- c(
  "prom_2000_2000"
)

###############
## Load data ##
###############

# Load methylation overdispersion estimates
# io$met.overdispersion <- "/Users/ricard/Downloads/filt_prom_2000_2000_vb.txt.gz"
# met_overdispersion.dt <- fread(io$met.overdispersion) %>%
#   .[anno%in%opts$anno]# %>%
#   # setnames("feature_name","ens_id")
# 
# met_overdispersion.dt %>% setnames("Feature","ens_id") %>% setnames("epsilon_median","epsilon")
# 
# # Load RNA overdispersion estimates
# rna_overdispersion.dt <- fread(io$rna.overdispersion)
# 
# # Load gene metadata
# gene_metadata <- fread(io$gene.metadata) %>% .[,c("ens_id","symbol")]

##########
## Link ##
##########


dt <- fread("/Users/ricard/data/gastrulation/metrna/overdispersion/metrna_variability_basics.csv")
# dt <- merge(
#   met_overdispersion.dt,
#   rna_overdispersion.dt
# ) %>% merge(gene_metadata,by="ens_id")

# Define groups
dt[,color:="black"]
dt[epsilon>0.65 & rna_epsilon>1.23,color:="red"]
dt[epsilon<0.2 & rna_epsilon>2.4,color:="blue"]

# dt[,color:="black"]
# dt[epsilon>0.75 & bio.disp>3,color:="red"]
# dt[epsilon<0.25& bio.disp>4,color:="blue"]

##########
## Plot ##
##########

# genes.to.highlight <- c(
#   "Pim2","Slc7a3","Cldn6","Lefty2","Trh","Aplnr","Zic3","Dll3","Krt18","Krt8","Id3","Mesp1","Mixl1",
#   "Dppa5a","Dppa4","Peg3","Cbx7","Slc38a4","Tet1","Tex19.1","Dppa2","Spp1","Mest"
# )

p <- ggplot(dt, aes(x=epsilon, y=rna_epsilon, alpha=color, fill=color, size=color)) +
  labs(x="Met. overdispersion", y="RNA overdispersion") +
  geom_point(shape=21, stroke=0.15, color="black") +
  ggrastr::geom_point_rast(shape=21, stroke=0.15, color="black") +
  # geom_point() +
  scale_fill_manual(values=c("black"="gray80", "red"="#4E934D", "blue"="#C400AD")) +
  scale_size_manual(values=c("black"=0.5, "red"=2, "blue"=2)) +
  scale_alpha_manual(values=c("black"=0.7, "red"=0.9, "blue"=0.9)) +
  ggrepel::geom_text_repel(data=dt[color=="red"], aes(x=epsilon, y=rna_epsilon, label=symbol), size=5, color="#4E934D") +
  ggrepel::geom_text_repel(data=dt[color=="blue"], aes(x=epsilon, y=rna_epsilon, label=symbol), size=5, color="#C400AD") +
  coord_cartesian(xlim=c(-0.9,1.55)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/met_vs_rna_overdispersion_scatterplot.pdf"), width=11, height=4, useDingbats = F)
print(p)
dev.off()


##############################
## Plot Met/RNA correlation ##
##############################

# dt.cor <- fread("/Users/ricard/data/gastrulation/metrna/cor/metrna_cor_promoters_v2.txt.gz") %>%
#   .[,c("id","anno","padj_fdr","r","sig")] %>%
#   setnames(c("id"),c("ens_id")) %>% 
#   merge(dt, by=c("ens_id","anno"))
# 
# 
# dt.cor[,median(r),by="color"]
# 
# my_comparisons <- list(c("red","blue"))
# 
# p <- ggplot(dt.cor, aes(x=color, y=r, fill=color)) +
#   labs(x="", y="Met/RNA correlation") +
#   geom_violin(alpha=0.6, size=0.25) +
#   geom_boxplot(alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
#   geom_jitter(alpha=0.6, size=1.5, stroke=0.1, shape=21, color="black", width=0.1) +
#   scale_fill_manual(values=c("black"="black", "red"="#4E934D", "blue"="#C400AD")) +
#   ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") +
#   theme_classic() +
#   theme(
#     axis.ticks.x = element_blank(),
#     legend.position = "none",
#     axis.text.y = element_text(color="black", size=rel(0.9)),
#     axis.text.x = element_blank()
#   )
# 
# 
# # geom_point(shape=21, stroke=0.15, color="black") +
# #   # geom_point() +
# #   scale_fill_manual(values=c("black"="gray80", "red"="#4E934D", "blue"="#C400AD")) +
# #   scale_size_manual(values=c("black"=0.5, "red"=2, "blue"=2)) +
# #   scale_alpha_manual(values=c("black"=0.7, "red"=0.9, "blue"=0.9)) +
# #   
# pdf(paste0(io$outdir,"/boxplots_overdispersion_vs_MetRNAcor.pdf"), width=4, height=4, useDingbats = F)
# print(p)
# dev.off()
# 