
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
met_overdispersion.dt <- fread(io$met.overdispersion) %>%
  .[anno%in%opts$anno] %>%
  setnames("feature_name","ens_id")

# Load RNA overdispersion estimates
rna_overdispersion.dt <- fread(io$rna.overdispersion)

# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>% .[,c("ens_id","symbol")]

##########
## Link ##
##########

dt <- merge(
  met_overdispersion.dt,
  rna_overdispersion.dt
) %>% merge(gene_metadata,by="ens_id")

# Define groups
# dt[,color:="black"]
# dt[gamma>0.3 & bio.disp>2,color:="red"]
# dt[gamma<0.2 & bio.disp>4,color:="blue"]

dt[,color:="black"]
dt[epsilon>0.75 & bio.disp>3,color:="red"]
dt[epsilon<0.25& bio.disp>4,color:="blue"]

##########
## Plot ##
##########

# genes.to.highlight <- c(
#   "Pim2","Slc7a3","Cldn6","Lefty2","Trh","Aplnr","Zic3","Dll3","Krt18","Krt8","Id3","Mesp1","Mixl1",
#   "Dppa5a","Dppa4","Peg3","Cbx7","Slc38a4","Tet1","Tex19.1","Dppa2","Spp1","Mest"
# )

p <- ggplot(dt, aes(x=epsilon, y=bio.disp, alpha=color, fill=color, size=color)) +
  labs(x="Met. overdispersion", y="RNA overdispersion") +
  geom_point(shape=21, stroke=0.15, color="black") +
  # geom_point() +
  scale_fill_manual(values=c("black"="gray80", "red"="#4E934D", "blue"="#C400AD")) +
  scale_size_manual(values=c("black"=0.5, "red"=2, "blue"=2)) +
  scale_alpha_manual(values=c("black"=0.7, "red"=0.9, "blue"=0.9)) +
  ggrepel::geom_text_repel(data=dt[color=="red"], aes(x=epsilon, y=bio.disp, label=symbol), size=5, color="#4E934D") +
  ggrepel::geom_text_repel(data=dt[color=="blue"], aes(x=epsilon, y=bio.disp, label=symbol), size=5, color="#C400AD") +
  coord_cartesian(xlim=c(-0.9,1.6)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/met_vs_rna_overdispersion_scatterplot.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()


###############
## Plot hits ##
###############

# io$script <- "/Users/ricard/scnmt_gastrulation/metrna/dynamics_individual_examples/boxplots_argparse.R"
# 
# args <- list()
# args$met.anno <- "prom_2000_2000"
# # args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
# args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Mesoderm")
# args$outdir <- "/Users/ricard/data/scnmt_gastrulation/metrna/plot_individual_examples"
# 
# # to.plot <- dt[color=="red"]
# # to.plot <- dt[color=="blue"]
# # to.plot <- dt[gamma>0.3 & bio.disp<(-2)]
# to.plot <- dt[epsilon<(-1)]
# 
# for (i in 1:nrow(to.plot)) {
#   # args$gene <- args$met.id  <- i
#   args$met.id <- to.plot[i,ens_id]
#   args$gene <- to.plot[i,symbol]
#   
#   cmd <- sprintf("Rscript %s --gene %s --met.id %s --met.anno %s --stage_lineage %s --outdir %s", 
#                  io$script, args$gene, args$met.id, args$met.anno, paste(args$stage_lineage,collapse=" "), args$outdir)
#   system(cmd)
# }


##############################
## Plot Met/RNA correlation ##
##############################

dt.cor <- fread("/Users/ricard/data/gastrulation/metrna/cor/metrna_cor_promoters_v2.txt.gz") %>%
  .[,c("id","anno","padj_fdr","r","sig")] %>%
  setnames(c("id"),c("ens_id")) %>% 
  merge(dt, by=c("ens_id","anno"))


dt.cor[,median(r),by="color"]

my_comparisons <- list(c("red","blue"))

p <- ggplot(dt.cor, aes(x=color, y=r, fill=color)) +
  labs(x="", y="Met/RNA correlation") +
  geom_violin(alpha=0.6, size=0.25) +
  geom_boxplot(alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
  geom_jitter(alpha=0.6, size=1.5, stroke=0.1, shape=21, color="black", width=0.1) +
  scale_fill_manual(values=c("black"="black", "red"="#4E934D", "blue"="#C400AD")) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.text.x = element_blank()
  )


# geom_point(shape=21, stroke=0.15, color="black") +
#   # geom_point() +
#   scale_fill_manual(values=c("black"="gray80", "red"="#4E934D", "blue"="#C400AD")) +
#   scale_size_manual(values=c("black"=0.5, "red"=2, "blue"=2)) +
#   scale_alpha_manual(values=c("black"=0.7, "red"=0.9, "blue"=0.9)) +
#   
pdf(paste0(io$outdir,"/boxplots_overdispersion_vs_MetRNAcor.pdf"), width=4, height=4, useDingbats = F)
print(p)
dev.off()

##########
## TEST ##
##########

foo <- fread("/Users/ricard/data/betabinomial/argelaguet2019/met/variability/betabinomial/bayes/txt/allcells_prom_2000_2000_vb.txt.gz") %>%
  .[,c("gauss_mean","gauss_var","Feature","anno")] %>% setnames("Feature","ens_id") %>%
  merge(rna_overdispersion.dt, by="ens_id") %>% merge(gene_metadata,by="ens_id")

foo[,color:="black"]
foo[ens_id%in%dt[color=="red",ens_id],color:="red"]
foo[ens_id%in%dt[color=="blue",ens_id],color:="blue"]

p1 <- ggplot(foo, aes(x=gauss_var, y=var, fill=color, size=color)) +
  labs(x="Methylation variance (gaussian)", y="RNA variance (gaussian)") +
  geom_point(shape=21, stroke=0.25, color="black") +
  # scale_fill_manual(values=c("black"="black", "red"="red", "blue"="blue")) +
  scale_fill_manual(values=c("black"="black", "red"="#4E934D", "blue"="#C400AD")) +
  scale_size_manual(values=c("black"=0.5, "red"=1.5, "blue"=1.5)) +
  ggrepel::geom_text_repel(data=foo[color=="red" & symbol%in%genes.to.highlight], aes(x=gauss_var, y=var, label=symbol), size=4, color="#4E934D") +
  ggrepel::geom_text_repel(data=foo[color=="blue" & symbol%in%genes.to.highlight], aes(x=gauss_var, y=var, label=symbol), size=4, color="#C400AD") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/met_vs_rna_gaussian_variance.pdf"), width=6, height=4, useDingbats = F)
print(p1)
dev.off()

foo[,binomial_var:=10*gauss_mean*(1-gauss_mean)]

p2 <- ggplot(foo, aes(x=binomial_var, y=var, fill=color, size=color)) +
  labs(x="Methylation variance (binomial, N=10)", y="RNA variance (gaussian)") +
  geom_point(shape=21, stroke=0.25, color="black") +
  # scale_fill_manual(values=c("black"="black", "red"="red", "blue"="blue")) +
  scale_fill_manual(values=c("black"="black", "red"="#4E934D", "blue"="#C400AD")) +
  scale_size_manual(values=c("black"=0.5, "red"=1.5, "blue"=1.5)) +
  ggrepel::geom_text_repel(data=foo[color=="red" & symbol%in%genes.to.highlight], aes(x=binomial_var, y=var, label=symbol), size=4, color="#4E934D") +
  ggrepel::geom_text_repel(data=foo[color=="blue" & symbol%in%genes.to.highlight], aes(x=binomial_var, y=var, label=symbol), size=4, color="#C400AD") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/met_vs_rna_binomial_variance.pdf"), width=6, height=4, useDingbats = F)
print(p2)
dev.off()

p <- cowplot::plot_grid(plotlist = list(p1,p2), nrow=2)

pdf(paste0(io$outdir,"/met_vs_rna_variance.pdf"), width=10, height=8, useDingbats = F)
print(p)
dev.off()
