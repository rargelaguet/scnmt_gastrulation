
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
  # "H3K27ac_distal_E7.5_union_intersect12_500",
  # "H3K27ac_distal_E7.5_union_intersect12",
  "prom_2000_2000"
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_Ect_intersect12"
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

# rna.diff <- fread("/Users/ricard/data/scnmt_gastrulation/rna/results/differential/E4.5E5.5Epiblast_vs_E7.5EctodermEndodermMesoderm.txt.gz") %>%
#   .[abs(logFC)>3 & padj_fdr<0.01]

##########
## Link ##
##########

dt <- merge(
  met_overdispersion.dt,
  rna_overdispersion.dt
) %>% merge(gene_metadata,by="ens_id")

table(dt$foo)
table(dt$bar)

dt[,color:="black"]
dt[gamma>0.3 & bio.disp>2,color:="red"]
dt[gamma<0.2 & bio.disp>4,color:="blue"]

p <- ggplot(dt, aes(x=gamma, y=bio.disp, fill=color, size=color)) +
  labs(x="Met. overdispersion", y="RNA overdispersion") +
  geom_point(shape=21, stroke=0.25, color="black") +
  scale_fill_manual(values=c("black"="black", "red"="red", "blue"="blue")) +
  scale_size_manual(values=c("black"=0.5, "red"=1.5, "blue"=1.5)) +
  ggrepel::geom_text_repel(data=dt[color=="red"], aes(x=gamma, y=bio.disp, label=symbol), size=4, color="red") +
  ggrepel::geom_text_repel(data=dt[color=="blue"], aes(x=gamma, y=bio.disp, label=symbol), size=4, color="blue") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(paste0(io$outdir,"/met_vs_rna_overdispersion.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()


###############
## Plot hits ##
###############

io$script <- "/Users/ricard/scnmt_gastrulation/metrna/dynamics_individual_examples/boxplots_argparse.R"

args <- list()
args$met.anno <- "prom_2000_2000"
args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
args$outdir <- "/Users/ricard/data/scnmt_gastrulation/metrna/plot_individual_examples"

# to.plot <- dt[gamma>0.25 & bio.disp<(-2.5)]
to.plot <- dt[gamma>0.3 & bio.disp>2]
# to.plot <- dt[gamma<0.2 & bio.disp>4]

for (i in 1:nrow(to.plot)) {
  # args$gene <- args$met.id <- args$acc.id <- i
  args$met.id <- args$acc.id <- to.plot[i,ens_id]
  args$gene <- to.plot[i,symbol]
  
  cmd <- sprintf("Rscript %s --gene %s --met.id %s --met.anno %s --stage_lineage %s --outdir %s", 
                 io$script, args$gene, args$met.id, args$met.anno, paste(args$stage_lineage,collapse=" "), args$outdir)
  system(cmd)
}


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
  geom_jitter(alpha=0.5, size=1, stroke=0.2, shape=21, color="black", width=0.1) +
  scale_fill_manual(values=c("black"="black", "red"="red", "blue"="blue")) +
  ggpubr::stat_compare_means(comparisons = my_comparisons) +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.text.x = element_blank()
  )
p

pdf(paste0(io$outdir,"/boxplots_overdispersion_vs_MetRNAcor.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()