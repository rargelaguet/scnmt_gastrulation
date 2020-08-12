#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/settings.R")

# Define I/O
# io$mm10.genome <- "/Users/ricard/data/mm10_sequence/mm10.genome"
io$outdir <- paste0(io$basedir,"/met/results/stats/pdf")

# Define options
# opts$chr <- c(paste0("chr",1:19),"X")
# opts$chr <- c(1:19,"X")

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[!is.na(id_met) & !is.na(lineage10x_2)] %>%
  droplevels

sample_metadata[,length(unique(embryo)),by=c("stage_lineage","sex")]

############################
## Load precomputed stats ##
############################

stats <- fread(io$met.stats_per_chr) %>%
  .[chr%in%opts$chr] %>%
  .[,chr:=factor(chr,levels=opts$chr)] %>%
  # .[,mean_coverage:=mean(coverage),by="id_met"] %>%
  # .[,relative_coverage:=coverage/mean_coverage,by="id_met"] %>%
  merge(sample_metadata, by="id_met") %>%
  .[,sex:=factor(sex, levels=c("Female","Male"))]


# Sanity check
# all(sort(unique(stats$id_met)) == sort(sample_metadata$id_met))

# Load chromosome length
mm10.genome <- fread(io$mm10.genome) %>%
  setnames(c("chr","chr_length")) %>%
  .[chr%in%opts$chr] %>% .[,chr:=factor(chr,levels=opts$chr)]

####################################
## Plot relative coverage per chr ##
####################################

to.plot <- stats %>%
  .[,.(coverage=as.double(sum(coverage))), by=c("embryo","chr")] %>%
  merge(mm10.genome, by="chr") %>% 
  .[,coverage:=coverage/as.double(chr_length), by="embryo"]# %>%
  # .[,norm_coverage:=coverage/mean(coverage),by="embryo"]

ggscatter(to.plot, x="embryo", y="coverage") +
  # geom_hline(yintercept=1, linetype="dashed") +
  facet_wrap(~chr, scales="fixed") +
  labs(x="", y="normalised DNAm coverage") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###################################
## Plot methylation rate per chr ##
###################################

to.plot <- stats %>%
  .[sex%in%c("Female","Male")] %>% droplevels %>% 
  .[,.(rate=mean(mean)), by=c("embryo","chr","sex","stage_lineage")]


for (i in unique(to.plot$stage_lineage)) {
  p <- ggboxplot(to.plot[stage_lineage==i], x="sex", y="rate", fill="sex") +
    # geom_hline(yintercept=1, linetype="dashed") +
    facet_wrap(~chr, scales="fixed") +
    labs(x="", y="Global DNA methylation (%)") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(sprintf("%s/methylation_per_chr_%s.pdf",io$outdir,i), width=8, height=10)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
