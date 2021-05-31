###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/acc/results/motifs/test")


# Options

opts$annos <- c("multiome_peaks") # H3K27ac_distal_E7.5_union_intersect12
opts$motif.annotation <- "jaspar2020"

opts$stage_lineage <- c(
  # "E4.5_Epiblast",
  # "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm"
)

opts$to.merge <- c(
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Intermediate_mesoderm" = "Mature_mesoderm",
  "Pharyngeal_mesoderm" = "Mature_mesoderm",
  "Paraxial_mesoderm" = "Mature_mesoderm",
  "Somitic_mesoderm" = "Mature_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

opts$remove.small.lineages <- TRUE

opts$min.CpGs <- 1


########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_accQC==TRUE & stage_lineage%in%opts$stage_lineage] %>%
  .[,lineage10x:=stringr::str_replace_all(lineage10x,opts$to.merge)] %>%
  .[,c("sample","id_acc","stage","lineage10x","lineage10x_2","stage_lineage")]
  # .[,lineage10x_2:=stringr::str_replace_all(lineage10x_2,"_"," ")] %>% 

if (opts$remove.small.lineages) {
  opts$min.cells <- 5
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("lineage10x")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

table(sample_metadata$lineage10x)


#################
## Load motifs ##
#################

motifOverlap.se <- readRDS(sprintf("%s/%s_%s.rds",io$motifs.dir,opts$anno,opts$motif.annotation))
opts$motifs <- colnames(motifOverlap.se$motifMatches)# %>% head(n=2)


###############
## Load data ##
###############

acc_dt <- opts$annos %>% map(function(i) {
  opts$motifs %>% map(function(j) {
    file <- sprintf("%s/%s_%s_%s.tsv.gz",io$acc_data_motifs,i,opts$motif.annotation,j)
    if (file.exists(file)) {
      fread(file, header=F, colClasses = c("factor","character","factor","integer","integer","integer")) %>%
        .[V1%in%sample_metadata$id_acc] %>%
        setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>%
        .[Ntotal>=opts$min.CpGs] %>%
        .[,.(rate=100*(sum(Nmet)/sum(Ntotal)), Nmet=sum(Nmet), N=sum(Ntotal)),by=c("id_acc","anno")]
    }
  }) %>% rbindlist
}) %>% rbindlist


fwrite(acc_dt, paste0(io$outdir,"/motif_acc.txt.gz"), sep="\t")


##################
## Prepare data ##
##################

opts$min.observations <- 25

to.plot <- acc_dt %>%
  .[N>=opts$min.observations] %>%
  merge(sample_metadata, by="id_acc") %>%
  .[,stage_lineage:=factor(stage_lineage, levels=names(opts$stagelineage.colors))]  %>%
  .[,motif:=gsub(sprintf("%s_%s_",opts$annos,opts$motif.annotation),"",anno)]


# Regress out global rates
foo <- fread(io$acc.stats) %>% .[,c("id_acc","mean")]
to.plot <- to.plot %>% merge(foo, by="id_acc") %>%
  .[,rate:=lm(formula=rate~mean)[["residuals"]], by=c("motif","anno")]


##########
## Plot ##
##########

for (i in unique(to.plot$motif)) {
  p <- ggplot(to.plot[motif==i], aes(x=stage_lineage, y=rate)) +
    geom_violin(aes(fill=stage_lineage), outlier.shape=NA, coef=1) +
    geom_boxplot(aes(fill=stage_lineage)) +
    scale_fill_manual(values=opts$stagelineage.colors) +
    labs(x="", y="Accessibility (%)") +
    # facet_wrap(~motif, scales="fixed") +
    # coord_cartesian(ylim=c(8,93)) +
    guides(x = guide_axis(angle = 90)) +
    stat_summary(fun.data = give.n, geom = "text", size=2) +
    theme_bw() +
    # guides(color=F, fill=F)# +
    theme(
      # plot.title = element_text(size=rel(1.2), color="black", hjust=0.5),
      # axis.text.x = element_text(size=rel(1.4), color="black", angle=50, vjust=1, hjust=1),
      axis.text = element_text(size=rel(0.5), color="black"),
      axis.title.y = element_text(size=rel(0.7), color="black"),
      # strip.background = element_rect(fill="#F37A71"),
      strip.text = element_text(size=rel(1), color="black"),
      legend.position="none",
      legend.title = element_blank()
    )
  pdf(sprintf("%s/%s_boxplots_acc_lineage2.pdf",io$outdir,i), width=6, height=3)
  print(p)
  dev.off()
}


for (i in unique(to.plot$motif)) {
  p <- ggplot(to.plot[motif==i], aes(x=lineage10x, y=rate)) +
    geom_violin() +
    geom_boxplot(outlier.shape=NA, coef=1) +
    labs(x="", y="Accessibility (%)") +
    geom_hline(yintercept=0, linetype="dashed") +
    # facet_wrap(~motif, scales="fixed") +
    # coord_cartesian(ylim=c(8,93)) +
    guides(x = guide_axis(angle = 90)) +
    stat_summary(fun.data = give.n, geom = "text", size=2) +
    theme_bw() +
    theme(
      axis.text = element_text(size=rel(0.5), color="black"),
      axis.title.y = element_text(size=rel(0.7), color="black"),
      strip.text = element_text(size=rel(1), color="black"),
      legend.position="none",
      legend.title = element_blank()
    )
  pdf(sprintf("%s/%s_boxplots_acc_lineage1.pdf",io$outdir,i), width=8, height=3)
  print(p)
  dev.off()
}
