here::here("metacc/profiles/calculate_metacc_profiles.R")

suppressMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",  help='Cell metadata')
# p$add_argument('--anno',    type="character",    help='Genomic annotation')
p$add_argument('--file',    type="character",    help='Precomputed file')
p$add_argument('--outdir',  type="character",    help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
args$file  <- file.path(io$basedir,"results_new/metacc/profiles/multiome_peaks/precomputed_metacc_multiome_peaks.txt.gz")
# args$markers_file <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"
args$markers_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"
args$outdir  <- file.path(io$basedir,"results_new/metacc/profiles/multiome_peaks")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

# Options
opts$celltypes = c(
  "Surface_ectoderm",
  # "Gut",
  "Pharyngeal_mesoderm",
  "Endothelium",
  "Haematoendothelial_progenitors",
  "Blood_progenitors",
  "Erythroid"
)

opts$rename.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "early_Erythroid" = "Erythroid",
  "late_Erythroid" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm"
)


###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[,celltype:=stringr::str_replace_all(celltype.mapped,opts$rename.celltypes)] %>%
  .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>% .[,class:=factor(class,levels=c("WT","TET-TKO"))] %>%
  .[(pass_metQC==TRUE | pass_accQC==TRUE) & celltype%in%opts$celltypes]

table(sample_metadata$celltype,sample_metadata$class)

#######################
## Load marker peaks ##
#######################

opts$min_marker_score <- 0.75

marker_peaks.dt <- fread(args$markers_file) %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[celltype%in%opts$celltypes] %>%
  .[,.(score=mean(score)), by=c("celltype","idx")] %>%
  .[score>=opts$min_marker_score]

table(marker_peaks.dt$celltype)

###########################
## Load precomputed data ##
###########################

# metacc.dt <- fread(args$file) %>%
#   .[cell%in%sample_metadata$cell & id%in%unique(marker_peaks.dt$idx)]
# fwrite(metacc.dt, file.path(args$outdir,"precomputed_metacc_multiome_peaks_filt.txt.gz"))
metacc.dt <- fread(file.path(args$outdir,"precomputed_metacc_multiome_peaks_filt.txt.gz")) %>%
  .[cell%in%sample_metadata$cell & id%in%unique(marker_peaks.dt$idx)]

###########################################
## Plot TSS profiles one class at a time ##
###########################################

# classes.to.plot <- unique(sample_metadata$class)

opts$window_size <- max(metacc.dt$dist)

celltypes.to.plot <- opts$celltypes

# i <- celltypes.to.plot[1]
for (i in celltypes.to.plot) {
  print(i)
  
  # Note that we only use cells that pass quality control for both met and acc
  
  # Methylation
  to.plot.met <- metacc.dt %>%
    .[cell%in%sample_metadata[pass_metQC==TRUE & celltype==i,cell]] %>%
    merge(marker_peaks.dt[,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=T) %>%
    # .[id%in%marker_peaks.dt[celltype==i,idx]] %>%
    .[,.(rate=mean(rate), N=sum(N)),by=c("dist","context","cell","celltype_marker")]# %>%
    # .[N>=25]
  
  to.plot.acc <- metacc.dt %>%
    .[cell%in%sample_metadata[pass_accQC==TRUE & celltype==i,cell]] %>%
    merge(marker_peaks.dt[,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=T) %>%
    # .[id%in%marker_peaks.dt[celltype==i,idx]] %>%
    .[,.(rate=mean(rate), N=sum(N)),by=c("dist","context","cell","celltype_marker")]# %>%
    # .[N>=50]
  
  to.plot <- rbind(to.plot.met,to.plot.acc) %>% merge(sample_metadata[,c("cell","class")])
  to.plot[,tmp:=sprintf("(%s) %s",class,celltype_marker)]
  
  to.plot.lines <- sample_metadata[celltype==i,.(CG=mean(met_rate,na.rm=T), GC=mean(acc_rate,na.rm=T)), by=c("class")] %>%
    melt(id.vars=c("class"), variable.name="context", value.name="rate")
  
  p <- ggplot(to.plot, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    stat_summary(geom="ribbon", fun.data="mean_sd", alpha=1, color="black") +
    # stat_summary(geom="line", fun.data="mean_se") +
    # geom_line(size=2) +
    facet_wrap(~celltype_marker+class, scales="fixed", ncol=2) +
    geom_hline(aes(yintercept=rate, color=context), linetype="dashed", alpha=0.75, size=0.75, data=to.plot.lines) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(10,100)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill="none", color="none", linetype="none") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(1.0), colour="black"),
      axis.text.y = element_text(size=rel(1.1), colour="black")
    )
  
  pdf(file.path(args$outdir,sprintf("per_class/%s.pdf",i)), width=4, height=8)
  print(p)
  dev.off()
}

