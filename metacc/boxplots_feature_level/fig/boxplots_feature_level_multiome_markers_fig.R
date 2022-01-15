here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

suppressMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",  help='Cell metadata')
p$add_argument('--met_file',  type="character",  help='DNA methylation file with the feature level quantification')
p$add_argument('--acc_file',  type="character",  help='Chr. accessibility file with the feature level quantification')
p$add_argument('--anno',      type="character",  help='Genomic annotation')
p$add_argument('--markers_file',      type="character",  help='Celltype marker peaks file from the 10x Multiome atlas')
p$add_argument('--outdir',    type="character",  help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$anno <- "multiome_peaks"
# args$met_file <- file.path(io$basedir,sprintf("processed/met/feature_level/%s.tsv.gz",args$anno))
# args$acc_file <- file.path(io$basedir,sprintf("processed/acc/feature_level/%s.tsv.gz",args$anno))
# args$markers_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"
# args$outdir  <- file.path(io$basedir,"results_new/metacc/boxplots_feature_level/markers")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"absolute"), showWarnings = F)
dir.create(file.path(args$outdir,"relative"), showWarnings = F)

# Options
opts$min_observations <- 50
opts$min_cells <- 10
opts$min_marker_score <- 0.75

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>% 
  .[!is.na(celltype.mapped) & (pass_metQC==TRUE | pass_accQC==TRUE)]

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

###############
## Load data ##
###############

met.dt <- fread(args$met_file) %>% 
  setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_met%in%opts$met_cells]

acc.dt <- fread(args$acc_file) %>% 
  setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_acc%in%opts$acc_cells]

# Add common cell identifier
met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
  .[,id_met:=NULL] %>% .[,context:="CG"]
acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
  .[,id_acc:=NULL] %>% .[,context:="GC"]

# Merge
metacc.dt <- rbind(met.dt,acc.dt)

#######################
## Load marker peaks ##
#######################

marker_peaks.dt <- fread(args$markers_file)

scnmt.peaks <- unique(metacc.dt$id)
multiome.peaks <- unique(marker_peaks.dt$idx)

length(scnmt.peaks)
length(multiome.peaks)

# length(intersect(unique(marker_peaks.dt$idx),unique(metacc.dt$id)))

############################
## Calculate global rates ##
############################

global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>% 
  setnames(c("cell","CG","GC")) %>%
  melt(id.vars="cell", variable.name="context", value.name="global_rate")

# Normalise rates by genome-wide levels
# to.plot <- metacc.dt %>%
#   .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
#   .[Ntotal>=opts$min_observations] %>%
#   .[,rate:=100*(Nmet/Ntotal)] %>%
#   merge(sample_metadata[,c("cell","sample","celltype.mapped")], by="cell") %>%
#   merge(global_rates.dt, by=c("cell","context")) %>%
#   .[,rate_norm:=rate/global_rate]

##########################
## Plot absolute levels ##
##########################

celltypes.to.plot <- unique(marker_peaks.dt$celltype)# %>% head(n=3)

celltypes.to.subset <- sample_metadata[,.N,by=c("celltype.mapped","class")] %>% .[N>=opts$min_cells] %>% .[,.N,by="celltype.mapped"] %>% .[N>1] %>% .$celltype.mapped

for (i in celltypes.to.plot) {
  
  # sum(scnmt.peaks %in% marker_peaks.dt[celltype==i & score>=opts$min_marker_score,idx])
  
  tmp <- sample_metadata[celltype.mapped%in%celltypes.to.subset,c("cell","class","sample","celltype.mapped")]
  
  to.plot <- metacc.dt %>%
    .[id%in%marker_peaks.dt[celltype==i & score>=opts$min_marker_score,idx]]  %>%
    .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
    .[Ntotal>=opts$min_observations] %>%
    .[,rate:=100*(Nmet/Ntotal)] %>%
    merge(tmp, by="cell")
  
  p.met <- ggplot(to.plot[context=="CG"], aes(x = class, y = rate, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 1.25, alpha = 0.75, shape=21) +
    facet_wrap(~celltype.mapped) +
    labs(x="", y="Methylation levels (%)") +
    geom_hline(yintercept=mean(global_rates.dt[cell%in%tmp$cell & context=="CG",global_rate]), linetype="dashed") +
    # scale_fill_manual(values=opts$context.colors) +
    coord_cartesian(ylim=c(0,100)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p.acc <- ggplot(to.plot[context=="GC"], aes(x = class, y = rate, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 1.25, alpha = 0.75, shape=21) +
    facet_wrap(~celltype.mapped) +
    geom_hline(yintercept=mean(global_rates.dt[cell%in%tmp$cell & context=="GC",global_rate]), linetype="dashed") +
    labs(x="", y="Chr. accessibility levels (%)") +
    # scale_fill_manual(values=opts$context.colors) +
    coord_cartesian(ylim=c(8,60)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(file.path(args$outdir,sprintf("absolute/boxplots_metacc_%s_markers.pdf",i)), width=10, height=6)
  print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 1))
  dev.off()
}


##########################
## Plot relative levels ##
##########################

for (i in celltypes.to.plot) {
  
  to.plot <- metacc.dt %>%
    .[id%in%marker_peaks.dt[celltype==i & score>=opts$min_marker_score,idx]]  %>%
    .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
    .[Ntotal>=opts$min_observations] %>%
    .[,rate:=100*(Nmet/Ntotal)] %>%
    merge(sample_metadata[celltype.mapped%in%celltypes.to.subset,c("cell","class","sample","celltype.mapped")], by="cell") %>%
    merge(global_rates.dt, by=c("cell","context")) %>%
    .[,rate_norm:=rate/global_rate]
  
  p.met <- ggplot(to.plot[context=="CG"], aes(x = class, y = rate_norm, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 1.25, alpha = 0.75, shape=21) +
    facet_wrap(~celltype.mapped) +
    geom_hline(yintercept=1, linetype="dashed") +
    labs(x="", y="Methylation levels (%)") +
    # scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(0,100)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p.acc <- ggplot(to.plot[context=="GC"], aes(x = class, y = rate_norm, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 1.25, alpha = 0.75, shape=21) +
    facet_wrap(~celltype.mapped) +
    labs(x="", y="Chr. accessibility levels (%)") +
    geom_hline(yintercept=1, linetype="dashed") +
    # scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(8,60)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(file.path(args$outdir,sprintf("relative/boxplots_metacc_%s_markers.pdf",i)), width=10, height=6)
  print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 1))
  dev.off()
}
