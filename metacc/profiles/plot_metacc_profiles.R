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
args$outdir  <- file.path(io$basedir,"results_new/metacc/profiles/multiome_peaks")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_cell"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

# Options

###########################
## Load precomputed data ##
###########################

metacc.dt <- fread(args$file)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>%
  .[cell%in%unique(metacc.dt$cell)]

##########################################
## Plot TSS profiles one cell at a time ##
##########################################

cells.to.plot <- unique(metacc.dt$cell)# %>% head(n=5)

opts$window_size <- max(metacc.dt$dist)

for (i in cells.to.plot) {
  print(i)
  
  to.plot <- metacc.dt[cell==i] %>%
    # .[,.(rate=mean(rate)),by=c("cell","dist","context")] %>%
    merge(sample_metadata[cell==i,c("cell","met_rate","acc_rate")] %>% setnames(c("cell","global_met_rate","global_acc_rate")))
  
  p <- ggplot(to.plot, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    stat_summary(geom="ribbon", fun.data="mean_se", alpha=1, color="black") +
    geom_hline(yintercept=to.plot$global_met_rate, color="#F37A71", linetype="dashed", alpha=0.75, size=1) +
    geom_hline(yintercept=to.plot$global_acc_rate, color="#00BFC4", linetype="dashed", alpha=0.75, size=1) +
    # stat_summary(geom="line", fun.data="mean_se") +
    # geom_line(size=2) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(0,80)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill="none", color="none", linetype="none") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(1.0), colour="black"),
      axis.text.y = element_text(size=rel(1.1), colour="black")
    )
  # print(p_list[[i]])
  
  pdf(file.path(args$outdir,sprintf("per_cell/%s.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}

###########################################
## Plot TSS profiles one class at a time ##
###########################################

# plates.to.plot <- sample_metadata[pass_metQC==TRUE & pass_accQC==TRUE,plate] %>% unique
classes.to.plot <- unique(sample_metadata$class)

opts$window_size <- max(metacc.dt$dist)

for (i in classes.to.plot) {
  print(i)
  
  # Note that we only use cells that pass quality control for both met and acc
  to.plot <- metacc.dt %>%
    .[cell%in%sample_metadata[pass_metQC==TRUE & pass_accQC==TRUE & class==i,cell]] %>%
    .[,.(rate=mean(rate), N=sum(N)),by=c("cell","dist","context")] %>%
    .[N>=50]
  
  p <- ggplot(to.plot, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    stat_summary(geom="ribbon", fun.data="mean_sd", alpha=1, color="black") +
    geom_hline(yintercept=mean(sample_metadata[class==i,met_rate],na.rm=T), color="#F37A71", linetype="dashed", alpha=0.75, size=1) +
    geom_hline(yintercept=mean(sample_metadata[class==i,acc_rate],na.rm=T), color="#00BFC4", linetype="dashed", alpha=0.75, size=1) +
    # stat_summary(geom="line", fun.data="mean_se") +
    # geom_line(size=2) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(0,80)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill="none", color="none", linetype="none") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(1.0), colour="black"),
      axis.text.y = element_text(size=rel(1.1), colour="black")
    )
  
  pdf(file.path(args$outdir,sprintf("per_class/%s.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}


to.plot <- metacc.dt %>%
  merge(sample_metadata[,c("cell","class")]) %>%
  .[,.(rate=mean(rate), N=sum(N)),by=c("cell","dist","context","class")] %>%
  .[N>=50]

p <- ggplot(to.plot, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
  stat_summary(geom="ribbon", fun.data="mean_sd", alpha=1, color="black") +
  geom_hline(aes(yintercept=V1), color="#F37A71", linetype="dashed", alpha=0.75, size=1, data=sample_metadata[,mean(met_rate),by="class"]) +
  geom_hline(aes(yintercept=V1), color="#00BFC4", linetype="dashed", alpha=0.75, size=1, data=sample_metadata[,mean(acc_rate),by="class"]) +
  facet_wrap(~class, nrow=1) +
  labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
  coord_cartesian(ylim=c(10,80)) +
  xlim(-opts$window_size, opts$window_size) +
  guides(fill="none", color="none", linetype="none") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1.15), colour="black"),
    axis.text.x = element_text(size=rel(1.0), colour="black"),
    axis.text.y = element_text(size=rel(1.1), colour="black")
  )

pdf(file.path(args$outdir,"per_class/metacc_profiles.pdf"), width=7, height=5)
print(p)
dev.off()

###########################################################
## Scatterplots or CG vs GC promoter-to-background ratio ##
###########################################################

# Note that this is only applied for TSS

opts$min_gc_ratio <- 1.5
opts$max_cg_ratio <- 0.5

to.plot <- metacc.dt %>%
  .[,.(ratio=mean(.SD[abs(dist)<=50,rate])/mean(.SD[abs(dist)>2500,rate])), by=c("cell","context")] %>%
  .[,context:=paste0("ratio_",context)] %>%
  dcast(cell~context, value.var="ratio") %>%
  merge(sample_metadata[,c("cell","pass_metQC","pass_accQC","nCG","nGC")], by="cell") %>%
  .[,c("log_nCG","log_nGC"):=list(log(nCG),log(nGC))]

p <- ggscatter(to.plot, x="ratio_CG", y="ratio_GC", fill="log_nCG", shape=21, size=1.5) +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  coord_cartesian(xlim = c(0.01,0.75), ylim=c(0.5,5)) +
  geom_hline(yintercept=opts$min_gc_ratio, linetype="dashed", color="black") +
  geom_vline(xintercept=opts$max_cg_ratio, linetype="dashed", color="black")

pdf(file.path(args$outdir,"CG_ratio_vs_GC_ratio_coloured_by_nCG.pdf"), width=6, height=5)
print(p)
dev.off()

p <- ggscatter(to.plot, x="ratio_CG", y="ratio_GC", fill="log_nGC", shape=21, size=1.5) +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  coord_cartesian(xlim = c(0.01,0.75), ylim=c(0.5,5)) +
  geom_hline(yintercept=opts$min_gc_ratio, linetype="dashed", color="black") +
  geom_vline(xintercept=opts$max_cg_ratio, linetype="dashed", color="black")

pdf(file.path(args$outdir,"CG_ratio_vs_GC_ratio_coloured_by_nGC.pdf"), width=6, height=5)
print(p)
dev.off()