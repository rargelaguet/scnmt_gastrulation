here::here("metrna/coupling/local_coupling/plot_metrna_coupling.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--file',    type="character",    help='Precomputed file')
p$add_argument('--outdir',  type="character",    help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
# args <- list()
# args$metadata <- io$metadata# file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$file  <- file.path(io$basedir,"results/metrna/coupling/precomputed_metrna_coupling.txt.gz")
# args$outdir  <- file.path(io$basedir,"results/metrna/coupling")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_cell"), showWarnings = F)

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_metQC==TRUE & pass_rnaQC==TRUE]

###########################
## Load precomputed data ##
###########################

metrna_coupling.dt <- fread(args$file)

###################
## Plot per cell ##
###################

# cells.to.plot <- unique(metrna_coupling.dt$cell)# %>% head(n=5)
# 
# for (i in cells.to.plot) {
#   
#   to.plot <- metrna_coupling.dt[cell==i]
# 
#   p <- ggplot(to.plot, aes(x=window_center, y=r)) +
#     stat_summary(fun.data=mean_sd, geom="smooth", alpha=0.2, size=1.0, color="black", fill="black") +
#     geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
#     geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5) +
#     labs(x="Genomic distance from TSS (bp)", y="Met/RNA correlation") +
#     # coord_cartesian(ylim=c(-0.75,0.5)) +
#     theme_classic() +
#     theme(
#       axis.text = element_text(color="black", size=rel(0.8))
#     )
#   
#   pdf(file.path(args$outdir,sprintf("per_cell/%s.pdf",i)), width=6, height=5)
#   print(p)
#   dev.off()
# }

#####################
## Plot per sample ##
#####################

to.plot <- metrna_coupling.dt %>% 
  merge(sample_metadata[,c("cell","sample")]) %>%
  .[,.(r=mean(r,na.rm=T)), by=c("window_center","sample")]

p <- ggplot(to.plot, aes(x=window_center, y=r, color=sample)) +
  geom_line(size=1.25) +
  geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5) +
  scale_color_brewer(palette="Dark2") +
  labs(x="Genomic distance from TSS (bp)", y="Met/RNA correlation") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(file.path(args$outdir,"metrna_coupling_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()
