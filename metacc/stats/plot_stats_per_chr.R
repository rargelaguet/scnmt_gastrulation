here::i_am("metacc/stats/plot_stats_per_chr.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--stats',    type="character",  help='Stats file')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
# args$stats <- file.path(io$basedir,"results/met/stats/met_stats_per_chr.txt.gz")
# args$outdir <- file.path(io$basedir,"results/met/stats/pdf")
# args$context <- "CG"
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# I/O
dir.create(args$outdir, showWarnings=F)

# Optioms
opts$chr <- paste0("chr",c(1:19,"X"))

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) 

################
## Load stats ##
################

stats.dt <- fread(args$stats) %>% 
  .[chr%in%opts$chr] %>% .[,chr:=factor(chr,levels=opts$chr)]

# Merge with cell metadata
if (args$context=="CG") {
  to.plot <- stats.dt %>%
    merge(sample_metadata[pass_metQC==TRUE,c("cell","id_met","sample","stage","celltype")],by="id_met") %>%
    setnames(c("nCG","met_rate"),c("N","rate"))
} else if (args$context=="GC") {
  to.plot <- stats.dt %>%
    merge(sample_metadata[pass_accQC==TRUE,c("cell","id_acc","sample","stage","celltype")],by="id_acc") %>%
    setnames(c("nGC","acc_rate"),c("N","rate"))
}

##################################################
## Boxplots with rate per chromosome and sample ##
##################################################

p <- ggboxplot(to.plot[N>=100], x = "sample", y = "rate", outlier.shape=NA, fill="sample", alpha=0.5) +
  facet_wrap(~chr, scales="fixed") +
  scale_fill_brewer(palette="Dark2") +
  labs(x="", y="Rate") +
  theme(
    axis.text.y = element_text(size=rel(0.80)),
    axis.text.x = element_blank(),
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,sprintf("%s_rate_per_chr.pdf",args$context)), width=11, height=9)
print(p)
dev.off()


