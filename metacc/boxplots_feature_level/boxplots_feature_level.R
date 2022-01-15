here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",  help='Cell metadata')
p$add_argument('--met_file',  type="character",  help='DNA methylation file with the feature level quantification')
p$add_argument('--acc_file',  type="character",  help='Chr. accessibility file with the feature level quantification')
p$add_argument('--anno',      type="character",  help='Genomic annotation')
p$add_argument('--outdir',    type="character",  help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Define settings ##
###################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$anno <- "prom_200_200"
# args$met_file <- file.path(io$basedir,sprintf("processed/met/feature_level/%s.tsv.gz",args$anno))
# args$acc_file <- file.path(io$basedir,sprintf("processed/acc/feature_level/%s.tsv.gz",args$anno))
# args$outdir  <- file.path(io$basedir,"results/metacc/boxplots_feature_level")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

# Options
opts$min_observations <- 250

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

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

###########
## Merge ##
###########

# Add common cell identifier
met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
  .[,id_met:=NULL] %>% .[,context:="CG"]
acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
  .[,id_acc:=NULL] %>% .[,context:="GC"]

# Merge
metacc.dt <- rbind(met.dt,acc.dt)

###############################
## Prepare data for plotting ##
###############################

global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>% 
  setnames(c("cell","CG","GC")) %>%
  melt(id.vars="cell", variable.name="context", value.name="global_rate")

# Normalise rates by genome-wide levels
to.plot <- metacc.dt %>%
  .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
  .[Ntotal>=opts$min_observations] %>%
  .[,rate:=100*(Nmet/Ntotal)] %>%
  merge(sample_metadata[,c("cell","celltype")], by="cell") %>%
  merge(global_rates.dt, by=c("cell","context")) %>%
  .[,rate_norm:=rate/global_rate]

##########################
## Plot absolute levels ##
##########################

p <- ggplot(to.plot, aes(x = celltype, y = rate, fill=context)) +
  # facet_wrap(~context) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.5, alpha = 0.4, shape=21) +
  facet_wrap(~context, scales="free_y") +
  labs(x="", y="Met/Acc levels (%)") +
  scale_color_manual(values=opts$context.colors) +
  scale_fill_manual(values=opts$context.colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black")
  )

pdf(file.path(args$outdir,sprintf("boxplots_metacc_%s.pdf",args$anno)), width=8, height=6)
print(p)
dev.off()


##########################
## Plot relative levels ##
##########################

p <- ggplot(to.plot, aes(x = celltype, y = rate_norm, fill=context)) +
  # facet_wrap(~context) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.5, alpha = 0.4, shape=21) +
  facet_wrap(~context, scales="free_y") +
  geom_hline(yintercept = 1, linetype="dashed") +
  labs(x="", y="Met/Acc levels (relative to background)") +
  scale_color_manual(values=opts$context.colors) +
  scale_fill_manual(values=opts$context.colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black")
  )

pdf(file.path(args$outdir,sprintf("boxplots_metacc_%s_relative.pdf",args$anno)), width=8, height=6)
print(p)
dev.off()
