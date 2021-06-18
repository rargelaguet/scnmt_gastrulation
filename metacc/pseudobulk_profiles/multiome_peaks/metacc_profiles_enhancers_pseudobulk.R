
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
}

opts$test <- TRUE

if (opts$test) opts$stage_lineage <- opts$stage_lineage %>% head(n=2)

##############################
## Load genomic annotations ##
##############################

# Define genomic contexts
opts$annos <- c(
  "prom_2000_2000" = "Promoters"
)

# Define window positions and characteristics
opts$positions <- c(
  "prom_2000_2000"="center"
)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_annotations.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_annotations.R")
}
anno_df.met <- anno_df
anno_df.acc <- anno_df

#################################################################
## Load genome-wide global methylation and accessibility stats ##
#################################################################

####################################
## Load pseudobulk scNMT-seq data ##
####################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_pseudobulk_data.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_pseudobulk_data.R")
}

# Merge
data <- rbind(
  met.dt[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")],
  acc.dt[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")]
)

# Rename genomic annotations
# data[,anno:=stringr::str_replace_all(anno,opts$annos)]

# Load precomputed data
saveRDS(data, paste0(io$outdir,"/pseudobulk_data.rds"))
data <- readRDS(paste0(io$outdir,"/pseudobulk_data.rds"))

###########################
## Split by marker peaks ##
###########################

# Load celltype marker peaks
multiome_marker_peaks.dt <- fread(io$multiome.atac.marker_peaks)

##########
## Plot ##
##########

# data[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]
# stats[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]

for (i in opts$stage_lineage) {
  print(i)
  
  tmp <- data[stage_lineage==i]
  
  p <- ggplot(tmp, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    facet_wrap(~anno, nrow=1, scales="fixed") +
    stat_summary(geom="ribbon", fun.data="mean_se", alpha=1) +
    stat_summary(geom="line", fun.data="mean_se") +
    geom_hline(yintercept=stats[context=="CG" & stage_lineage==i,median(mean,na.rm=T)], color="#F37A71", linetype="dashed", alpha=0.75, size=0.75) +
    geom_hline(yintercept=stats[context=="GC" & stage_lineage==i,median(mean,na.rm=T)], color="#00BFC4", linetype="dashed", alpha=0.75, size=0.75) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(0,100)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill=FALSE, color=FALSE, linetype=FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(0.8), colour="black"),
      axis.text.y = element_text(size=rel(1.2), colour="black")
    )

  pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8.5, height=5)
  print(p)
  dev.off()
}

