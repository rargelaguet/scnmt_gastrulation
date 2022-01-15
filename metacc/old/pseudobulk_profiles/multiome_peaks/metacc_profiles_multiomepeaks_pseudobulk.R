
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
}

opts$test <- FALSE

if (opts$test) opts$stage_lineage <- opts$stage_lineage %>% head(n=2)

##############################
## Load genomic annotations ##
##############################

# Define genomic contexts
opts$annos <- c(
  "prom_2000_2000" = "Promoters",
  "multiome_peaks" = "Multiome peaks"
)

# Define window positions and characteristics
opts$positions <- c(
  "prom_2000_2000"="center",
  "multiome_peaks"="center"
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
  met.dt[,c("stage_lineage","id","anno","dist","rate","context")],
  acc.dt[,c("stage_lineage","id","anno","dist","rate","context")]
)

# Rename genomic annotations
# data[,anno:=stringr::str_replace_all(anno,opts$annos)]

# Load precomputed data
saveRDS(data, paste0(io$outdir,"/pseudobulk_data.rds"))
data <- readRDS(paste0(io$outdir,"/pseudobulk_data.rds"))

data[,rate:=as.integer(rate/100)]

#################################################################
## Load genome-wide global methylation and accessibility stats ##
#################################################################

met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>%
  merge(fread(io$metadata) %>% .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% .[,c("stage_lineage","id_met")], by="id_met") %>% 
  .[!grepl("NA",stage_lineage)] %>%
  .[,.(mean=mean(mean)),by="stage_lineage"] %>%
  .[,context:="CG"]

acc.stats <- fread(io$acc.stats) %>% .[,c("id_acc","mean")] %>%
  merge(fread(io$metadata) %>% .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% .[,c("stage_lineage","id_acc")], by="id_acc") %>% 
  .[!grepl("NA",stage_lineage)] %>%
  .[,.(mean=mean(mean)),by="stage_lineage"] %>%
  .[,context:="GC"]

stats <- rbind(
  met.stats[,c("stage_lineage","mean","context")],
  acc.stats[,c("stage_lineage","mean","context")]
)

###########################
## Split by marker peaks ##
###########################

# Load celltype marker peaks
multiome_marker_peaks.dt <- fread(io$multiome.atac.marker_peaks) %>%
  .[,class:=as.factor(sprintf("%s markers",celltype))] %>%
  setnames("idx","id")

# Rename peaks
tmp <- fread(paste0(io$features.dir,"/multiome_peaks.bed.gz"), select=c(1,2,3,5)) %>%
  setnames(c("chr","start","end","id")) %>%
  .[,correct_id:=sprintf("chr%s:%s-%s",chr,start,end)]
data <- data %>% merge(tmp[,c("id","correct_id")],by="id") %>% .[,id:=NULL] %>% setnames("correct_id","id")

##########
## Plot ##
##########

celltypes.to.plot <- c(
  "Epiblast",
  "Primitive_Streak",
  # "Caudal_epiblast",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Somitic_mesoderm",
  "Cardiomyocytes",
  "Haematoendothelial_progenitors",
  "Endothelium",
  # "Blood_progenitors_1",
  "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  # "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
  # "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$stage_lineage <- c(
  # "E3.5_ICM",
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  # "E6.5_Primitive_Streak",
  "E7.5_Ectoderm",
  "E7.5_Mesoderm",
  "E7.5_Endoderm"
)

data2 <- data %>% 
  .[stage_lineage%in%opts$stage_lineage] %>%
  merge(multiome_marker_peaks.dt[celltype%in%celltypes.to.plot,c("class","id")], allow.cartesian=TRUE) %>%
  .[,.(rate=mean(rate)), by=c("dist","context","anno","stage_lineage","class")] %>%
  .[,rate_scaled:=minmax.normalisation(rate),by=c("context")]

for (i in opts$stage_lineage) {
  
  to.plot <- data2[stage_lineage==i]
  
  p <- ggplot(to.plot, aes(x=dist, y=rate_scaled, group=context, fill=context, color=context)) +
    # geom_line(size=0.5, color="black") +
    geom_ribbon(aes(ymin = rate_scaled - 0.01, ymax = rate_scaled + 0.01), color="black") +
    facet_wrap(~class, ncol=4, scales="fixed") +
    # stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.1,  fun.args=list(mult = 1)) +
    # stat_summary(geom="ribbon", fun.data="mean_se", alpha=0.2) +
    # stat_summary(geom="line", fun.data="mean") +
    # geom_hline(yintercept=stats[context=="CG" & stage_lineage==i,mean], color="#F37A71", linetype="dashed", alpha=0.75, size=0.75) +
    # geom_hline(yintercept=stats[context=="GC" & stage_lineage==i,mean], color="#00BFC4", linetype="dashed", alpha=0.75, size=0.75) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (scaled)") +
    # coord_cartesian(ylim=c(0,100)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill=FALSE, color=FALSE, linetype=FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(0.8), colour="black"),
      axis.text.y = element_text(size=rel(1.0), colour="black")
    )

  pdf(sprintf("%s/%s_pseudobulk.pdf",io$outdir,i), width=13, height=9)
  print(p)
  dev.off()
}


for (i in unique(data2$class)) {
  
  to.plot <- data2[class==i]
  
  p <- ggplot(to.plot, aes(x=dist, y=rate_scaled, group=context, fill=context, color=context)) +
    geom_ribbon(aes(ymin = rate_scaled - 0.01, ymax = rate_scaled + 0.01), color="black") +
    facet_wrap(~stage_lineage, nrow=2, scales="fixed") +
    labs(x="Distance from center (bp)", y="Met/Acc levels (scaled)") +
    xlim(-opts$window_size, opts$window_size) +
    guides(fill=FALSE, color=FALSE, linetype=FALSE) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(0.8), colour="black"),
      axis.text.y = element_text(size=rel(1.0), colour="black")
    )
  
  pdf(sprintf("%s/%s_pseudobulk.pdf",io$outdir,i), width=11, height=7)
  print(p)
  dev.off()
}

