## ----load_modules, echo=FALSE, include=FALSE----------------------------------

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/multiome_peaks/load_settings.R")
}

opts$test <- TRUE

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[,c("sample","id_acc","id_met","stage","lineage10x_2")] %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[,c("lineage10x_2"):=NULL] %>%
  .[id_met%in%opts$met.cells | id_acc%in%opts$acc.cells]

sample_metadata[,c("stage","stage_lineage"):=list(factor(stage,levels=opts$stages), factor(stage_lineage,levels=opts$stage_lineage))]
sample_metadata[,c("sample","id_acc","id_met"):=list(as.factor(sample), as.factor(id_acc), as.factor(id_met))]

# Merge E5.5 and E6.5 epiblast
# sample_metadata %>% .[,stage_lineage:=ifelse(stage_lineage=="E5.5 Epiblast","E6.5 Epiblast",stage_lineage)]

# Subset cells
if (opts$test) {
  opts$ncells <- 5
  opts$filt.cells <- sample_metadata[,head(unique(sample),n=opts$ncells),by="stage_lineage"] %>% .$V1
  sample_metadata <- sample_metadata[sample %in% opts$filt.cells]
  opts$met.cells <- sample_metadata$id_met
  opts$acc.cells <- sample_metadata$id_acc
}


##############################
## Load genomic annotations ##
##############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_annotations.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_annotations.R")
}
anno_df.met <- anno_df
anno_df.acc <- anno_df


# Load celltype marker peaks
multiome_marker_peaks.dt <- fread(io$multiome.atac.marker_peaks)

#########################
## Load scNMT-seq data ##
#########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_data.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metacc/pseudobulk_profiles/load_data.R")
}

# Add sample metadata
met <- met %>% merge(sample_metadata, by="id_met") %>% droplevels()
acc <- acc %>% merge(sample_metadata, by="id_acc") %>% droplevels()

# Merge
data <- rbind(
  met[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")],
  acc[,c("sample","stage","stage_lineage","id","anno","dist","rate","context")]
)
# data[,rate:=rate*100]

# Rename genomic annotations
# data[,anno:=stringr::str_replace_all(anno,opts$annos)]


# Load precomputed data
# saveRDS(data, "/Users/ricard/data/gastrulation/metacc/pseudobulked_profiles/lineage_enhancers/data.rds")
# data <- readRDS("/Users/ricard/data/gastrulation/metacc/pseudobulked_profiles/lineage_enhancers/data.rds")


#################################################################
## Load genome-wide global methylation and accessibility stats ##
#################################################################

met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>%
  merge(sample_metadata[,.(sample,id_met)], by="id_met") %>% .[,context:="CG"]

acc.stats <- fread(io$acc.stats) %>% .[,c("id_acc","mean")] %>%
  merge(sample_metadata[,.(sample,id_acc)], by="id_acc") %>% .[,context:="GC"]

stats <- rbind(
  met.stats[,c("sample","mean","context")],
  acc.stats[,c("sample","mean","context")]
) %>% merge(sample_metadata[,c("sample","stage","stage_lineage")],by="sample") %>%
  .[,.(mean=mean(mean)),by=c("stage_lineage","context")]


##########
## Plot ##
##########

# data[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]
# stats[,stage_lineage:=stringr::str_replace_all(stage_lineage,"_"," ")]

p_list <- list()

for (i in opts$stage_lineage) {
  print(i)
  
  tmp <- data[stage_lineage==i]
  
  p_list[[i]] <- ggplot(tmp, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
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
  # print(p_list[[i]])

  pdf(file=sprintf("%s/%s.pdf",io$pdfdir,i), width=8.5, height=5)
  print(p_list[[i]])
  dev.off()
}

