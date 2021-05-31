###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/met/results/motifs/umap")


# Options

opts$annos <- c("multiome_peaks") # H3K27ac_distal_E7.5_union_intersect12
opts$motif.annotation <- "jaspar2020"

opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
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

opts$min.observations <- 25

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[pass_metQC==TRUE & stage_lineage%in%opts$stage_lineage] %>%
  .[,lineage10x:=stringr::str_replace_all(lineage10x,opts$to.merge)] %>%
  .[,c("sample","id_met","stage","lineage10x","lineage10x_2","stage_lineage")]
  # .[,lineage10x_2:=stringr::str_replace_all(lineage10x_2,"_"," ")] %>% 

#################
## Load motifs ##
#################

motifOverlap.se <- readRDS(sprintf("%s/%s_%s.rds",io$motifs.dir,opts$anno,opts$motif.annotation))
opts$motifs <- colnames(motifOverlap.se$motifMatches)# %>% head(n=2)

###############
## Load data ##
###############

# met_dt <- opts$annos %>% map(function(i) {
#   opts$motifs %>% map(function(j) {
#     file <- sprintf("%s/%s_%s_%s.tsv.gz",io$met_data_motifs,i,opts$motif.annotation,j)
#     if (file.exists(file)) {
#       fread(file, header=F, colClasses = c("factor","character","factor","integer","integer","integer")) %>%
#         .[V1%in%sample_metadata$id_met] %>%
#         setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>%
#         .[Ntotal>=opts$min.CpGs] %>%
#         .[,.(rate=100*(sum(Nmet)/sum(Ntotal)), Nmet=sum(Nmet), N=sum(Ntotal)),by=c("id_met","anno")]
#     }
#   }) %>% rbindlist
# }) %>% rbindlist

# Save
# fwrite(met_dt, paste0(io$outdir,"/motif_met.txt.gz"), sep="\t")

# Load pre-computed data
io$motifs.met <- paste0(io$basedir,"/met/results/motifs/motif_met.txt.gz")
met_dt <- fread(io$motifs.met) %>%
  .[N>=opts$min.observations] %>%
  .[,motif:=gsub(sprintf("%s_%s_",opts$annos,opts$motif.annotation),"",anno)]

###########################
## Load UMAP coordinates ##
###########################

io$umap <- paste0(io$basedir,"/metaccrna/mofa/all_stages/umap_coordinates.txt")
umap.dt <- fread(io$umap) %>% setnames("sample","id_met")

sample_metadata <- sample_metadata %>% 
  merge(umap.dt, by="id_met", all.x = TRUE)

##################
## Prepare data ##
##################

# Regress out global rates
foo <- fread(io$met.stats) %>% .[,c("id_met","mean")]
met_dt <- met_dt %>% merge(foo, by="id_met") %>%
  .[,rate:=lm(formula=rate~mean)[["residuals"]], by=c("motif","anno")]

to.plot <- met_dt %>%
  merge(sample_metadata, by="id_met")

##########
## Plot ##
##########

for (i in unique(to.plot$motif)) {

  p <- ggplot(to.plot[motif==i], aes(x=V1, y=V2, color=rate)) +
    geom_point(alpha=0.75, size=2.5) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.8))
    )
  
  # p <- p + scale_colour_gradientn(colours = terrain.colors(10))
  p <- p + scale_color_gradient2(low = "blue", mid="gray90", high = "red")
  # if (view=="motif_met") {
  #   p <- p + scale_colour_gradientn(colours = brewer.pal(9, "OrRd"))
  # } else if (view=="motif_acc") {
  #   p <- p + scale_colour_gradientn(colours = rev(brewer.pal(9, "Blues")))
  # }

  pdf(sprintf("%s/%s_met_umap.pdf",io$outdir,i), width=6, height=5)
  print(p)
  dev.off()
}

