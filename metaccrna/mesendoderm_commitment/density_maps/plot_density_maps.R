library(data.table)
library(purrr)
library(ggplot2)

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation/metaccrna/differential/load_settings.R")
io$outdir <- "/Users/ricard/gastrulation/metaccrna/mesendoderm_commitment/density_maps/out"

opts$comparisons <- c(
    "E6.5Epiblast_vs_E6.5E7.5Primitive_Streak",
    "E6.5E7.5Primitive_Streak_vs_E7.5Mesoderm",
    "E6.5E7.5Primitive_Streak_vs_E7.5Endoderm"
  )

###################################
## Select lineage-defining sites ##
###################################

# Lineage-defining elements are defined as ChIP-seq peaks that show differential activity during germ layer commitment

io$diff.met <- "/Users/ricard/data/gastrulation/met/results/differential"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/results/differential"

opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.acc.diff <- 5
opts$min.met.diff <- 5

source("/Users/ricard/gastrulation/metaccrna/mesendoderm_commitment/density_maps/load_data.R")

###############################################
## Load differential DNA methylation results ##
###############################################

# Load precomputed differential results
diff.met.ps <- lapply(opts$comparisons, function(i) 
  lapply(names(opts$met.annos), function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$diff.met,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]

######################################################
## Load differential chromatin accessiblity results ##
######################################################

diff.acc.ps <- lapply(opts$comparisons, function(i) 
  lapply(names(opts$acc.annos), function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$diff.acc,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]

###################################
## Subset lineage-defining sites ##
###################################

# Methylation
diff.met.ps <- diff.met.ps %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id%in%diff.met[sig==T & anno==y,id]]) %>%
  rbindlist

# Accessibility
diff.acc.ps <- diff.acc.ps %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id%in%diff.acc[sig==T & anno==y,id]]) %>%
  rbindlist

#################################################
## Merge methylation and accessibility results ##
#################################################

diff.met.ps <- diff.met.ps %>% 
  .[,anno:=stringr::str_replace_all(anno,opts$met.annos)] %>%
  .[,c("id","anno","diff","sig","comparison")]
  
diff.acc.ps <- diff.acc.ps %>% 
  .[,anno:=stringr::str_replace_all(anno,opts$acc.annos)] %>%
  .[,c("id","anno","diff","sig","comparison")]

diff.metacc <- rbind(
  diff.met.ps[,assay:="met"], 
  diff.acc.ps[,assay:="acc"]
) %>% data.table::dcast(id+anno+comparison~assay, value.var=c("diff","sig"))


#############################
## Bivariate density plots ##
#############################

opts$ymin <- -30
opts$ymax <- 30
opts$xmin <- -55
opts$xmax <- 55

to.plot <- diff.metacc %>%  .[complete.cases(.)]

p <- ggplot(to.plot, aes(x=diff_met, y=diff_acc)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
  # scale_fill_viridis_c(option="inferno", end=.9, direction=-1) +
  scale_fill_gradient(low="white", high="purple") +
  # scale_fill_distiller(palette = "Spectral") +
  
  stat_density2d(contour=T, size=0.2, alpha=0.8, color="black") +
  
  geom_segment(aes(x=opts$xmin, xend=opts$xmax, y=0, yend=0), size=0.25, color="black") +
  geom_segment(aes(x=0, xend=0, y=opts$ymin, yend=opts$ymax), size=0.25, color="black") +
  
  coord_cartesian(xlim=c(opts$xmin,opts$xmax), ylim=c(opts$ymin,opts$ymax)) +
  
  facet_wrap(~comparison+anno, ncol=3, scales="fixed") +
  
  labs(x="Differential methylation rate", y="Differential accessibility rate") +
  
  theme(
    axis.text.x = element_text(size=rel(1.1), color='black'),
    axis.text.y = element_text(size=rel(1.1), color='black'),
    axis.title.x = element_text(size=rel(1.2), color='black'),
    axis.title.y = element_text(size=rel(1.2), color='black'),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_blank()
  )

# pdf(sprintf("%s/mesendoderm_density_enhancers.pdf",io$outdir), width=10, height=11, useDingbats = F)
print(p)
# dev.off()
