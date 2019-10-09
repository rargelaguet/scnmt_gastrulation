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

io$diff.met <- "/Users/ricard/data/gastrulation/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/gastrulation/acc/differential/feature_level"

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
    fread(cmd=sprintf("zcat < %s/%s_%s.txt.gz",io$diff.met,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]

######################################################
## Load differential chromatin accessiblity results ##
######################################################

diff.acc.ps <- lapply(opts$comparisons, function(i) 
  lapply(names(opts$acc.annos), function(j)
    fread(cmd=sprintf("zcat < %s/%s_%s.txt.gz",io$diff.acc,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]

###################################
## Subset lineage-defining sites ##
###################################

# Methylation
diff.met.ps <- diff.met.ps %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id%in%diff.met[sig==T & anno==y,id]]) %>%
  rbindlist %>% droplevels()

# Accessibility
diff.acc.ps <- diff.acc.ps %>% split(.$anno) %>%
  map2(.,names(.), function(x,y) x[id%in%diff.acc[sig==T & anno==y,id]]) %>%
  rbindlist %>% droplevels()

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

########################################
## Plot fraction of differential hits ##
########################################

# tmp <- diff.metacc %>%
#   .[,.(Nmet=mean(met,na.rm=T), Nacc=mean(acc,na.rm=T)), by=c("anno","comparison")] %>%
#   melt(id.vars=c("anno","comparison"), variable.name="assay", value.name="N")

# for (i in opts$comparisons) {
#   p <- ggplot(tmp[comparison==i], aes(x=anno, y=N, group=assay)) +
#     geom_bar(aes(fill=assay), stat="identity", position="dodge", color="black", size=0.25) +
#     scale_fill_manual(values=c("Nmet"="#F37A71", "Nacc"="#00BFC4")) +
#     labs(x="", y="Fraction of differential sites", title=i) +
#     coord_cartesian(ylim=c(0,0.35)) +
#     theme_pub() + theme(legend.position = "none")
#   
#   # pdf(paste0(io$outdir,"/mes_fractionsigcor.pdf"), width=8, height=5)
#   print(p)
#   # dev.off()
# }

# tmp %>% .[,comparison:=factor(comparison, levels=opts$comparisons)]
# 
# p <- ggplot(tmp, aes(x=comparison, y=N, group=assay)) +
#   geom_bar(aes(fill=assay), stat="identity", position="dodge", color="black", size=0.25) +
#   scale_fill_manual(values=c("Nmet"="#F37A71", "Nacc"="#00BFC4")) +
#   facet_wrap(~anno, nrow=2) +
#   labs(x="", y="Fraction of differential sites") +
#   coord_cartesian(ylim=c(0,0.35)) +
#   theme_pub() + theme(legend.position = "none", axis.text.x = element_blank())
#   
# # pdf(paste0(io$outdir,"/mes_fractionsigcor.pdf"), width=8, height=5)
# print(p)
# # dev.off()

  
#######################################
## Plot number of differential hits  ##
#######################################

# tmp <- rbind(diff.met, diff.acc) %>%
#   .[,.(number_positive_hits=sum(sig==T & diff>0, na.rm=T), 
#        number_negative_hits=-sum(sig==T & diff<0, na.rm=T)), by=c("assay","anno","comparison")] %>%
#   melt(id.vars=c("anno","comparison","assay"))
# 
# 
# ylim <- c(min(tmp$value), max(tmp$value))
# 
# for (i in opts$comparisons) {
#   # p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim) +
#   #   theme(axis.text.x = element_blank())
# 
#   p <- ggplot(tmp[comparison==i], aes(x=anno, y=value)) +
#     geom_bar(aes(fill=assay), color="black", stat="identity", position="dodge", size=0.25) +
#     scale_fill_manual(values=c("met"="#F37A71", "acc"="#00BFC4")) +
#     geom_hline(yintercept=0, color="black") +
#     scale_y_continuous(limits=c(ylim[1],ylim[2])) +
#     facet_wrap(~anno, nrow=1, scales="free_x") +
#     labs(x="", y="Number of hits") +
#     theme_bw() +
#     theme(
#       plot.title = element_text(size=11, face='bold', hjust=0.5),
#       axis.text = element_text(size=rel(0.9), color='black'),
#       axis.text.x = element_blank(),
#       # axis.text.x = element_text(size=rel(0.75), angle=30, hjust=1, vjust=1, color="black"),
#       axis.ticks.x = element_blank(),
#       axis.title = element_text(size=rel(1.0), color='black'),
#       axis.line = element_line(color="black"),
#       legend.position="none"
#     )
#   p
#   
#   pdf(sprintf("%s/%s_barplotsN.pdf",io$outdir,i), width=6, height=2.5)
#   print(p)
#   dev.off()
# }


# a <- diff.met[comparison=="E6.5E7.5Primitive_Streak_vs_E6.5E7.5Mesoderm" & anno=="Mes- enhancers"] %>% .[sig==T]
# b <- diff.met[comparison=="E6.5E7.5Primitive_Streak_vs_E7.5Endoderm" & anno=="Mes- enhancers"] %>% .[sig==T]
# foo <- VennDiagram::venn.diagram(
#   x = list("Mesoderm"=a$id, "Endoderm"=b$id),
#   filename=NULL,
#   col="transparent", fill=c("#CD3278","#43CD80"), alpha = 0.60, cex = 1.5,
#   fontfamily = "serif", fontface = "bold")
# pdf(file=sprintf("%s/venn.pdf",io$outdir))
# grid.draw(foo)
# dev.off()


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
print(p)

pdf(sprintf("%s/mesendoderm_density_enhancers.pdf",io$outdir), width=10, height=11, useDingbats = F)
print(p)
dev.off()
