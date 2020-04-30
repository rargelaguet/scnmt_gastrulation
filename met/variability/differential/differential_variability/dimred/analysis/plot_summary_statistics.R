library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$indir <- "/Users/ricard/data/gastrulation/met/results/variability/differential/differential_variability/dimred/txt"
  io$outdir <- "/Users/ricard/data/gastrulation/met/results/variability/differential/differential_variability/dimred"
} else {
  stop("Computer not recognised")
}


#############
## Options ##
#############

opts <- list()

# Define genomic contexts
opts$anno <- c(
  "H3K27ac_distal_E7.5_union_intersect12" = "Distal H3K27ac",
  "H3K27ac_distal_E7.5_union_intersect12_500" = "Distal H3K27ac (+500bp)",
  "prom_2000_2000" = "Promoters"
)

# Number of features
opts$number.hvg <- seq(25,500,by=25)

# Metrics to plot
opts$metrics <- c(
  # "silhouette" = "Cluster silhouette",
  "purity" = "Cluster purity",
  # "prediction_accuracy" = "Prediction accuracy",
  # "r2" = "Variance explained (%)",
  "adj_rand_index" = "Adjusted rand index",
  "jaccard_index" = "Jaccard index"
)

###############################
## Load pre-computed results ##
###############################

foo <- list()
for (i in names(opts$anno)) {
  for (j in opts$number.hvg) {
    file <- sprintf("%s/diffvar_%s_%d_clustering.txt", io$indir,i,j)
    if (file.exists(file)) {
      foo[[file]] <- fread(file,sep=",")
    } else {
      print(file)
    }
  }
}

  

##########
## Plot ##
##########

to.plot <- rbindlist(foo) %>%
  .[,anno:=stringr::str_replace_all(anno, opts$anno)] %>%
  melt(id.vars=c("anno","hvg"), variable.name="metric") %>%
  setorder(anno,hvg) %>%
  .[,value2:=frollmean(value, align="right", n=4), by=c("anno","metric")] %>%
  .[is.na(value2),value2:=value]

for (i in names(opts$metrics)) {
  p <- ggline(to.plot[metric==i], x="hvg", y="value", color = "anno", plot_type="l", size=0.5, point_size=0)  +
    scale_color_brewer(palette="Dark2") +
    # facet_wrap(~metric, nrow=2, scales="fixed") +
    labs(x="Number of HVF", y=opts$metrics[i]) +
    theme(
      legend.position = "right",
      axis.text = element_text(size=rel(0.7))
    )

  pdf(sprintf("%s/%s.pdf",io$outdir,i), width=6, height=4)
  print(p)
  dev.off()
}




# for (i in names(opts$metrics)) {
#   p <- ggplot(to.plot[metric == i], aes(x = hvg, y = value2, colour = model)) +
#     # geom_jitter(size = 2.1, width = 0.15, height = -0.1, shape = 21) +
#     # geom_point(data = dt2, aes(x = x, y = y), position = pd, size = 3, shape = 21, fill = "white") +
#     geom_smooth(aes(fill = model), span = 0.6, method = "loess", se = TRUE, size = 1, alpha = 0.5) +
#     facet_wrap(~anno, nrow = 1, scales = "fixed") +
#     labs(x = "Number of HVF", y = opts$metrics[i], title = NULL) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
#     scale_fill_manual(values = c("#E69F00", "#999999", "#56B4E9")) +
#     scale_colour_manual(values = c("#E69F00", "#999999", "#56B4E9")) +
#     theme_classic() +
#     theme(
#       legend.position = "none", #c(0.9, 0.168), #"right",
#       legend.margin=margin(t=-.1, r=0, b=-0.2, l=0, unit="cm"),
#       legend.title = element_blank(),
#       legend.text = element_text(color = "black", size = rel(1.1)),
#       axis.text = element_text(color = "black", size = rel(0.8)),
#       axis.title = element_text(color = "black", size = rel(1.1))
#     )
#   pdf(sprintf("%s/%s.pdf",io$outdir,i), width = 8, height = 2.5)
#   print(p)
#   dev.off()
# }
