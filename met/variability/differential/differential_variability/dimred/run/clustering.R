suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(clues))

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

dt <- sample_metadata %>% 
  .[,c("id_met","lineage10x_2")] %>% setnames("id_met","sample") %>%
  merge(Z.mofa,by="sample")

################
## Clustering ##
################

# Define true labels
true_clusters <- as.factor( dt$lineage10x_2 )
ntrue_clusters <- length(unique(true_clusters))

# Run k-means clustering on the latent space defined by MOFA
clustering <- cluster_samples(model, k=ntrue_clusters, iter.max=25)$cluster %>% as.factor
dt <- dt %>% merge( data.table(sample=names(clustering), cluster=clustering), by="sample")

# Display summary table of clustering results vs cell type labels
# table(dt$cluster, true_clusters)

########################################################################
## Plot dimensionality reduction coloured by predicted cluster labels ##
########################################################################

# p <- plot_dimred(dt, color.by="cluster")

# pdf(sprintf("%s_cluster.pdf", args$outprefix), useDingbats=F, width=5.5, height=4)
# print(p)
# dev.off()

###############################
## Get clustering statistics ##
###############################

tmp <- data.table(
  anno = paste(args$anno, collapse=" "),
  hvg = args$hvg,
  # r2 = model@cache[["variance_explained"]]$r2_total[[1]][[1]],
  adj_rand_index = as.numeric(NA),
  jaccard_index = as.numeric(NA),
  silhouette = as.numeric(NA),
  purity = as.numeric(NA)
  # prediction_accuracy = as.numeric(NA)
)

# Calculate silhouette (from 0 to 1)
sh <- silhouette(as.numeric(dt$cluster), dist(Z.mofa[,-"sample"]))
tmp$silhouette <- summary(sh)$avg.width %>% round(3)

# Calculate purity (from 0 to 1)
tmp$purity <- ClusterPurity(true_clusters, dt$cluster) %>% round(3)

# Calculate rand index 
tmp$adj_rand_index <- adjustedRand(cl1=as.numeric(true_clusters), cl2=as.numeric(dt$cluster), randMethod="HA")

# Calculate jaccard index 
tmp$jaccard_index <- adjustedRand(cl1=as.numeric(true_clusters), cl2=as.numeric(dt$cluster), randMethod="Jaccard")

# Random Forest prediction
# df <- get_factors(model, as.data.frame = T) %>% as.data.table %>%
#   merge(sample_metadata, by="sample") %>%
#   dcast(sample+lineage10x_2~factor, value.var="value") %>%
#   .[,lineage10x_2:=as.factor(lineage10x_2)] %>%
#   tibble::column_to_rownames("sample")
# 
# rf <- randomForest(formula=lineage10x_2~., ntree=250, data=df)
# accuracy <- mean(rf$predicted==df$lineage10x_2)
# 
# tmp$prediction_accuracy <- round(accuracy,3)

#################
## Save output ##
#################

fwrite(tmp, sprintf("%s_clustering.txt", args$outprefix), quote=TRUE)

