library(RColorBrewer)

#######################
## Plotting settings ##
#######################

theme_pub <- function() {
  theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.2)),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# min/max M-value for DNA methylation
opts$min_met <- (-3.5)   
opts$max_met <- 3.5

# min/max M-value for chromatin accessibility
opts$min_acc <- (-2.5)
opts$max_acc <- 0.5

##################################
## Run dimensionality reduction ##
##################################

MOFAobject <- run_tsne(MOFAobject, factors="all")
dimred <- as.data.table(MOFAobject@dim_red$TSNE) %>%
  setnames(c("sample","V1","V2"))

###############
## Plot UMAP ##
###############

# p <- plot_dimred(MOFAobject, method="TSNE", color_by="technology", stroke = 0.2, dot_size = 2.5) +
#   scale_fill_brewer(palette="Dark2") +
p <- plot_dimred(MOFAobject, method="TSNE", color_by="stage", stroke = 0.2, dot_size = 2.5) +
  scale_fill_brewer(palette="Set1") +
  labs(x="", y="") +
  theme_pub() +
  theme(
    legend.position = "top"
  )
# print(p)

pdf(sprintf("%s/pdf/umap_stage.pdf",io$outdir), width=3.5, height=4, useDingbats = F)
print(p)
dev.off()


#################################
## Plot methylation imputation ##
#################################

factor <- 5
view <- "met_genebody"
features <- names(head(abs(sort(get_weights(MOFAobject, views=view, factor=factor)[[1]][,1])), n=100))

# met <- unlist(lapply(MOFAobject@imputed_data[[view]], function(x) colMeans(x[features,], na.rm=T)))
met <- unlist(lapply(MOFAobject@data[[view]], function(x) colMeans(x[features,], na.rm=T)))
met[met<(opts$min_met)] <- opts$min_met
met[met>(opts$max_met)] <- opts$max_met

# Convert M-values to B-values
met <- 100*2**met/(1+2**met)

to.plot <- dimred %>% 
  merge(data.table(sample = samples_names(MOFAobject)[[1]], met = met), by="sample")

# Subsample
to.plot[sample(nrow(to.plot), size=1200),met:=NA]

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(fill="grey", alpha=0.4, size=1, data=to.plot[is.na(met)]) +
  geom_point(aes(fill=met), shape=21, stroke=0.15, size=2, data=to.plot[!is.na(met)]) +
  scale_fill_gradientn(colours = brewer.pal(9, "OrRd"), limits=c(0,100)) +
  labs(x="", y="") +
  theme_pub() 
# print(p)

pdf(sprintf("%s/pdf/met_not_imputed.pdf",io$outdir), width=4, height=3.5, useDingbats = F)
print(p)
dev.off()


###################################
## Plot accessibility imputation ##
###################################

factor <- 5
view <- "acc_genebody"
features <- names(head(abs(sort(get_weights(MOFAobject, views=view, factor=factor)[[1]][,1])), n=100))

acc <- unlist(lapply(MOFAobject@data[[view]], function(x) colMeans(x[features,], na.rm=T)))
# acc <- unlist(lapply(MOFAobject@imputed_data[[view]], function(x) colMeans(x[features,], na.rm=T)))
acc[acc<(opts$min_acc)] <- opts$min_acc
acc[acc>(opts$max_acc)] <- opts$max_acc

# Convert M-values to B-values
acc <- 100*2**acc/(1+2**acc)

to.plot <- dimred %>% 
  merge(data.table(sample = samples_names(MOFAobject)[[1]], acc = acc), by="sample")

# Subsample
to.plot[sample(nrow(to.plot), size=1200),acc:=NA]

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(fill="grey", alpha=0.4, size=1, data=to.plot[is.na(acc)]) +
  geom_point(aes(fill=acc), shape=21, stroke=0.15, size=2, data=to.plot[!is.na(acc)]) +
  scale_fill_gradientn(colours = rev(brewer.pal(9, "Blues")), limits=c(10,60)) +
  labs(x="", y="") +
  theme_pub() 
# print(p)

pdf(sprintf("%s/pdf/acc_not_imputed.pdf",io$outdir), width=4, height=3.5, useDingbats = F)
print(p)
dev.off()