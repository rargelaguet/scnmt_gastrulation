    # mean imputation
    file <- sprintf("%s/imputed_data_mean_N%s_seed%s.rds",io$outdir,i,j)
    if (file.exists(file)) {
      imputed_data_mean <- readRDS(file)
    } else {
      imputed_data_mean <- list()
      for (k in views.to.impute) {
        tmp <- get_data(model, views=k)[[1]][[1]]
        imputed_data_mean[[k]] <- apply(tmp, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }) %>% t
      }
      saveRDS(imputed_data_mean, file)
    }







    stop()

#################################################
## Compare predictions to pseudobulk estimates ##
#################################################

met.dt <- fread("/Users/ricard/data/gastrulation/met/feature_level/motifs/prom_2000_2000_jaspar2020_MSGN1_412.tsv.gz") %>%
  setnames(c("id_met","feature_id","anno","met_reads","total_reads","rate"))

sample_metadata <- fread(io$metadata) %>% 
  .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[pass_metQC==T & stage_lineage%in%opts$stage_lineage] %>%
  droplevels


umap.dt <- fread(io$umap) %>%
  .[sample%in%unique(data$sample)]

sample_metadata <- sample_metadata %>% merge(umap.dt,by="sample")

data <- data[sample%in%unique(umap.dt$sample)]
  
met.dt <- met.dt %>% merge(sample_metadata[,c("stage_lineage","id_met")], by="id_met")
met_pseudobulk.dt <- met.dt %>%
  .[,.(met_reads=sum(met_reads), total_reads=sum(total_reads)),by="stage_lineage"] %>%
  .[,rate:=100*(met_reads/total_reads)] %>%
  .[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] 

hist(original_data$motif_met$single_group["MSGN1_412_met",])
grep("MSGN",rownames(original_data$motif_met$single_group), value=T)


imput.knn <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/quantification_imputation/imputed_data_knn_N100_seed1.rds")
imput.mefisto <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/quantification_imputation/imputed_data_mefisto_N100_seed1.rds")
imput.mean <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/quantification_imputation/imputed_data_mean_N100_seed1.rds")
imput.mofa <- readRDS("/Users/ricard/data/gastrulation/metaccrna/mefisto/quantification_imputation/imputed_data_mofa_N100_seed1.rds")

grep("CREM",rownames(original_data$motif_met[[1]]),value=T)
feature_id <- "HNF1A_85_met"

# to.plot <- data.table(
#   id_met = colnames(original_data$motif_met[[1]]),
#   original = original_data$motif_met[[1]][feature_id,],
#   mefisto = imput.mefisto$motif_met[feature_id,],
#   mean = imput.mean$motif_met[feature_id,],
#   knn = imput.knn$motif_met[feature_id,],
#   mofa = imput.mofa$motif_met[[1]][feature_id,]
# ) %>% 
#   melt(id.vars=c("id_met","original"), variable.name="method") %>%
#   merge(sample_metadata[,c("stage_lineage","id_met")], by="id_met") 

features.to.plot <- rownames(original_data$motif_met[[1]])# %>% head(n=25)
# features.to.plot <- rownames(original_data$motif_acc[[1]])# %>% head(n=25)

view <- "motif_met"
results_per_feature.dt <- features.to.plot %>% map(function(i) {
  data.table(
    feature = i,
    id_met = samples_masked,
    view = view,
    original = original_data[[view]][[1]][i,samples_masked],
    mefisto = imput.mefisto[[view]][i,samples_masked],
    mean = imput.mean[[view]][i,samples_masked],
    knn = imput.knn[[view]][i,samples_masked],
    mofa = imput.mofa[[view]][[1]][i,samples_masked],
    random = sample(imput.mofa[[view]][[1]][i,], size = length(samples_masked), replace=TRUE)
  ) %>% melt(id.vars=c("feature","view","id_met","original"), variable.name="method") %>%
  .[,.(mae=mean(abs(value-original),na.rm=T)),by=c("method","view","feature")]
}) %>% rbindlist


to.plot <- results_per_feature.dt# %>% 
  # merge(sample_metadata[,c("stage_lineage","id_met")], by="id_met") 

ggboxplot(to.plot, x="method", y="mae", fill="method") +
  facet_wrap(~view) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="Number of sampled masked", y="Mean absolute error") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75)),
    legend.text = element_text(size=rel(0.8))
  )
  
# ggscatter(to.plot, x="original", y="value", color="stage_lineage") +
#   geom_abline(slope=1, intercept=0) +
#   scale_color_manual(values=opts$stagelineage.colors) +
#   facet_wrap(~method, scales="fixed")