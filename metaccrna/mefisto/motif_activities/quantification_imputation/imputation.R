suppressMessages(library(RColorBrewer))
suppressMessages(library(MOFA2))

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else if (grepl("argelagr",Sys.info()['nodename'])) {
    source("/Users/argelagr/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else {
  stop()
}

# Options
opts$number_samples_to_mask <- c(100,150,200,250) %>% as.character
opts$seed <- c(1,2,3) %>% as.character

# opts$number_samples_to_mask <- c(100) %>% as.character
# opts$seed <- c(1) %>% as.character


# I/O
io$original.data <- paste0(io$basedir,"/metaccrna/mefisto/quantification_imputation/original_data.rds")
io$mefisto.models.dir <- paste0(io$basedir,"/metaccrna/mefisto/quantification_imputation")
io$outdir <- paste0(io$basedir,"/metaccrna/mefisto/quantification_imputation")

# tmp <- expand.grid(opts$seed,opts$number_samples_to_mask)
# io$mefisto.models <- sprintf("%s/mefisto_imputation_N%s_seed%s.rds",io$mefitsto.models.dir,tmp[,2],tmp[,1])

########################
## Load original data ##
########################

original_data <- readRDS(io$original.data)

################
## Load model ##
################

mefisto_models <- list()
for (i in opts$number_samples_to_mask) {
  mefisto_models[[i]] <- list()
  for (j in opts$seed) {
    mefisto_models[[i]][[j]] <-  readRDS(sprintf("%s/mefisto_imputation_N%s_seed%s.rds",io$mefisto.models.dir,i,j))
  }
}

# Get covariates
# covariates.dt <- get_covariates(mefisto, as.data.frame = T) %>% as.data.table %>% dcast(sample~covariate)

################
## Imputation ##
################

views.to.impute <- c("motif_acc", "motif_met")

# i <- opts$number_samples_to_mask[1]
# j <- opts$seed[1]
results.dt <- opts$number_samples_to_mask %>% map(function(i) { 
  opts$seed %>% map(function(j) {
  
    samples_masked <- mefisto_models[[i]][[j]]$samples_masked
    model <- mefisto_models[[i]][[j]]$model
    
    # Standard MOFA imputation
    file <- sprintf("%s/imputed_data_mofa_N%s_seed%s.rds",io$outdir,i,j)
    if (file.exists(file)) {
      imputed_data_mofa <- readRDS(file)
    } else {
      model <- impute(model, views=views.to.impute)
      imputed_data_mofa <- model@imputed_data
      saveRDS(imputed_data_mofa, file)
    }
    
    # MEFISTO imputation using GP means
    file <- sprintf("%s/imputed_data_mefisto_N%s_seed%s.rds",io$outdir,i,j)
    if (file.exists(file)) {
      imputed_data_mefisto <- readRDS(file)
    } else {
      model <- interpolate_factors(model, model@covariates[[1]])
      Z_interpol <- t(get_interpolated_factors(model, only_mean = TRUE)[[1]]$mean)
      imputed_data_mefisto <- list()
      for (k in views.to.impute) {
        imputed_data_mefisto[[k]] <- tcrossprod(Z_interpol, get_weights(model,k)[[1]]) %>% t
        colnames(imputed_data_mefisto[[k]]) <- colnames(model@data[[k]][[1]])
      }
      saveRDS(imputed_data_mefisto, sprintf("%s/imputed_data_mefisto_N%s_seed%s.rds",io$outdir,i,j))
    }
    
    # kNN imputation
    file <- sprintf("%s/imputed_data_knn_N%s_seed%s.rds",io$outdir,i,j)
    if (file.exists(file)) {
      imputed_data_knn <- readRDS(file)
    } else {
      imputed_data_knn <- list()
      for (k in views.to.impute) {
        imputed_data_knn[[k]] <- smoother_aggregate_nearest_nb(
          mat = get_data(model, views=k)[[1]][[1]], 
          D = pdist(t(model@covariates[[1]])), 
          k = 15
        )
        colnames(imputed_data_knn[[k]]) <- colnames(model@data[[k]][[1]])
      }
      saveRDS(imputed_data_knn, file)
    }
    
    # Store results
    # views.to.impute %>% map(function(k) {
    #   data.table(
    #     N = i,
    #     seed = j,
    #     view = k,
    #     mefisto = mean(abs(imputed_data_mefisto[[k]][,samples_masked] - original_data[[k]][[1]][,samples_masked]),na.rm=T),
    #     mofa = mean(abs(imputed_data_mofa[[k]][[1]][,samples_masked] - original_data[[k]][[1]][,samples_masked]),na.rm=T),
    #     knn = mean(abs(imputed_data_knn[[k]][,samples_masked] - original_data[[k]][[1]][,samples_masked]),na.rm=T)
    #   )
    # }) %>% rbindlist %>% return
    views.to.impute %>% map(function(view) {
        rownames(original_data[[view]][[1]]) %>% map(function(feature_id) {
        data.table(
          view = view,
          feature = feature_id,
          id_met = samples_masked,
          original = original_data[[view]][[1]][feature_id,samples_masked],
          mefisto = imputed_data_mefisto[[view]][feature_id,samples_masked],
          knn = imputed_data_knn[[view]][feature_id,samples_masked],
          mofa = imputed_data_mofa[[view]][[1]][feature_id,samples_masked],
          random = sample(imputed_data_mofa[[view]][[1]][feature_id,], size = length(samples_masked), replace=TRUE)
        ) %>% melt(id.vars=c("feature","view","id_met","original"), variable.name="method") %>%
          .[,.(mae=mean(abs(value-original),na.rm=T)),by=c("method","view","feature")] %>%
          .[,c("N","seed") := list(i,j)]
      }) %>% rbindlist
    }) %>% rbindlist
    
    
  }) %>% rbindlist
}) %>% rbindlist

# Load precomputed results
fwrite(results.dt, sprintf("%s/imputation_evaluation.txt.gz",io$outdir))
results.dt <- fread(sprintf("%s/imputation_evaluation.txt.gz",io$outdir))
       
##########
## Plot ##
##########

to.plot <- results.dt

to.plot[,facet:=sprintf("%s (N=%s)",c("motif_met"="Motif methylation", "motif_acc"="Motif accessibility")[view], N)]
p <- ggboxplot(to.plot[seed=="1"], x="method", y="mae", fill="method", group="seed", outlier.shape=NA) + # add = "dotplot"
  facet_wrap(~facet, nrow=2) +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim=c(0,0.4)) +
  labs(x="", y="Mean absolute error") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(size=rel(0.75)),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(size=rel(0.75)),
    # axis.ticks.x = element_blank(),
    legend.text = element_text(size=rel(0.8))
  )


pdf(sprintf("%s/imputation_evaluation.pdf",io$outdir), width=10, height=6)
print(p)
dev.off()

