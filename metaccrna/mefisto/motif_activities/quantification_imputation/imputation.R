suppressMessages(library(RColorBrewer))
suppressMessages(library(MOFA2))

###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/metaccrna/mefisto/motif_activities/load_settings.R")
} else {
  stop()
}

# Options
opts$number_samples_to_mask <- c(100)

# I/O
io$original.data <- paste0(io$outdir,"/original_data.rds")
io$mefitsto.models.dir <- paste0(io$basedir,"/metaccrna/mefisto/quantification_imputation")
io$outdir <- paste0(io$basedir,"/metaccrna/mefisto/quantification_imputation")
io$mefisto.models <- sprintf("%s/mefisto_imputation_N%s.rds",io$mefitsto.models.dir,opts$number_samples_to_mask)

########################
## Load original data ##
########################

original_data <- readRDS(io$original.data)

################
## Load model ##
################

# mefisto <- readRDS(io$mefisto.model)
mefisto_models <- io$mefisto.models %>% map(readRDS)
names(mefisto_models) <- as.character(opts$number_samples_to_mask)

# Get covariates
# covariates.dt <- get_covariates(mefisto, as.data.frame = T) %>% as.data.table %>% dcast(sample~covariate)

################
## Imputation ##
################

views.to.impute <- c("motif_acc", "motif_met")

for (i in 1:length(mefisto_models)) {
  
  samples_masked <- mefisto_models[[i]]$samples_masked
  model <- mefisto_models[[i]]$model
  
  # Standard MOFA imputation
  model <- impute(model, views=views.to.impute)
  imputed_data_mofa <- model@imputed_data[]
  
  # MEFISTO imputation using GP means
  model <- interpolate_factors(model, model@covariates[[1]])
  Z_interpol <- t(get_interpolated_factors(model, only_mean = TRUE)[[1]]$mean)
  imputed_data_mefisto <- list()
  for (j in views.to.impute) {
    # model@imputed_data[[view]][[1]] <- tcrossprod(Z_interpol,get_weights(model, view)[[1]]) %>% t
    imputed_data_mefisto[[j]] <- tcrossprod(Z_interpol,get_weights(model, j)[[1]]) %>% t
    colnames(imputed_data_mefisto[[j]]) <- colnames(model@data[[j]][[1]])
  }
  
  # mean imputation
  imputed_data_mean <- list()
  for (j in views.to.impute) {
    tmp <- get_data(model, views=j)[[1]][[1]]
    imputed_data_mean[[j]] <- apply(tmp, 1, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }) %>% t
  }
  
  
  # kNN imputation
  imputed_data_knn <- list()
  for (j in views.to.impute) {
    imputed_data_knn[[j]] <- smoother_aggregate_nearest_nb(
      mat = get_data(model, views=j)[[1]][[1]], 
      D = pdist(t(model@covariates[[1]])), 
      k = 25
    )
    colnames(imputed_data_knn[[j]]) <- colnames(model@data[[j]][[1]])
  }
}


###################
## Calculate MSE ##
###################

# foo <- data.table(
#   sample = unlist(samples_names(mefisto)),
#   imputed = colMeans(data_imputed,na.rm=T),
#   non_imputed = colMeans(data,na.rm=T),
#   interpolated = colMeans(data_imputed_GP,na.rm=T)
# )

i <- "motif_met"
N <- 50
data.table(
  N = N,
  mefisto = sum(abs(imputed_data_mefisto[[i]][,samples_masked] - original_data[[i]][[1]][,samples_masked]),na.rm=T),
  mofa = sum(abs(imputed_data_mofa[[i]][[1]][,samples_masked] - original_data[[i]][[1]][,samples_masked]),na.rm=T),
  mean = sum(abs(imputed_data_mean[[i]][,samples_masked] - original_data[[i]][[1]][,samples_masked]),na.rm=T),
  knn = sum(abs(imputed_data_knn[[i]][,samples_masked] - original_data[[i]][[1]][,samples_masked]),na.rm=T)
)

##########
## Plot ##
##########

for (view in c("motif_met","motif_acc")) {
  
  features.to.plot <- features_names(mefisto)[view]
  # features.to.plot <- list(paste0("MSGN1_412",gsub("motif","",view))); names(features.to.plot) <- view
  data <- get_data(mefisto, views=view, features=features.to.plot)[[1]][[1]]
  data_imputed <- get_imputed_data(mefisto, views=view, features=features.to.plot)[[1]][[1]]
  data_imputed_mefisto <- imputed_data_mefisto[[view]][features.to.plot[[1]],,drop=F]
  data_imputed_knn <- imputed_data_knn[[view]][features.to.plot[[1]],,drop=F]
  
  for (i in features.to.plot[[view]]) {
  # for (i in paste0("FOXD1_2",gsub("motif","",view))) {
    print(i)
    
    to.plot <- data.table(
      sample = unlist(samples_names(mefisto)),
      non_imputed = data[i,],
      # imputed = data_imputed[i,],
      knn = data_imputed_knn[i,],
      interpolated = data_imputed_mefisto[i,]
    ) %>% 
      # setnames(c("sample",sprintf("%s (original)",i),sprintf(" %s (imputed)",i))) %>%
      setnames(c("sample",sprintf("%s (original)",i),sprintf(" %s (kNN imputed)",i),sprintf(" %s (MEFISTO imputed)",i))) %>%
      melt(id.vars="sample", value.name="value") %>% 
      merge(covariates.dt, by="sample")
    
    max.value <- max(to.plot[variable==sprintf(" %s (MEFISTO imputed)",i),value])
    min.value <- min(to.plot[variable==sprintf(" %s (MEFISTO imputed)",i),value])
    to.plot[value>max.value,value:=max.value]
    to.plot[value<min.value,value:=min.value]
    
    # to.plot[,beta_value:=100*2**m_value/(1+2**m_value)]
    
    p <- ggplot(to.plot, aes(x=V1, y=V2, fill=value)) +
      geom_point(alpha=0.9, size=2.0, shape=21, stroke=0.15) +
      facet_wrap(~variable) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.8))
      )
    
    # p <- p + scale_colour_gradientn(colours = terrain.colors(10))
    if (view=="motif_met") {
      p <- p + scale_fill_gradient2(low = "blue", mid="gray90", high = "red")
      # p <- p + scale_fill_gradientn(colours = brewer.pal(9, "OrRd"))
    } else if (view=="motif_acc") {
      p <- p + scale_fill_gradient2(low = "yellow", mid="gray90", high = "purple")
      # p <- p + scale_fill_gradientn(colours = rev(brewer.pal(9, "Blues")))
    }
    
    pdf(sprintf("%s/imputation/knn/%s_imputed.pdf",io$outdir,i), width=8, height=3)
    print(p)
    dev.off()
  }
}



###############################
## Compare to kNN imputation ##
###############################

# kNN imputation
mtx <- get_data(mefisto, views="motif_met")[[1]][[1]]
mtx.smoothed <- smoother_aggregate_nearest_nb(mat=mtx, D=pdist(t(mefisto@covariates[[1]])), k=25)
