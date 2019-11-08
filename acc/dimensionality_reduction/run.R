library(data.table)
library(purrr)
library(MOFA)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation/acc/dimensionality_reduction/load_settings.R")
} else {
  source("/homes/ricard/gastrulation/acc/dimensionality_reduction/load_settings.R")
}

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(names(opts$annos), function(n) 
  fread(sprintf("%s/%s.tsv.gz",io$data.dir,n), select=c(1,2,3,5,6)) %>% 
    setnames(c("id_acc","id","anno","N","rate")) %>% .[N>=opts$min.GpCs] %>% .[,N:=NULL]
  ) %>% rbindlist 

####################################################
##  Merge accessibility data with sample metadata ##
####################################################

acc_dt <- acc_dt %>% merge(sample_metadata, by="id_acc") %>% droplevels()

#######################################
## Calculate M value from Beta value ##
#######################################

acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#################
## Filter data ##
#################

# Filter features by minimum number of GpC sites
acc_dt <- acc_dt

# Filter features by coverage
nsamples <- length(unique(acc_dt$id_acc))
acc_dt <- acc_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

############################
## Regress out covariates ##
############################

# Global accessibility rate
foo <- acc_dt[,.(mean=mean(m)),by=c("id_acc","anno")]
acc_dt <- acc_dt %>% merge(foo, by=c("id_acc","anno")) %>%
  .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")] %>%
  .[,mean:=NULL]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

#######################################
## Create matrix from the data.table ##
#######################################

dmatrix_list <- acc_dt %>% split(.$anno) %>% 
  map(~
    dcast(.[,c("id_acc","id","m")], formula=id_acc~id, value.var="m") %>% as.data.frame %>% 
      tibble::column_to_rownames("id_acc") %>% as.matrix %>% t
     )
lapply(dmatrix_list,dim)

####################
## Fit MOFA model ##
####################

MOFAobject <- createMOFAobject(dmatrix_list)

# Set options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 2

TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 1000
TrainOptions$tolerance <- 0.10
TrainOptions$DropFactorThreshold <- 0.00
TrainOptions$seed <- 42

# Prepare
MOFAmodel <- prepareMOFA(
  MOFAobject,
  ModelOptions = ModelOptions, 
  TrainOptions = TrainOptions
)

# Train the model
model <- runMOFA(MOFAmodel, io$outfile)
