library(MOFA)
library(data.table)
library(purrr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation/met/dimensionality_reduction/load_settings.R")
} else {
  source("/homes/ricard/gastrulation/met/dimensionality_reduction/load_settings.R")
}

###############################
## Load DNA methylation data ##
###############################

met_dt <- lapply(names(opts$annos), function(n)
  fread(sprintf("%s/%s.tsv.gz",io$data.dir,n), select=c(1,2,3,5,6)) %>%
  setnames(c("id_met","id","anno","N","rate")) %>% .[N>=opts$min.CpGs] %>% .[,N:=NULL]
) %>% rbindlist

################################################
## Merge methylation data and sample metadata ##
################################################

met_dt <- met_dt %>% merge(sample_metadata, by="id_met") %>% droplevels()

#######################################
## Calculate M value from Beta value ##
#######################################

met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#################
## Filter data ##
#################

# Filter features by coverage
nsamples <- length(unique(met_dt$id_met))
met_dt <- met_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

############################
## Regress out covariates ##
############################

# Global methylation rate
foo <- met_dt[,.(mean=mean(m)),by=c("id_met","anno")]
met_dt <- met_dt %>% merge(foo, by=c("id_met","anno")) %>%
  .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")]

# Filter features by variance
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n = opts$nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

#######################################
## Create matrix from the data.table ##
#######################################

dmatrix_list <- met_dt %>% split(.$anno) %>% 
  map(~
    dcast(.[,c("id_met","id","m")], formula=id_met~id, value.var="m") %>% as.data.frame %>% 
      tibble::column_to_rownames("id_met") %>% as.matrix %>% t
     )
lapply(dmatrix_list,dim)

####################
## Fit MOFA model ##
####################

# Create MOFAobject
MOFAobject <- createMOFAobject(dmatrix_list)

# Set options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 2

TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 1000
TrainOptions$tolerance <- 0.25
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

