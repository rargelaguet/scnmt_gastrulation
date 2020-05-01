library(MOFA)

#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/met/dimensionality_reduction/within_stages/load_settings.R")

###############################
## Load DNA methylation data ##
###############################

met_dt <- lapply(names(opts$annos), function(n)
  fread(sprintf("%s/%s.tsv.gz",io$data.dir,n), select=c(1,2,3,5,6)) %>%
  setnames(c("id_met","id","anno","N","rate")) %>% .[N>=opts$min.CpGs] %>% .[,N:=NULL]
) %>% rbindlist %>% merge(sample_metadata, by="id_met")

################
## Parse data ##
################

# Calculate M value from Beta value 
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Filter features by coverage
nsamples <- length(unique(met_dt$id_met))
met_dt <- met_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

############################
## Regress out covariates ##
############################

# (Optional) Regress out global methylation rate (per feature and stage_lineage)
# foo <- met_dt[,.(mean=mean(m)),by=c("id_met","anno")]
# met_dt <- met_dt %>% merge(foo, by=c("id_met","anno")) %>%
#   .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")]

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

# NOTE THAT MOFA V1 IS NOW OUTDATED AND YOU SHOULD USE MOFA+ 
# (https://github.com/bioFAM/MOFA2)

# Create MOFAobject
MOFAobject <- createMOFAobject(dmatrix_list)

# Set options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 2

TrainOptions <- getDefaultTrainOptions()
TrainOptions$seed <- 42

# Prepare
MOFAobject <- prepareMOFA(MOFAobject,
  ModelOptions = ModelOptions, 
  TrainOptions = TrainOptions
)

# Train the model
model <- runMOFA(MOFAobject, io$outfile)

