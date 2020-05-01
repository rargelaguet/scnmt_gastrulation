library(MOFA)

source("/Users/ricard/scnmt_gastrulation/acc/dimensionality_reduction/within_stages/load_settings.R")

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(names(opts$annos), function(n) 
  fread(sprintf("%s/%s.tsv.gz",io$data.dir,n), select=c(1,2,3,5,6)) %>% 
    setnames(c("id_acc","id","anno","N","rate")) %>% .[N>=opts$min.GpCs] %>% .[,N:=NULL]
  ) %>% rbindlist %>% merge(sample_metadata, by="id_acc")

################
## Parse data ##
################

# Calculate M value from Beta value
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Filter features by coverage
nsamples <- length(unique(acc_dt$id_acc))
acc_dt <- acc_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

# Regress out differences in global accessibility rate (per feature and stage_lineage)
# foo <- acc_dt[,.(mean=mean(m)),by=c("id_acc","anno")]
# acc_dt <- acc_dt %>% merge(foo, by=c("id_acc","anno")) %>%
#   .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")] %>%
#   .[,mean:=NULL]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% 
  map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$nfeatures) %>% .$id)
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

# NOTE THAT MOFA V1 IS NOW OUTDATED AND YOU SHOULD USE MOFA+ 
# (https://github.com/bioFAM/MOFA2)

MOFAobject <- createMOFAobject(dmatrix_list)

# Set options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 3

TrainOptions <- getDefaultTrainOptions()
TrainOptions$seed <- 42

# Prepare
MOFAmodel <- prepareMOFA(
  MOFAobject,
  ModelOptions = ModelOptions, 
  TrainOptions = TrainOptions
)

# Train the model
model <- runMOFA(MOFAmodel, io$outfile)
