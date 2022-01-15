
###########################
## Load methylation data ##
###########################

cat("Loading DNA methylation data...")

met.dt <- opts$stage_lineage %>% map(function(i) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_pseudobulk_raw,i), showProgress=F, select=c(1,2,5)) %>%
    .[!is.na(pos)] %>%
    .[,c("chr","pos","rate"):=list(as.factor(chr),as.integer(pos),as.integer(rate))] %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    # .[,id:=as.character(id)] %>%
    .[,dist:=as.integer(ifelse(strand %in% c("+","*"),bp-center,center-bp))] %>% 
    .[,dist:=opts$met.tile*round(dist/opts$met.tile)] %>%
    .[,list(rate=as.integer(mean(rate)), n=.N),by=c("id","dist","anno")] %>%
    .[,stage_lineage:=as.factor(i)]
# }) %>% rbindlist %>% .[!is.na(bp)] %>%
}) %>% rbindlist %>%
  .[,c("context"):=list(as.factor("CG"))]
  
#############################
## Load accessibility data ##
#############################

cat("Loading chromatin accessibility data...")

acc.dt <- opts$stage_lineage %>% map(function(i) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_pseudobulk_raw,i), showProgress=F, select=c(1,2,5)) %>%
    .[!is.na(pos)] %>%
    .[,c("chr","pos","rate"):=list(as.factor(chr),as.integer(pos),as.integer(rate))] %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    # .[,id:=as.character(id)] %>%
    .[,dist:=as.integer(ifelse(strand %in% c("+","*"),bp-center,center-bp))] %>% 
    .[,dist:=opts$acc.tile*round(dist/opts$acc.tile)] %>%
    .[,list(rate=as.integer(mean(rate)), n=.N),by=c("id","dist","anno")] %>%
    .[,stage_lineage:=as.factor(i)]
}) %>% rbindlist %>%
  .[,c("context"):=list(as.factor("GC"))]
