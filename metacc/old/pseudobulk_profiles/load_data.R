
###########################
## Load methylation data ##
###########################

cat("Loading DNA methylation data...")

met.dt <- opts$met.cells %>% map(function(i) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), showProgress=F, select=c(1,2,5)) %>%
    # .[!is.na(pos)] %>%
    .[,c("chr","pos","rate"):=list(as.factor(chr),as.integer(pos),as.integer(rate))] %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    # .[,id:=as.character(id)] %>%
    .[,dist:=as.integer(ifelse(strand %in% c("+","*"),bp-center,center-bp))] %>% 
    .[,dist:=opts$met.tile*round(dist/opts$met.tile)] %>%
    .[,list(rate=as.integer(mean(rate)), n=.N),by=c("id","dist","anno")] %>%
    .[,cell:=as.factor(i)]
  # }) %>% rbindlist %>% .[!is.na(bp)] %>%
}) %>% rbindlist %>%
  .[,c("context"):=list(as.factor("CG"))]


# met_list <- list()
# for (cell in opts$met.cells) {
#   tmp <- fread(sprintf("%s/%s.tsv.gz",io$met.dir,cell), showProgress=F) %>%
#     .[,c("chr","pos","rate")] %>% .[,id_met:=cell] %>% 
#     .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
#     
#     foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
#     .[,id:=as.character(id)] %>%
#     .[,dist:=as.integer(ifelse(strand %in% c("+","*"),bp-center,center-bp))] %>% 
#     .[, dist:=opts$met.tile*round(dist/opts$met.tile)] %>%
#     .[,list(rate=as.integer(100*mean(rate)), n=.N),by=.(id_met,id,dist,anno)]
#   met_list[[cell]] <- tmp
# }
# met <- rbindlist(met_list) %>%
#   .[,c("id_met","id","context"):=list(as.factor(id_met),as.factor(id),as.factor("CG"))]
#   
# rm(met_list)

#############################
## Load accessibility data ##
#############################

cat("Loading chromatin accessibility data...")

acc.dt <- opts$cells %>% map(function(i) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), showProgress=F, select=c(1,2,5)) %>%
    # .[!is.na(pos)] %>%
    .[,c("chr","pos","rate"):=list(as.factor(chr),as.integer(pos),as.integer(rate))] %>%
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    # .[,id:=as.character(id)] %>%
    .[,dist:=as.integer(ifelse(strand %in% c("+","*"),bp-center,center-bp))] %>% 
    .[,dist:=opts$acc.tile*round(dist/opts$acc.tile)] %>%
    .[,list(rate=as.integer(mean(rate)), n=.N),by=c("id","dist","anno")] %>%
    .[,cell:=as.factor(i)]
}) %>% rbindlist %>%
  .[,c("context"):=list(as.factor("GC"))]

# acc_list <- list()
# for (cell in opts$acc.cells) {
#   tmp <- fread(sprintf("%s/%s.tsv.gz",io$acc.dir,cell), showProgress=F) %>%
#     .[,c("chr","pos","rate")] %>% .[,id_acc:=cell] %>% 
#     .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
#     
#     foverlaps(.,anno_df.acc, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
#     .[,id:=as.character(id)] %>%
#     .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
#     .[,dist:=as.integer(opts$acc.tile*round(dist/opts$acc.tile))] %>%
#     .[,list(rate=as.integer(100*mean(rate)), n=.N),by=.(id_acc,id,dist,anno)]
#   acc_list[[cell]] <- tmp
# }
# acc <- rbindlist(acc_list) %>%
#   .[,c("id_acc","id","context"):=list(as.factor(id_acc),as.factor(id),as.factor("GC"))]
#   
# rm(acc_list)

#####################################
## Merge data with sample metadata ##
#####################################

# met <- met %>% merge(sample_metadata, by="id_met") %>% droplevels()
# acc <- acc %>% merge(sample_metadata, by="id_acc") %>% droplevels()

