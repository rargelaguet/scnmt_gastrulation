# Sanity checks
# (...)

###########################
## Load methylation data ##
###########################

print("Loading methylation...")

met_list <- list()
for (i in opts$met.cells) {
  print(i)
  # met_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), select = c(1,2,3,4), colClasses = c("chr"="factor", "start"="integer", "end"="integer", "rate"="numeric")) %>% 
  met_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), colClasses = c("chr"="character", "pos"="integer", "rate"="numeric")) %>% 
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    .[,id_met:=as.factor(i)] %>% 
    .[,c("start","end"):=pos] %>%
    setkey("chr","start","end") %>%
    .[,bp:=start] %>%
    foverlaps(.,anno_df, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=args$met_tile*round(dist/args$met_tile)] %>%
    .[,list(rate=100*mean(rate), N=.N),by=.(id_met,id,dist,anno)]
}
met.dt <- rbindlist(met_list) %>%
  .[,c("id_met","id","context"):=list(as.factor(id_met),as.factor(id),as.factor("CG"))]
  
rm(met_list)

print(object.size(met.dt), units="auto")


#############################
## Load accessibility data ##
#############################

print("Loading accessibility...")

acc_list <- list()
for (i in opts$acc.cells) {
  print(i)
  # acc_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), select = c(1,2,3,4), colClasses = c("chr"="factor", "start"="integer", "end"="integer", "rate"="numeric")) %>% 
  acc_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), colClasses = c("chr"="character", "pos"="integer", "rate"="numeric")) %>% 
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    .[,id_acc:=as.factor(i)] %>% 
    .[,c("start","end"):=pos] %>%
    setkey("chr","start","end") %>%
    .[,bp:=start] %>%
    foverlaps(.,anno_df, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=args$acc_tile*round(dist/args$acc_tile)] %>%
    .[,list(rate=100*mean(rate), N=.N),by=.(id_acc,id,dist,anno)]
}
acc.dt <- rbindlist(acc_list) %>%
  .[,c("id_acc","id","context"):=list(as.factor(id_acc),as.factor(id),as.factor("GC"))]
  
rm(acc_list)

print(object.size(acc.dt), units="auto")
