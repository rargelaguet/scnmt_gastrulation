
#############################
## Load metadata enhancers ##
#############################

# Load genomic features from bed file
feature_metadata <- lapply(names(opts$annos), function(n) fread(sprintf("%s/%s.bed",io$features.dir,n), showProgress=F)) %>% rbindlist  
colnames(feature_metadata) <- c("chr", "start", "end", "strand", "id", "anno")


########################
## Load enhancer data ##
########################

enh <- fread(io$enh, stringsAsFactors=T, showProgress=F)
colnames(enh) <- c("probe","chr","start","end","strand","feature","id","description","feature_strand","type","feature_location","distance","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")

enh <- enh[,c("chr","start","end","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")] %>%   
  .[,start:=start-1] %>%
  .[,length:=end-start] %>%
  .[,range:=paste(start, end, sep="_")] 

enh <- merge(enh, feature_metadata[,c("start","end","chr","id","anno")], by=c("start","end","chr")) %>%
  .[!duplicated(.$range), ] %>%
  #.[id %in% met.diff.enh | id %in% acc.diff.enh] %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)] #%>%
#.[anno=="Ectoderm enhancers"]

################
## Parse data ##
################

# log transform
list.enh_log <- NULL
for (i in 1:length(enh)) {
  foo <- as.matrix(enh[,..i])
  if(is.numeric(foo)){
    median <- matrixStats::colMedians(foo)
    #print(median)
    if(median<100){
      foo <- as.data.table(log2(foo))
    }
    else{
      foo <- as.data.table(foo)
    }
  }
  else{
    foo <- as.data.table(foo)
    #print("not numeric")
  }
  list.enh_log[[colnames(enh[,..i])]] <- foo
}
enh_log <- as.data.table(list.enh_log)
colnames(enh_log) <- c("start","end","chr","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine","length","range","id","anno")

# replace Inf and NA values with 0 in data table
enh_log <- do.call(data.table, lapply(enh_log, function(x) {
  replace(x, is.infinite(x) | is.na(x), 0)
})
)


###########################################
## Identify cell type specific enhancers ##
###########################################

# ESC enhancers
enh.esc <- fread(io$enh.esc, stringsAsFactors=T, showProgress=F)
colnames(enh.esc) <- c("probe","chr","start","end","strand","feature","id","description","feature_strand","type","feature_location","distance","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")

enh.esc <- enh.esc[,c("chr","start","end","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")] %>%   
  .[,start:=start-1] %>%
  .[,length:=end-start] %>%
  .[,range:=paste(start, end, sep="_")] 

enh.esc <- merge(enh.esc, feature_metadata[,c("start","end","chr","id","anno")], by=c("start","end","chr")) %>%
  .[!duplicated(.$range), ] %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)] %>%
  .[,id]

# brain enhancers
enh.brain <- fread(io$enh.brain, stringsAsFactors=T, showProgress=F)
colnames(enh.brain) <- c("probe","chr","start","end","strand","feature","id","description","feature_strand","type","feature_location","distance","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")

enh.brain <- enh.brain[,c("chr","start","end","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")] %>%   
  .[,start:=start-1] %>%
  .[,length:=end-start] %>%
  .[,range:=paste(start, end, sep="_")] 

enh.brain <- merge(enh.brain, feature_metadata[,c("start","end","chr","id","anno")], by=c("start","end","chr")) %>%
  .[!duplicated(.$range), ] %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)] %>%
  .[,id]

# gut enhancers
enh.gut <- fread(io$enh.gut, stringsAsFactors=T, showProgress=F)
colnames(enh.gut) <- c("probe","chr","start","end","strand","feature","id","description","feature_strand","type","feature_location","distance","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")

enh.gut <- enh.gut[,c("chr","start","end","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")] %>%   
  .[,start:=start-1] %>%
  .[,length:=end-start] %>%
  .[,range:=paste(start, end, sep="_")] 

enh.gut <- merge(enh.gut, feature_metadata[,c("start","end","chr","id","anno")], by=c("start","end","chr")) %>%
  .[!duplicated(.$range), ] %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)] %>%
  .[,id]

# heart enhancers
enh.heart <- fread(io$enh.heart, stringsAsFactors=T, showProgress=F)
colnames(enh.heart) <- c("probe","chr","start","end","strand","feature","id","description","feature_strand","type","feature_location","distance","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")

enh.heart <- enh.heart[,c("chr","start","end","serum","E10.5_midbrain","E10.5_heart","E12.5_intestine")] %>%   
  .[,start:=start-1] %>%
  .[,length:=end-start] %>%
  .[,range:=paste(start, end, sep="_")] 

enh.heart <- merge(enh.heart, feature_metadata[,c("start","end","chr","id","anno")], by=c("start","end","chr")) %>%
  .[!duplicated(.$range), ] %>%
  .[,anno:=stringr::str_replace_all(anno,opts$annos)] %>%
  .[,id]


########################################################
## Select enhancers that show tissue-specific marking ##
########################################################

enh_marked_log_all <- enh_log %>%
  .[,marked_esc:=id %in% enh.esc] %>%
  .[,marked_brain:=id %in% enh.brain] %>%
  .[,marked_gut:=id %in% enh.gut] %>%
  .[,marked_heart:=id %in% enh.heart]

enh_marked_log_all$nr_marked <- enh_marked_log_all$marked_esc + enh_marked_log_all$marked_brain + enh_marked_log_all$marked_gut + enh_marked_log_all$marked_heart 

enh_marked <- merge(enh, enh_marked_log_all[,c("id","nr_marked","marked_esc","marked_brain","marked_gut","marked_heart")], by="id") %>%
  .[nr_marked>=1 & nr_marked<3]

enh_marked_log <- enh_marked_log_all[nr_marked>=1 & nr_marked<3]


#############################################################################
## Identify enhancers that are differentially marked between ESC and brain ##
#############################################################################

enh_ect <-  enh_marked_log_all %>%
  .[anno=="Ectoderm enhancers"] %>%
  .[nr_marked<3] %>%
  #.[serum!=0] %>%
  .[marked_esc==T | marked_brain==T] %>%
  .[,ratio:=E10.5_midbrain/serum]

enh_brain <- dplyr::top_n(enh_ect, opts$nr.diff, ratio) 
enh_brain$class <- "brain"
enh_esc <- dplyr::top_n(enh_ect, -opts$nr.diff, ratio)
enh_esc$class <- "esc"
enh_diff <- as.data.table(rbind(enh_brain, enh_esc))

