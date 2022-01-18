
# Sanity checks
# (...)

# Load genomic annotations
anno_list <- list()
for (i in args$anno) {
  tmp <- fread(sprintf("%s/%s.bed.gz",io$features.dir,i), 
               select=1:5, colClasses=c("V1"="character", "V2"="integer", "V3"="integer", "V4"="factor", "V5"="character")) %>%
    setnames(c("chr","start","end","strand","id")) %>%
    .[,anno:=factor(i)] %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    .[,chr:=factor(chr,levels=opts$chr)]
  
  # Define central position for the window approach
  if (opts$positions[i] == "start") {
    tmp <- rbind(tmp[strand=="+",.(chr,start,strand,id,anno)] %>% .[,center:=start] %>% .[,c("start"):=NULL], 
                 tmp[strand=="-",.(chr,end,strand,id,anno)] %>% .[,center:=end] %>% .[,c("end"):=NULL]) 
  }
  if (opts$positions[i] == "center") {
    stopifnot(all(tmp[,end] > tmp[,start]))
    tmp <- tmp[,.(chr,start,end,strand,id,anno)][,center:=round(end+start)/2][,c("start","end"):=NULL]
  }
  if (opts$positions[i] == "end") {
    tmp <- rbind(tmp[strand=="+",.(chr,end,strand,id,anno)][,center:=end][,c("end"):=NULL], 
                 tmp[strand=="-",.(chr,start,strand,id,anno)][,center:=start][,c("start"):=NULL])
  }
  anno_list[[i]] <- tmp %>% .[, c("start","end") := list(center-args$window_size,center+args$window_size)]
}

# Concatenate
anno_df <- rbindlist(anno_list) %>% 
  # .[,chr:=as.factor(sub("chr","",chr))] %>%
  setkey(chr,start,end)

rm(anno_list,tmp)

