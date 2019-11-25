
# Load genomic annotations
anno_list <- list()
for (anno in names(opts$annos)) {
  tmp <- fread(sprintf("%s/%s.bed",io$annos_dir,anno))[,c(1,2,3,4,5,6)]
  colnames(tmp) <- c("chr","start","end","strand","id","anno")
  
  # Define central position for the window approach
  if (opts$positions[anno] == "start") {
    tmp <- rbind(tmp[strand=="+",.(chr,start,strand,id,anno)] %>% .[,center:=start] %>% .[,c("start"):=NULL], 
                 tmp[strand=="-",.(chr,end,strand,id,anno)] %>% .[,center:=end] %>% .[,c("end"):=NULL]) 
  }
  if (opts$positions[anno] == "center") {
    stopifnot(all(tmp[,end] > tmp[,start]))
    tmp <- tmp[,.(chr,start,end,strand,id,anno)][,center:=round(end+start)/2][,c("start","end"):=NULL]
  }
  if (opts$positions[anno] == "end") {
    tmp <- rbind(tmp[strand=="+",.(chr,end,strand,id,anno)][,center:=end][,c("end"):=NULL], 
                 tmp[strand=="-",.(chr,start,strand,id,anno)][,center:=start][,c("start"):=NULL])
  }
  anno_list[[anno]] <- tmp %>% .[, c("start","end") := list(center-opts$window_size,center+opts$window_size)]
}

anno_df <- rbindlist(anno_list) %>% 
  .[,chr:=sub("chr","",chr)] %>%
  setkey(chr,start,end)

rm(anno_list)
