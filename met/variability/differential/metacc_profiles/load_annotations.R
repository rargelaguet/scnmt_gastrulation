
# Load genomic annotations
anno_list <- list()
for (anno in names(opts$annos)) {
  tmp <- fread(
    file = sprintf("%s/%s.bed.gz",io$features.dir,anno),
    colClasses = c("V1"="factor", "V2"="integer", "V3"="integer", "V4"="factor", "V5"="factor", "V6"="factor")
  ) %>% setnames(c("chr","start","end","strand","id","anno"))
  
  # Define central position for the window approach
  stopifnot(all(tmp[,end] > tmp[,start]))
  anno_list[[anno]] <- tmp %>%
    .[,c("chr","start","end","strand","id","anno")] %>%
    .[,center:=round(end+start)/2] %>%
    .[,c("start","end"):=NULL] %>%
    .[, c("start","end") := list(center-opts$window_size,center+opts$window_size)]
}

anno_df <- rbindlist(anno_list) %>% 
  .[,chr:=as.factor(sub("chr","",chr))]

integer.cols <- c("start","end","center")
anno_df %>%  .[,(integer.cols):=lapply(.SD, as.integer),.SDcols=(integer.cols)]
  
anno_df %>% setkey(chr,start,end)

rm(anno_list)
