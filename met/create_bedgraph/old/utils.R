######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_cpgs.x), met_cpgs.x := 0L] %>%
    .[is.na(met_cpgs.y), met_cpgs.y := 0L] %>%
    .[is.na(nonmet_cpgs.x), nonmet_cpgs.x := 0L] %>%
    .[is.na(nonmet_cpgs.y), nonmet_cpgs.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_cpgs=met_cpgs.x+met_cpgs.y, nonmet_cpgs=nonmet_cpgs.x+nonmet_cpgs.y)]
}

fread_and_merge <- function(dt, file){
  fread(file, colClasses=list(factor=1L), select=c(1,2,5)) %>% 
    setnames(c("chr","pos","rate")) %>%
    merge_and_sum(dt)
}