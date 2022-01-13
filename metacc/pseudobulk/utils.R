######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_sites.x), met_sites.x := 0L] %>%
    .[is.na(met_sites.y), met_sites.y := 0L] %>%
    .[is.na(nonmet_sites.x), nonmet_sites.x := 0L] %>%
    .[is.na(nonmet_sites.y), nonmet_sites.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_sites=met_sites.x+met_sites.y, nonmet_sites=nonmet_sites.x+nonmet_sites.y)]
}

fread_and_merge <- function(dt, file){
    # tmp <- fread(file, select=c(1,2,4)) %>% setnames(c("chr","pos","rate"))
    tmp <- fread(file) %>% .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))]
    stopifnot(tmp$rate%in%c(0,1))
    tmp %>%
        .[,met_sites:=rate] %>%
        .[,nonmet_sites:=+(!rate)] %>%
        .[,rate:=NULL] %>%
    merge_and_sum(dt)
}