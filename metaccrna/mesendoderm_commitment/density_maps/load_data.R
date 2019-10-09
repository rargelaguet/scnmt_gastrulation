
###############################################
## Load differential DNA methylation results ##
###############################################

# Load precomputed differential results
diff.met <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$diff.met,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]


######################################################
## Load differential chromatin accessiblity results ##
######################################################

diff.acc <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j)
    fread(sprintf("%s/%s_%s.txt.gz",io$diff.acc,i,j))
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist %>% .[,sig:=padj_fdr<opts$min.fdr & abs(diff)>opts$min.diff]




## Merge methylation and accessibility results ##
diff.met <- diff.met %>% 
  .[,anno:=stringr::str_replace_all(anno,opts$met.annos)] %>%
  .[,c("id","anno","diff","sig","comparison")]
  
diff.acc <- diff.acc %>% 
  .[,anno:=stringr::str_replace_all(anno,opts$acc.annos)] %>%
  .[,c("id","anno","diff","sig","comparison")]


diff.metacc <- rbind(
  diff.met[,type:="met"], 
  diff.acc[,type:="acc"]
) %>% dcast(id+gene+lineage+anno~type, value.var=c("diff","sig")) %>% .[complete.cases(.)]

