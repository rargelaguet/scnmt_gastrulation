dt <- fread("/Users/ricard/data/gastrulation/met/cpg_level/pseudobulk/E4.5_Epiblast.tsv.gz") %>%
  .[chr=="1",c(1,2,5)] %>% 
  setnames("pos","start") %>%
  setnames("rate","value") %>%
  .[,end:=start] %>%
  .[,chr:=paste0("chr",chr)] %>%
  .[,c("chr","start","end","value")] %>% 
  .[value<10, value:=10]
fwrite(dt, "/Users/ricard/data/scnmt_gastrulation/met/bedgraph/test.bedgraph", col.names = F, sep="\t")
