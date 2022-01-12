# sample_metadata <- fread("")
library(data.table)
library(purrr)
plate_metadata.dt <- fread("/Users/argelagr/data/gastrulation/plate_metadata.txt")

tmp <- plate_metadata.dt[,c("plate","stage")] %>% unique
tmp %>% .[,id:=1:.N,by="stage"] %>% .[,sample:=sprintf("%s_sample_%s",stage,id)] %>% .[,c("id","stage"):=NULL]
plate_metadata.dt <- plate_metadata.dt %>% merge(tmp,by="plate")
fwrite(plate_metadata.dt, file.path("/Users/argelagr/data/gastrulation/plate_metadata.txt"), quote=F, na="NA", sep="\t")


# foo <- fread("/Users/argelagr/data/scnmt_gastrulation_argelaguet2019/processed/sample_metadata.txt.gz")
# bar <- foo %>% merge(plate_metadata.dt[,c("plate","embryo","sample","stage")],by=c("plate","embryo","stage")) %>%
#   .[,c("id_rna", "id_met", "id_acc", "sample", "plate", "embryo", "stage")]
# fwrite(bar, file.path("/Users/argelagr/data/scnmt_gastrulation_argelaguet2019/processed/sample_metadata.txt.gz"), quote=F, na="NA", sep="\t")
