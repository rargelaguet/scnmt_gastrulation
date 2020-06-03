source("/Users/ricard/scnmt_gastrulation/settings.R")

foo <- fread("/Users/ricard/data/gastrulation/rna/results/mapping/mapping.txt.gz") %>%
  .[,c("id_rna","celltype.mapped","celltype.score")] %>%
  setnames(c("id_rna","celltype_v2","score_v2"))

bar <- sample_metadata[pass_rnaQC==T & stage%in%c("E6.5","E7.5")] %>%
  .[,c("id_rna","lineage10x")] %>% setnames(c("id_rna","celltype_v1"))

foobar <- merge(foo,bar, by="id_rna")

baz <- foobar[celltype_v1!=celltype_v2]
