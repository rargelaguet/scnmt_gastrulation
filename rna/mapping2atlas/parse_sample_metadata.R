

foo <- mapping$mapping %>% as.data.table %>%
  .[,c("cell","celltype.mapped","stage.mapped")] %>%
  setnames("cell","id_rna")

bar <- fread(paste0(path2scNMT, "/sample_metadata.txt")) 

foobar <- merge(foo,bar, by="id_rna", all.y=T)  %>%
  .[is.na(celltype.mapped),c("celltype.mapped","stage.mapped"):=list(lineage,stage)]  %>%
  .[pass_rnaQC==F,c("lineage","celltype.mapped","stage.mapped"):=NA] %>%
  .[,lineage10x:=celltype.mapped] %>% .[,celltype.mapped:=NULL] %>%
  .[is.na(lineage10x),lineage10x:=lineage] %>% .[,lineage:=NULL] %>%
  .[,lineage10x:=stringr::str_replace_all(lineage10x," ","_")]

fwrite(foobar, file="/Users/ricard/data/gastrulation/sample_metadata2.txt", sep="\t", row.names=F, col.names=T, na="NA", quote=F)





library(data.table)
library(purrr)

sample_metadata <- fread("/Users/ricard/data/gastrulation/sample_metadata.txt") 

# Create lineage10x_2 by merging similar lineages
sample_metadata %>%
  .[,lineage10x_2:=lineage10x] %>%
  
  # Mature mesoderm
  .[lineage10x%in%c("Pharyngeal_mesoderm","Paraxial_mesoderm","ExE_mesoderm","Mesenchyme","Mixed_mesoderm","Blood_progenitors_2","Haematoendothelial_progenitors"), lineage10x_2:="Mature_mesoderm"] %>%
  .[lineage10x%in%c("Intermediate_mesoderm","Somitic_mesoderm","Caudal_mesoderm"), lineage10x_2:="Nascent_mesoderm"] %>%
  # Embryonic endoderm lineages
  .[lineage10x%in%c("Gut","Def._endoderm","Notochord"), lineage10x_2:="Embryonic_endoderm"] %>%
  # Extra-embryonic endoderm lineages
  .[lineage10x%in%c("Parietal_endoderm","ExE_endoderm"), lineage10x_2:="Visceral_endoderm"] %>%
  # Primitive streak lineages
  .[lineage10x%in%c("Caudal_epiblast","Anterior_Primitive_Streak"), lineage10x_2:="Primitive_Streak"] %>%
  # Ectoderm lineages
  .[lineage10x%in%c("Rostral_neurectoderm","Surface_ectoderm"), lineage10x_2:="Ectoderm"] %>%
  # PGC are not expected
  .[lineage10x%in%c("PGC"), lineage10x_2:=NA] %>%
  # At E6.5 we do not expect any mature ectoderm
  .[stage=="E6.5" & lineage10x%in%c("Rostral_neurectoderm","Surface_ectoderm"), lineage10x_2:="Epiblast"]

# Add sex information
sample_metadata <- sample_metadata %>%
  merge(fread("/Users/ricard/data/gastrulation/sample_metadata_sex.txt") %>% .[,c("id_rna","sex","pass_sexQC")], by=c("id_rna"))

## START TEST ##
sample_metadata %>%
  .[lineage10x_2%in%c("Mature_mesoderm","Nascent_mesoderm"), lineage10x_2:="Mesoderm"] %>%
  .[stage=="E7.5" & lineage10x_2%in%c("Embryonic_endoderm","Visceral_endoderm"), lineage10x_2:="Endoderm"]  %>%
  .[stage=="E5.5" & lineage10x_2%in%c("Endoderm"), lineage10x_2:="Visceral_endoderm"] %>%
  .[stage=="E6.5" & lineage10x_2%in%c("Endoderm"), lineage10x_2:="Visceral_endoderm"] 
## END TEST ##

table(sample_metadata[pass_metQC==T & stage=="E7.5",lineage10x_2])

fwrite(sample_metadata, "/Users/ricard/data/gastrulation/sample_metadata.txt", sep="\t", col.names=T, row.names=F, na="NA", quote=F)
