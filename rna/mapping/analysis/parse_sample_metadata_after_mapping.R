#########
## I/O ##
#########

source("/Users/ricard/scnmt_gastrulation/settings.R")

io$mapping.dir <- paste0(io$basedir,"/rna/results/mapping")
# io$metadata <- paste0(io$basedir,"/processed/sample_metadata.txt.gz")

stop()

###############
## Load data ##
###############

# Load metadata
metadata <- fread(io$metadata) %>% 
  .[,c("celltype.mapped","celltype.score","stage.mapped","closest.cell"):=NULL]

# Load mapping results
mapping.dt <- readRDS(file)$mapping %>% as.data.table %>%
        .[,c("cell","celltype.mapped","celltype.score","stage.mapped","closest.cell")] %>%
        .[,batch:=x]
# mapping.dt <- opts$batches %>% map(function(x) {
#   file <- sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x)
#   if (file.exists(file)) {
#     readRDS(file)$mapping %>% as.data.table %>%
#       .[,c("cell","celltype.mapped","celltype.score","stage.mapped","closest.cell")] %>%
#       .[,batch:=x]
#   }
# }
# ) %>% rbindlist

table(mapping.dt$celltype.mapped)
table(mapping.dt$stage.mapped)
unique(mapping.dt$batch)

#########################################
## Quick comparisons to other mappings ##
#########################################

# mapping.old <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/backups/sample_metadata_6No.txt.gz") %>%
#   .[,c("cell","batch","celltype.mapped","celltype.score")]
# foo <- merge(mapping.dt,mapping.old, by=c("cell","batch"))
# foo[celltype.mapped.x==celltype.mapped.y]

###########
## Merge ##
###########

metadata <- metadata %>% 
  merge(mapping.dt,by=c("cell","batch"),all.x=TRUE)

head(metadata)

#################
## Save output ##
#################

fwrite(metadata, io$metadata, sep="\t", na="NA", quote=F)


