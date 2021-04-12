suppressMessages(library(argparse))

#####################
## Argument parser ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--anno',         type="character",                         help='Genomic annotation')
p$add_argument('--motif_set',    type="character",   default="jaspar2020", help='Motif annotation')
p$add_argument('--width',        type="integer",     default=7,            help='Motif width (see motifmatchr::matchMotifs)')
p$add_argument('--extend',       type="integer",    default=50,           help='Extend window around the motif (bp)')
p$add_argument('--cutOff',       type="double",      default=5e-05,        help='Motif matching score cutoff (see motifmatchr::matchMotifs)')
p$add_argument('--outprefix',      type="character",                         help='Output file')
args <- p$parse_args(commandArgs(TRUE))


## START TEST ##
# args$anno <- "prom_2000_2000"
# args$motif_set <- "jaspar2020"
# args$width <- 7
# args$cutOff <- 5e-05
# args$outprefix <- sprintf("%s/features/motifs/%s_%s",io$basedir,args$anno,args$motif_set)
## END TEST ##

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/scnmt_gastrulation/motifs/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/scnmt_gastrulation/motifs/utils.R")
} else {
  stop("Computer not recognised")
}


# I/O

# Options

##############################
## Load genomic annotations ##
##############################

granges <- sprintf("%s/%s.bed.gz",io$features.dir,args$anno) %>% 
  load_annotation_as_granges_list %>% .[[1]]

############################
## Load motif annotations ##
############################

pwms <- load_motif_annotation(
  motifSet = args$motif_set,
  species = "Homo sapiens"
)[["motifs"]]

#################################
## Create motif overlap matrix ##
#################################

motifOverlap.se <- create_motif_overlap_matrix(
  pwms, 
  granges, 
  genome = "BSgenome.Mmusculus.UCSC.mm10", 
  cutOff = args$cutOff, 
  width = args$width
)

# Rename motifs
names(motifOverlap.se$motifPositions) <- names(pwms)

########################################################
## Create BED file with the motif matches coordinates ##
########################################################

names(motifOverlap.se$motifPositions) %>% walk(function(i) {
  motif_matches.dt <- as.data.table(motifOverlap.se$motifPositions[[i]]) %>% 
    .[,c(1,2,3,5)] %>% 
    setnames(c("chr","start","end","strand")) %>%
    .[,anno:=sprintf("%s_%s",args$anno,i)] %>%
    .[,chr:=gsub("chr","",chr)] %>%
    .[,c("start","end"):=list(start-args$extend,end+args$extend)] %>%
    .[,id:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end","strand","id","anno")] %>% 
    setorder(chr,start,end)
  
    fwrite(motif_matches.dt, sprintf("%s_%s.bed.gz",args$outprefix,i), sep="\t", quote=F, col.names = F)
})# %>% rbindlist %>% setorder(chr,start,end)


##########
## Save ##
##########

saveRDS(motifOverlap.se, paste0(args$outprefix,".rds"))
