# Function to load multipke methylation files and overlap them with a given genomic annotation
#   indir: 
#   cells: 
#   annotation: a data.table or data.frame with columns (chr,start,end,strand,id)
# loadMultipleMetData <- function(indir, cells, annotation=NULL) {
#   data <- cells %>% map(function(cell) { 
#     fread(sprintf("zcat < %s/%s.tsv.gz",indir,cell), sep="\t", showProgress=FALSE, stringsAsFactors=TRUE) %>% 
#       .[,id_met:=factor(cell, levels=levels(cells))] %>% .[,c("start","end"):=list(pos,pos)] %>% 
#       setkey("chr","start","end") %>%
#       foverlaps(annotation, nomatch=0) %>% 
#       .[,dist:=ifelse(strand %in% c("+","*"),pos-TSS,TSS-pos)] %>% 
#       .[,c("id_met","ens_id","dist","rate","type")]
#   }) %>% rbindlist
#   return(data)
# }
loadMultipleMetData <- function(indir, cells, annotation=NULL) {
  data_list <- list()
  for (cell in cells) {
    if (file.exists(sprintf("%s/%s.tsv.gz",indir,cell))) {
      # print(sprintf("Loading %s...",cell))
      data_list[[cell]] <- fread(sprintf("zcat < %s/%s.tsv.gz",indir,cell), sep="\t", stringsAsFactors=T, showProgress=F) %>% 
        .[,id_met:=cell] %>% .[,c("start","end"):=list(pos,pos)] %>% 
        setkey("chr","start","end") %>%
        foverlaps(annotation, nomatch=0) %>% 
        .[,dist:=ifelse(strand %in% c("+","*"),pos-TSS,TSS-pos)] %>% 
        .[,c("id_met","ens_id","dist","rate","type")]
    }
  }
  gc(reset=TRUE)
  return(rbindlist(data_list))
}

theme_pub <- function() {
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(size=rel(1.3),colour="black"),
    axis.text.y = element_text(size=rel(1.3), colour="black"),
    axis.title.x = element_text(size=rel(1.5)),
    axis.title.y = element_text(size=rel(1.5)),
    axis.line = element_line(size=rel(1.0)),
    axis.ticks.x = element_line(size=rel(1.1), color="black"),
    axis.ticks.y = element_line(size=rel(1.1), color="black"),
    legend.key = element_blank(),
    legend.position = "top",
    # legend.position = c(0.5,1.0),
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.key.width = unit(2.0,"line"),
    legend.key.height = unit(2.0,"line"),
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.3)),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
}

mean_sd <- function(x) { return(data.frame(y=mean(x), ymin=mean(x)-sd(x)/2, ymax=mean(x)+sd(x)/2)) }