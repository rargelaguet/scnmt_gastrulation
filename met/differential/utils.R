round_df <- function(df, digits) {
  nums <- names(which(vapply(df, is.numeric, FUN.VALUE = logical(1))))
  df[,(nums) := round(.SD,digits), .SDcols=nums]
  return(df)
}

scatter_theme <- function(){
  p <- theme(
      plot.title=element_text(size=28, face='bold', margin=margin(0,0,10,0), hjust=0.5),
      axis.text=element_text(size=rel(1.75), color='black'),
      axis.title=element_text(size=rel(1.95), color='black'),
      axis.title.y = element_text(margin=margin(0,10,0,0)),
      axis.title.x = element_text(margin=margin(10,0,0,0)),
      legend.position="none",
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}

gg_volcano_plot <- function(tmp, title = "", label=10, xlim=NA, ylim=NA) {
  negative_hits <- length(tmp[sig==TRUE & diff<0,id])
  positive_hits <- length(tmp[sig==TRUE & diff>0,id])
  all <- nrow(tmp)
  
  if (is.na(xlim)) xlim <- max(abs(tmp$diff), na.rm=T)
  if (is.na(ylim)) ylim <- max(-log10(tmp$p.value), na.rm=T)
  
  p <- ggplot(tmp, aes(x=diff, y=-log10(p.value))) +
    labs(title=title, x="Methylation rate difference", y=expression(paste("-log"[10],"(p.value)"))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=max(-log10(p.value)-1, na.rm=T)), color="orange") +
    geom_point(aes(color=sig), size=2) +
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-xlim-5,xlim+5)) +
    scale_y_continuous(limits=c(0,ylim+1)) +
    annotate("text", x=0, y=ylim+1, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-xlim+10, y=ylim+1, size=7, label=sprintf("%d (-)",negative_hits)) +
    annotate("text", x=xlim-10, y=ylim+1, size=7, label=sprintf("%d (+)",positive_hits)) +
    scatter_theme()
  
  if (label>0) {
    p <- p + ggrepel::geom_text_repel(data=head(tmp[sig==TRUE],n=label), aes(x=diff, y=-log10(p.value), label=gene), size=6)
  }
  return(p)
} 

boxplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,20,0)),
    axis.title.y = element_text(colour="black", size=20, vjust=1.5),
    axis.title.x = element_text(colour="black", size=20, vjust=1.5, margin=margin(15,0,0,0)),
    axis.text.x = element_text(colour="black",size=rel(1.6)),
    axis.text.y = element_text(colour="black",size=rel(1.6)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(0.8)),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position="none",
    legend.text=element_text(size=15),
    legend.title=element_blank(),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
}

round_df <- function(df, digits) {
  nums <- names(which(vapply(df, is.numeric, FUN.VALUE = logical(1))))
  df[,(nums) := round(.SD,digits), .SDcols=nums]
  return(df)
}
