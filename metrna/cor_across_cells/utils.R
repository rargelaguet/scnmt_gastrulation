
# Define ggplot2 theme for boxplots
boxplot_theme <- function() {
  p <- theme_classic() + theme(
    axis.title = element_text(colour="black", size=15, vjust=1.5),
    axis.text = element_text(colour="black",size=rel(1.2)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks = element_line(colour="black", size=rel(0.8)),
    legend.position="top",
    legend.direction = 'horizontal',
    legend.text=element_text(size=17),
    legend.title=element_blank()
  )
}

gg_volcano_plot <- function(to.plot, top_hits=10, xlim=NULL, ylim=NULL){
  
  # Compute positive and negative hits in terms of correlation
  negative_hits <- to.plot[sig==TRUE & r<0,id]
  positive_hits <- to.plot[sig==TRUE & r>0,id]
  all <- nrow(to.plot)
  
  if (is.null(xlim))
    xlim <- max(abs(to.plot$r), na.rm=T)
  if (is.null(ylim))
    ylim <- max(to.plot$log_padj_fdr, na.rm=T)
  
  max.yaxis <- max(to.plot$log_padj_fdr)
  
  p <- ggplot(to.plot, aes(x=r, y=log_padj_fdr)) +
    labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
    #geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim), color="orange", size=0.25) +
    # geom_point(aes(color=sig, size=sig)) +
    ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.75,1.25)) +
    scale_x_continuous(limits=c(-xlim-0.15,xlim+0.15)) +
    scale_y_continuous(limits=c(0,ylim+6)) +
    annotate("text", x=0, y=ylim+5, size=4, label=sprintf("(%d)", all)) +
    annotate("text", x=-xlim-0.1, y=ylim+5, size=4, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=xlim+0.1, y=ylim+5, size=4, label=sprintf("%d (+)",length(positive_hits))) +
    theme_classic() + 
    theme(
      axis.text=element_text(size=rel(0.75), color='black'),
      axis.title=element_text(size=rel(1), color='black'),
      legend.position="none"
    )
  
    if (top_hits>0) { 
      foo <- to.plot[sig == TRUE] %>% setkey(padj_fdr) %>% head(n=top_hits)      
      p <- p + ggrepel::geom_text_repel(data=foo, aes(x=r, y=log_padj_fdr, label=gene), size=4, color="red")
    }
  
  return(p)
} 

gg_qqplot = function(cor_res, perm_xs, ci=0.95, title = "Quantile-quantile plot of p-values"){
  xs <- cor_res$p
  cor_res <- cor_res[, expected := -log10(1:.N / .N)]
  N  <- length(xs)
  df <- data.frame(
    observed = -log10(sort(xs)),
    permuted = -log10(sort(perm_xs)),
    expected = -log10(1:N / N),
    cupper   = -log10(qbeta(ci,     1:N, N - 1:N + 1)),
    clower   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  )
  
  ggplot(df) +
    # geom_point(aes(expected, permuted), shape=3, size=1.5, color = "cornflowerblue") +
    ggrastr::geom_point_rast(aes(expected, permuted), shape=3, size=1.5, color = "cornflowerblue") +
    geom_abline(intercept=0, slope=1, alpha=0.5, color = "darkgrey") +
    geom_line(aes(expected, cupper), linetype=2, color = "darkgrey") +
    geom_line(aes(expected, clower), linetype=2, color = "darkgrey") +
    #geom_point(aes(expected, observed), size=2) +
    geom_point(data=cor_res, aes(expected, -log10(sort(p)), color = sig), size=2) +
    scale_color_manual(values=c("black","red")) +
    labs(x=expression(paste("Expected -log"[10],"(p.value)")), y=expression(paste("Observed -log"[10],"(p.value)"))) +
    theme_classic() + 
    theme(
      axis.text=element_text(size=rel(1), color='black'),
      axis.title=element_text(size=rel(1), color='black'),
      legend.position="none"
    )
}

round_df <- function(df, digits) {
  nums <- names(which(vapply(df, is.numeric, FUN.VALUE = logical(1))))
  df[,(nums) := round(.SD,digits), .SDcols=nums]
  return(df)
}