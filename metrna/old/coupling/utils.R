
# Define ggplot2 theme for boxplots

# boxplot_theme <- function() {
#   p <- theme(
#     plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,20,0)),
#     axis.title.y = element_text(colour="black", size=20, vjust=1.5),
#     axis.title.x = element_text(colour="black", size=20, vjust=1.5, margin=margin(15,0,0,0)),
#     axis.text.x = element_text(colour="black",size=rel(1.6)),
#     axis.text.y = element_text(colour="black",size=rel(1.6)),
#     axis.line = element_line(colour="black", size=rel(0.7)),
#     axis.ticks.x = element_line(colour="black", size=rel(0.8)),
#     axis.ticks.y = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     legend.position="top",
#     legend.direction = 'horizontal',
#     legend.text=element_text(size=17),
#     legend.title=element_blank(),
#     legend.background=element_blank(),
#     panel.border = element_blank()
#   )
# }

gg_volcano_plot <- function(cor_samples, title = "", label=F) {
  foo <- cor_samples[sig == TRUE] %>% setkey(padj_fdr) %>% head(n=10)
  p <- ggplot(cor_samples, aes(x=r, y=-log10(p))) +
    labs(title=title, x="Weighted Pearson correlation", y=expression(paste("-log"[10],"(",plain(p),")"))) +
    #geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=8.1), color="orange") +
    geom_point(aes(color=sig), size=2) +
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-1,1)) +
    scale_y_continuous(limits=c(0,8.5)) +
    annotate("text", x=0, y=8.47, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-0.5, y=8.47, size=7, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=0.5, y=8.47, size=7, label=sprintf("%d (+)",length(positive_hits))) +
    # geom_text(data=cor_samples[sig == TRUE], aes(x=r, y=log_padj_fdr, label=gene), vjust=-0.0, hjust=-0.3) +
    theme_classic()
  if (label==T) { 
    p <- p + ggrepel::geom_text_repel(data=foo, aes(x=r, y=-log10(p), label=gene), size=6, color="red")
  }
  return(p)
} 

gg_qqplot = function(cor_res, perm_xs, ci=0.95, title = "Quantile-quantile plot of p-values"){
  xs <- cor_res$p
  cor_res <- cor_res[, expected := -log10(1:.N / .N)]
  N  <- length(xs)
  df <- data.frame(observed = -log10(sort(xs)),
                   permuted = -log10(sort(perm_xs)),
                   expected = -log10(1:N / N),
                   cupper   = -log10(qbeta(ci,     1:N, N - 1:N + 1)),
                   clower   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
  
  log10Pe = expression(paste("Expected -log"[10], "(", plain(p), ")" ))
  log10Po = expression(paste("Observed -log"[10], "(", plain(p), ")" ))
  ggplot(df) +
    geom_point(aes(expected, permuted), shape=3, size=1.5, color = "cornflowerblue") +
    geom_abline(intercept=0, slope=1, alpha=0.5, color = "darkgrey") +
    geom_line(aes(expected, cupper), linetype=2, color = "darkgrey") +
    geom_line(aes(expected, clower), linetype=2, color = "darkgrey") +
    #geom_point(aes(expected, observed), size=2) +
    geom_point(data=cor_res, aes(expected, -log10(sort(p)), color = sig), size=2) +
    scale_color_manual(values=c("black","red")) +
    xlab(log10Pe) + ylab(log10Po) + labs(title=title) + 
    theme_classic()
}

round_df <- function(df, digits) {
  nums <- names(which(vapply(df, is.numeric, FUN.VALUE = logical(1))))
  df[,(nums) := round(.SD,digits), .SDcols=nums]
  return(df)
}