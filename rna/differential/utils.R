doDiffExpr <- function(sce, sample_metadata) {
  
  # Filter genes with high dropout rate in both groups
  opts$min.cdr <- 0.25
  sce.1 <- sce[,sce$group=="1"]
  cdr.1 <- rowMeans(exprs(sce.1)>0) > opts$min.cdr
  sce.0 <- sce[,sce$group=="0"]
  cdr.0 <- rowMeans(exprs(sce.0)>0) > opts$min.cdr
  sce <- sce[cdr.0 | cdr.1,]
  
  # Filter genes based on biological overdispersion
  # trend = scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
  # decomp = scran::decomposeVar(sce, fit = trend)
  # decomp = decomp[decomp$mean > 1,]
  # decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  # hvg <- rownames(decomp)[decomp$p.value < 0.10]
  # sce <- sce[hvg,]
  
  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce, type="edgeR")
  
  # Define design matrix (with intercept)
  cdr <- colMeans(exprs(sce)>0)
  design <- model.matrix(~cdr+sce$group)
  
  # Estimate dispersions
  sce_edger  <- estimateDisp(sce_edger,design)
  
  # Fit GLM
  fit <- glmQLFit(sce_edger,design)
  
  # Likelihood ratio test
  lrt <- glmQLFTest(fit)
  
  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
    setnames(c("id","logFC","logCPM","LR","p.value","padj_fdr")) %>%
    .[,c("logCPM","LR"):=NULL] %>%
    .[, c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
    .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>opts$min.logFC)]
}

MinMeanSEMMax <- function(x) {
  v <- c(mean(x) - sd(x)/sqrt(length(x)), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), mean(x) + sd(x)/sqrt(length(x)))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


gg_volcano_plot <- function(tmp, top_genes=10, xlim=NULL, ylim=NULL) {
  negative_hits <- tmp[sig==TRUE & logFC<0,id]
  positive_hits <- tmp[sig==TRUE & logFC>0,id]
  all <- nrow(tmp)
  
  if (is.null(xlim))
    xlim <- max(abs(tmp$logFC), na.rm=T)
  if (is.null(ylim))
    ylim <- max(-log10(tmp$p.value), na.rm=T)
  
  p <- ggplot(tmp, aes(x=logFC, y=-log10(p.value))) +
    labs(title="", x="Log Fold Change", y=expression(paste("-log"[10],"(p.value)"))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange") +
    ggrastr::geom_point_rast(aes(color=sig), size=1) +
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-xlim-2,xlim+2)) +
    scale_y_continuous(limits=c(0,ylim+1)) +
    annotate("text", x=0, y=ylim+1, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-8, y=ylim+1, size=7, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=8, y=ylim+1, size=7, label=sprintf("%d (+)",length(positive_hits))) +
    ggrepel::geom_text_repel(data=head(tmp[sig==T],n=top_genes), aes(x=logFC, y=-log10(p.value), label=symbol), size=5) +
    theme_bw() +
    theme(
      plot.title=element_text(size=28, face='bold', margin=margin(0,0,10,0), hjust=0.5),
      axis.text=element_text(size=rel(1.75), color='black'),
      axis.title=element_text(size=rel(1.95), color='black'),
      axis.title.y = element_text(margin=margin(0,10,0,0)),
      axis.title.x = element_text(margin=margin(10,0,0,0)),
      legend.position="none",
      # panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      # panel.background = element_blank()
    )
  return(p)
}

# scatter_theme <- function(){
#   p <- theme(
#       plot.title=element_text(size=28, face='bold', margin=margin(0,0,10,0), hjust=0.5),
#       axis.text=element_text(size=rel(1.75), color='black'),
#       axis.title=element_text(size=rel(1.95), color='black'),
#       axis.title.y = element_text(margin=margin(0,10,0,0)),
#       axis.title.x = element_text(margin=margin(10,0,0,0)),
#       legend.position="none",
#       panel.border=element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank()
#     )
# }

scatter_theme <- function() {
  p <- theme(
    axis.title.y = element_text(colour="black", size=rel(1.0)),
    axis.title.x = element_text(colour="black", size=rel(1.0)),
    axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.line = element_line(colour="black", size=rel(0.9)),
    axis.ticks = element_line(colour="black", size=rel(1.0)),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position="right",
    legend.text=element_text(size=15),
    legend.key = element_blank(),
    legend.title=element_text(size=17),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
}

boxplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=rel(1.2), hjust=0.5),
    axis.title.y = element_text(colour="black", size=rel(1.1), vjust=1.5),
    axis.title.x = element_text(colour="black", size=rel(1.1), vjust=1.5),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=40, hjust=1),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
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



matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

sample_colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}