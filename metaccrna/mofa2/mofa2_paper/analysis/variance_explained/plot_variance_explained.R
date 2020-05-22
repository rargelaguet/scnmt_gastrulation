io$pdfdir2 <- paste0(io$pdfdir,"/variance_explained")

for (i in factors(model)) {
# for (i in head(factors(model),n=1)) {
  p <- plot_variance_explained(model, x="group", y="view", factor=i, legend = T)
  
  pdf(sprintf("%s/%s_legend.pdf",io$pdfdir2,i), width=4, height=3)
  print(p)
  dev.off()
}
