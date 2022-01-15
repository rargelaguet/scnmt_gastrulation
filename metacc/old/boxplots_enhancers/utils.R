# Define lineage colors
colors <- c(
  "E4.5 Epiblast"="#C1CDCD",
  "E5.5 Epiblast"="#C1CDCD",
  "E6.5 Epiblast"="#C1CDCD",
  "E7.5 Epiblast"="#C1CDCD",
  "E7.5 Ectoderm"="steelblue",
  "E6.5 Primitive Streak"="sandybrown",
  "E7.5 Primitive Streak"="sandybrown",
  "E7.5 Endoderm"="#43CD80",
  "E7.5 Mesoderm"="#CD3278"
)

theme_pb <- function() {
  theme(
    plot.title = element_text(size=rel(1.2), color="black", hjust=0.5),
    axis.text.x = element_text(size=rel(1.4), color="black", angle=50, vjust=1, hjust=1),
    axis.text.y = element_text(size=rel(1.2), color="black"),
    axis.title.y = element_text(size=rel(1.4), color="black"),
    strip.background = element_rect(fill="#F37A71"),
    strip.text = element_text(size=rel(1.8), color="black"),
    legend.position="top",
    legend.title = element_blank()
  )
}