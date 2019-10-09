library(data.table)
library(purrr)

sample_metadata <- fread("/Users/ricard/data/gastrulation/sample_metadata.txt") %>%
  .[stage=="E7.5" & lineage10x=="Visceral_endoderm",lineage10x_2:="Visceral_endoderm"] %>%
  .[!is.na(lineage10x_2) ]

colors_lineage10x_2 <- c(
  "Epiblast" = "grey50",
  "Ectoderm" = "steelblue",
  "ExE ectoderm" = "steelblue",
  "Primitive Streak" = "sandybrown",
  "Mesoderm" = "#CD3278",
  "Primitive endoderm" = "#2E8B57",
  "Visceral endoderm" = "#2E8B57",
  "Endoderm" = "#43CD80"
)

colors_lineage10x <- c(
  "Epiblast" = "grey50",
  "Rostral neurectoderm" = "steelblue",
  "Surface ectoderm" = "steelblue",
  "ExE ectoderm" = "steelblue",
  "Caudal epiblast" = "sandybrown",
  "Anterior Primitive Streak" = "sandybrown",
  "Primitive Streak" = "sandybrown",
  "Nascent mesoderm" = "#FF82AB",
  "Mixed mesoderm" = "#CD3278",
  "Mature mesoderm" = "#CD3278",
  "Intermediate mesoderm" = "#CD3278",
  "Pharyngeal mesoderm" = "#CD3278",
  "Paraxial mesoderm" = "#CD3278",
  "Somitic mesoderm" = "#CD3278",
  "Haematoendothelial progenitors" = "#CD3278",
  "Caudal mesoderm" = "#CD3278",
  "ExE mesoderm" = "#CD3278",
  "Mesenchyme" = "#CD3278",
  "Primitive endoderm" = "#2E8B57",
  "Visceral endoderm" = "#2E8B57",
  "Parietal endoderm" = "#2E8B57",
  "Gut" = "#43CD80",
  "Notochord" = "#43CD80",
  "Embryonic endoderm" = "#43CD80"
)


to.plot <- sample_metadata[,.N, by=c("stage","lineage10x_2")] %>% 
  .[, lineage10x_2:=stringr::str_replace_all( lineage10x_2,"_"," ")] %>%
  .[, lineage10x_2:=factor( lineage10x_2,levels=names(colors_lineage10x_2))]
  
p <- ggplot(to.plot, aes(x= lineage10x_2, y=N)) +
  geom_bar(aes(fill= lineage10x_2), stat="identity", color="black") +
  scale_fill_manual(values=colors_lineage10x_2) +
  facet_wrap(~stage, nrow=1) +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.3)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

pdf("/Users/ricard/data/gastrulation/mapping_10x/mapping_stats/mapping_lineage10x_2.pdf", width=11, height=5)
print(p)
dev.off()
