library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)



io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$fimo <- paste0(io$data_dir, "/features/fimo/H3K27ac_distal_E7_5_Ect_intersect12.bed.fa/fimo.txt")
io$meta <- paste0(io$data_dir, "/sample_metadata_scNMT.txt")
io$pwm <- paste0(io$data_dir, "/features/sjc/PWM/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt")
io$met<- paste0(io$data_dir, "/met/parsed/H3K27ac_distal_E7.5_Ect_intersect12.tsv.gz")
io$cg_density <- paste0(io$data_dir, "/features/cg_density/cpg_density_perfeature.txt")

opts <- list()
opts$stage <- "E7.5"
opts$lineage <- "Ectoderm"
opts$q_cutoff <- 0.05


meta <- fread(io$meta) %>%
  .[stage == "E6.75", stage := "E6.5"] %>%
  .[pass_metQC == TRUE & KO_3b == "not"]

cells <- meta[stage %in% opts$stage & lineage %in% opts$lineage, sample]

pwm <- readLines(io$pwm) %>%
  .[grep("SP1", .)] %>%
  strsplit(" ") %>%
  map_chr(2)

ids_with_motif <- fread(io$fimo, select = c(1:2, 8)) %>%
  setnames(c("motif", "id", "q")) %>%
  .[q < opts$q_cutoff] %>%
  .[motif == pwm, id]
  
met <- fread(paste("zcat", io$met), select = c(1:2, 6)) %>%
  setnames(c("sample", "id", "rate")) %>%
  setkey(id)

met_stage <- merge(met, meta[, .(sample, stage, lineage)], by = "sample") %>%
  .[, .(mean_rate = mean(rate)), .(id, stage, lineage)]

early <- met_stage[stage == "E5.5" & lineage == "EPI" & mean_rate < 40, unique(id)]
early_low <- met_stage[stage == "E5.5" & lineage == "EPI" & mean_rate > 40, unique(id)] 
late <- met_stage[id %in% early_low & stage == "E7.5" & mean_rate < 40, unique(id)]
  

met_split <- met[sample %in% cells] %>%
  .[, motif := "without_motif"] %>%
  .[id %in% ids_with_motif, motif := "with_motif"] %>%
  .[id %in% early, early_late := "early"] %>%
  .[id %in% late, early_late := "late"] %>%
  .[, .(mean_rate = mean(rate)), .(id, motif,early_late)]


ggplot(met_split, aes(motif, mean_rate, fill = motif)) +
  geom_boxplot() +
  facet_wrap(~early_late)

cg_density <- fread(io$cg_density)

toplot <- merge(met_split, cg_density, by = "id")

ggplot(toplot, aes(cpg_density, fill = motif)) +
  geom_density(alpha=0.5)



