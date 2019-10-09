library(data.table)
library(purrr)
library(dplyr)

out_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data/rna/raw_bam_original"

dir.create(out_dir, recursive = TRUE)

file_names <- fread("/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata_scNMT.txt") %>%
  .[pass_rnaQC==TRUE & KO_3b == "not", .(id_rna, new_name = paste(stage, lineage, .I, sep="_"))]

dirs <- c("/bi/sequencing/Sample_4165_NMT_embryo_Transcriptome_Pool_1//Lane_5753_NMT_embryo_Transcriptome_Pool_1/Aligned/",
          "/bi/sequencing/Sample_4166_NMT_embryo_Transcriptome_Pool_2//Lane_5754_NMT_embryo_Transcriptome_Pool_2/Aligned/",
          "/bi/sequencing/Sample_4167_NMT_embryo_Transcriptome_Pool_3//Lane_5763_NMT_embryo_Transcriptome_Pool_3/Aligned/",
          "/bi/sequencing/Sample_4363_E4.5_E5.5_RNA/Lane_5998_E4.5_E5.5_RNA//Aligned"
          )

all_files <- map(dirs, dir, pattern = ".bam", full = TRUE) %>% unlist()

old_files <- map_chr(file_names[, id_rna], ~{
  all_files[grep(., all_files)]
})

#new_files <- paste0(out_dir, "/", file_names[, new_name], ".bam")
new_files <- paste0(out_dir, "/", basename(all_files))

file.copy(old_files, new_files)
