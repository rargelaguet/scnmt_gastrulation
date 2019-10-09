library(data.table)
library(purrr)
library(MotifDb)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(furrr)


# this script downloads position weight matrices (PWM) from an online repository
# (MotifDb) then scans the mouse genome, selecting the top hits for each PWM

# in/out
io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$annotation_out <- paste0(io$data_dir, "/features/sjc/motifs/motifdb_short")

# options
opts <- list()
opts$n_loci <- 1e5 # number of loci per annotation (selecting highest scoring)
opts$min_score <- "80%" # string of format "80%", selecting only sequences that are xx% of the maximum possible score for that PWM
opts$adjust_size <- FALSE # FALSE or integer size in bp of the resulting annotations

# load genome and remove shit chromosomes
chrs <- paste0("chr", c(1:19, "X", "Y"))
mm10 <- BSgenome.Mmusculus.UCSC.mm10
mm10@user_seqnames <- setNames(chrs, chrs)
mm10@seqinfo <- mm10@seqinfo[chrs]


# load all PWMs in mouse genome
pwm <- query(MotifDb, "mmusculus")




dir.create(io$annotation_out, recursive=TRUE)

# subset to only include files not already generated (useful if script crashes)
files <- paste0("Mmusculus-", dir(io$annotation_out)) %>%
  sub(".bed", "", .)
select <- names(pwm)[!names(pwm) %in% files]
pwm <- pwm[select]

plan(multiprocess)

future_map2(pwm, names(pwm), ~{
  anno <- sub("Mmusculus-", "", .y)
  file <- paste0(io$annotation_out, "/", anno, ".bed")
  
  dt <- matchPWM(.x, mm10, with.score=TRUE, opts$min_score) %>%
    as.data.frame() %>%
    setDT() %>%
    # select top sites by score
    .[order(-rank(score))] %>%
    .[1:opts$n_loci] %>%
    .[!is.na(start)] %>%
    .[, .(chr = sub("chr", "", seqnames), start, end, strand, id = paste0(anno, "_", .I), anno = anno, score)]

  if (opts$adjust_size!=FALSE) {
    dt[, mid := start+(end-start)/2 %>% round]
    dt[, c("start", "end", "mid") := .(mid-opts$adjust_size/2, 
                                       mid+opts$adjust_size/2,
                                       NULL)]
  }
  # save as bed file
  fwrite(dt, file, sep="\t")
  print(paste0("saving ", file))
})






