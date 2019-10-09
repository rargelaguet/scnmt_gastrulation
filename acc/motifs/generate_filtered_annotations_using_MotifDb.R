library(data.table)
library(purrr)
library(MotifDb)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)

# this script downloads position weight matrices (PWM) from an online repository
# (MotifDb) then scans the mouse genome, selecting the top hits for each PWM and 
# adjusts the start and end positions to a given window size then saves bed files

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$base_dir <- "/bi/group/reik/Stephen/gastrulation"
io$annotation_out <- paste0(io$data_dir, "/features/sjc/motifs/motifdb/filt/motif_only")
io$anno_subset <- paste0(io$base_dir, "/acc/motifs/profiles/out/anno_selection.rds")

opts <- list()
opts$min_score <- "80%"
opts$max_loci <- 100000
opts$extend_window <- FALSE # FALSE or numeric size of window
  

# load genome and remove shit chromosomes
chrs <- paste0("chr", c(1:19, "X", "Y"))
mm10 <- BSgenome.Mmusculus.UCSC.mm10
mm10@user_seqnames <- setNames(chrs, chrs)
mm10@seqinfo <- mm10@seqinfo[chrs]


# load all PWMs in mouse genome and filter for our subset
subset <- readRDS(io$anno_subset)
pwm <- query(MotifDb, "mmusculus") %>%
  .[subset]

# register parallel
cores <- detectCores()
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
pckgs <- c("data.table", "purrr", "MotifDb", "Biostrings", "BSgenome")

dir.create(io$annotation_out, recursive=TRUE)


foreach(i=seq_along(pwm), .packages=pckgs) %dopar% {
    anno <- names(pwm[i]) %>%
        sub("Mmusculus-", "", .)
    file <- paste0(io$annotation_out, "/", anno, ".bed")
    
    dt <- matchPWM(pwm[[i]], mm10, with.score=TRUE, min.score=opts$min_score) %>%
        as.data.frame() %>%
        setDT() %>%
        setnames("seqnames", "chr") %>%
        .[, c("string", "width") := NULL] %>%
        # select top sites by score
        .[order(-rank(score))] %>%
        .[1:min(.N, opts$max_loci)] %>%
        .[, chr := sub("chr", "", chr)] %>%
        .[, anno := anno] %>%
        .[, id := paste0(anno, "_", .I)] %>%
        setcolorder(c("chr", "start", "end", "strand", "id", "anno", "score"))
    # extend window size
    if (opts$extend_window!=FALSE){
      dt[, mid := start+(end-start/2)] %>%
        .[, c("start", "end", "mid") := .(mid-opts$extend_window/2, mid+opts$extend_window/2, NULL)]
    }
    
    # save as bed file
    fwrite(dt, file, sep="\t")
    anno
}

stopCluster(cl)

