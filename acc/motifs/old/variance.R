library(data.table)
library(purrr)
library(weights)
library(ggplot2)
library(cowplot)


io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/gastrulation_data"
io$anno_dir <- paste0(io$data_dir, "/features/sjc/motifs/motifdb/filt")

# acc <- data.table(id=paste0("id", 1:10000),
#                  rate=round(runif(10000, 0, 100)),
#                  weight=round(runif(10000, 1, 10)),
#                  sample=c("A", "B", "C", "D", "E")
#                  )
# 
# background_anno= data.table(id=paste0("id", 1:10000),
#                             chr=round(runif(10000, 1, 19)),
#                             start=round(runif(10000, 1e4, 1e8)),
#                             anno=c("A", "B", "C", "D", "E")
# ) %>% .[, end := start + 100]

acc <- dir(io$parsed_acc, full=TRUE, pattern=".tsv.gz") %>% 
    map(~fread(.) %>% .[, .(var=wtd.var(rate, weight), type="tfs"), .(id, anno)]
    )

# all loci are 100bp in size so background should be x sites of 100bp

background_anno <- dir(io$anno_dir, full=TRUE, pattern=".bed") %>% 
    map(fread) %>%
    rbindlist() %>% 
    .[, c("min", "max") := .(min(start), max(end)), chr] %>% 
    .[, .(start = runif(1, min, max)), .(id, anno)] %>% 
    .[, end := start+100]

acc_files <- dir(io$raw_acc, full=TRUE, pattern=".tsv.gz")
back_acc <- map(opts$cells, ~{
    file <- files[grep(., files)] %>% paste("zcat", .)
    fread(file) %>% 
        .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>% 
        setkey(chr, start, end) %>% 
        foverlaps(background_anno %>% setkey(chr, start, end), nomatch=0L) %>% 
        .[, .(rate=round(mean(rate*100)), weight=.N, sample=.), .(id, anno)]
}) %>% 
    rbindlist() %>% 
    .[, .(back_var=wtd.var(rate, weight), type="background"), .(id, anno)]


to_plot <-rbind(acc, back_acc)

ggplot(to_plot, aes(anno, var, fill=type)) +
    geom_boxplot(outlier.shape=NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.2))

top_var <- to_plot[, .(var=mean(var)), .(id, anno, type)] %>% 
    dcast(id+anno~type, value.var="var")


