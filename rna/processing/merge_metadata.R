library(data.table)
library(purrr)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

a <- fread("zcat < /Users/ricard/data/gastrulation/rna/counts/old/counts.txt.gz") %>% matrix.please
b <- fread("zcat < /Users/ricard/data/gastrulation/rna/counts/old/raw_counts_endoderm_p11_to_p16.txt.gz") %>% matrix.please

c <- intersect(rownames(a),rownames(b))

ab <- cbind(a[c,], b[c,])

write.table(ab, "/Users/ricard/data/gastrulation/rna/counts/counts.txt", col.names=T, row.names=T, sep="\t", na="NA", quote=F)

sample_metadata <- fread("/Users/ricard/data/gastrulation/sample_metadata_new_endoderm_dec2018.txt")
# sample_metadata[lineage=="",lineage:=NA]
sample_metadata[!is.na(pass_metQC),id_met:=id_rna]
sample_metadata[!is.na(pass_accQC),id_acc:=id_rna]
fwrite(sample_metadata, "/Users/ricard/data/gastrulation/sample_metadata_new_endoderm_dec2018.txt", col.names=T, row.names=F, sep="\t", na="NA", quote=F)


sample_metadata_good <- fread("/Users/ricard/data/gastrulation/sample_metadata_scNMT.txt")
colnames(sample_metadata_good)

foo <- plyr::rbind.fill(sample_metadata_good, sample_metadata)
fwrite(foo, "/Users/ricard/data/gastrulation/sample_metadata_scNMT2.txt", col.names=T, row.names=F, sep="\t", na="NA", quote=F)


lol <- fread("/Users/ricard/data/gastrulation/old_metadata/sample_metadata_scMT.txt")
lol[,.N,by=c("stage","lineage")]
