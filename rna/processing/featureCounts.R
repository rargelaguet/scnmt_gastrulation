##################################################
## Preprocess expression data: feature counting ##
##################################################

# This script is used to create the gene count matrix from the raw BAM files
# It performs step 1 of the following pipeline: 
# (1) Generate raw read count expression matrix
# (2) Normalise read count expression matrix
# (3) Filter genes and samples

# Input:
# - directory with bam files
# - metadata/annotation file in GTF or BED format

# Output
# - matrix of raw counts with dim (ngenes,nsamples)

suppressMessages(library(Rsubread))

####################
## Define options ##
####################

## I/O ##
io <- list()

# io$anno_infile <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/GTF/mRNA/Mus_musculus.GRCm38.87.gtf"
# io$anno_infile <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/GTF/mRNA_ERCC/Mus_musculus_ERCC.GRCm38.87.gtf"
# io$samples_indir <- "/hps/nobackup/stegle/users/ricard/gastrulation/rna/original/E7.5_scNMT_Plate1"
# io$outfile <- "/hps/nobackup/stegle/users/ricard/gastrulation/rna/original/test/raw_counts.txt"

io$anno_infile <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/GTF/mRNA/Mus_musculus.GRCm38.87.gtf"
io$samples_indir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_raw_data/rna/scNMT_symbolic"
io$outfile <- "/homes/ricard/counts_v2.txt"


# io$anno_infile <- "/bi/group/reik/Stephen/gastrulation/rna_preproc/Mus_musculus.GRCm38.87.gtf"
# io$samples_indir <- "/bi/sequencing/Sample_4363_E4.5_E5.5_RNA/Lane_5998_E4.5_E5.5_RNA/Aligned"
# io$outfile <- "/bi/group/reik/Stephen/gastrulation/rna_preproc/raw_counts_2018-02-20.txt"

## Options ##
opts <- list()

# Number of cores to parallelise work
opts$nthreads <- 10

# Samples to use
opts$samples <- "all"
if (opts$samples == "all") {
  opts$samples <- sub(".bam","",list.files(io$samples_indir, pattern=".bam"))
}

# featureCounts options
#opts$anno_format <- "BED6"
opts$anno_format <- "GTF"
opts$allowMultiOverlap = FALSE        # should a read be allowed to be assigned to more than one feature?
opts$isPairedEnd = TRUE               # paired end data? If TRUE, fragments will be counted instead of individual reads
opts$requireBothEndsMapped = FALSE    # should both ends from the same fragment be required to be aligned?
opts$strandSpecific = 0               # 0=unstranded, 1=stranded, 2=reversely stranded
opts$countChimericFragments = FALSE   # should we allow a fragment to have its two reads mapped to different chromosomes?
opts$countMultiMappingReads <- TRUE  # logical indicating if multi-mapping reads/fragments should be counted
if (opts$anno_format == "GTF") {
  opts$isGTFAnnotationFile = TRUE
} else {
  opts$isGTFAnnotationFile = FALSE
}


#####################
## Load annotation ##
#####################

# featureCounts require as input either GTF or SAF format: 
# SAF: five required columns in any order: GeneID, Chr, Start, End, Strand. 

if (opts$isGTFAnnotationFile == FALSE) {
  anno <- read.table(io$anno_infile, header=F,sep="\t",na.strings=c("","NA"))
  colnames(anno) <- c("chr","start","end","strand","ens_id")
  if (opts$anno_format == "BED6") {
    anno_saf <- anno[,c("ensembl_gene_id","chr","start","end","strand")]
    colnames(anno_saf) <- c("GeneID","Chr","Start","End","Strand")
    io$anno_infile <- anno_saf
  }
}


####################
## Filter samples ##
####################

# if(opts$filt_samples ==T) {
#   samples_df <- subset(samples_df,keep=="yes")
#   filenames = sub("^([^.]*).*", "\\1", basename(opts$in_files))
#   opts$in_files <- opts$in_files[filenames %in% samples_df$bam_file]
#   filenames_filt = sub("^([^.]*).*", "\\1", basename(opts$in_files))
#   sample_names = samples_df[match(filenames_filt,samples_df$bam_file),"sample_name"]
# }


#######################
## Run featureCounts ##
#######################


fc <- featureCounts(files=sprintf("%s/%s.bam",io$samples_indir,opts$samples), annot.inbuilt=NULL, annot.ext=io$anno_infile, 
                    isGTFAnnotationFile=opts$isGTFAnnotationFile, useMetaFeatures=TRUE, allowMultiOverlap=opts$allowMultiOverlap, 
                    isPairedEnd=opts$isPairedEnd, requireBothEndsMapped=opts$requireBothEndsMapped, 
                    nthreads=opts$nthreads, strandSpecific=opts$strandSpecific, 
                    countChimericFragments=opts$countChimericFragments)
counts <- as.data.frame(fc$counts)
colnames(counts) <- opts$samples
counts <- data.frame(ens_id=rownames(counts), counts[,gtools::mixedsort(colnames(counts))])


####################
## Analyse output ##
####################

# fc$counts: data matrix containing read counts for each feature or meta-feature for each library
#            dim (nfeatures,nsamples)
# View(fc$counts)

# fc$annotation: data frame with six columns including GeneID, Chr, Start, End and Length. 
#                When read summarization was performed at meta-feature level, each row in the data frame is 
#                a meta-feature and columns in the data frame give the annotation information for the features
#                For each meta-feature, the Length column gives the total length of genomic regions covered 
#                by features included in that meta-feature
#                dim (nfeatures,6)
# View(fc$annotation)

# targets: a character vector giving sample information.
# fc$targets

# stat: dataframe with stats: (eg. ambiguity, multi-mapping, secondary alignment, mapping quality,
#       fragment length, chimera, read duplicate, non-junction and so on),
# fc$stat

##################
## Save results ##
##################

write.table(counts,io$outfile,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

cells = colnames(counts)[2:385]
find_name = function(x){
    splt = strsplit(x, split="_")
    paste(unlist(splt)[4:5], collapse="_")
}
meta = data.frame(id_rna = cells, 
                  sample = sapply(cells, find_name))

write.table(meta,paste0(dirname(io$outfile), "/meta.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
