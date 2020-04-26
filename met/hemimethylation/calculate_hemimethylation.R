library(BSgenome.Mmusculus.UCSC.mm10)
library(data.table)
library(purrr)

#####################
## Define settings ##
#####################


################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation"
} else {
  stop("Computer not recognised")
}  
io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$data <- paste0(io$basedir,"/met/cpg_level")
io$features <- paste0(io$basedir,"/features/filt")
io$outdir <- paste0(io$basedir,"/met/results/hemimethylation")

####################
## Define options ##
####################

opts <- list()

# Define genomic contexts (use NULL for no genomic context filtering)
opts$annos <- c(
  "genebody",
  "prom_2000_2000_cgi",
  "prom_2000_2000_noncgi",
  # "prom_2000_2000"
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500",
  # "exons",
  # "introns",
  # "CGI",
  # "LINE",
  "LTR" = "LTR"
)
# opts$annos <- NULL

# Define which cells to use
opts$cells <- fread(io$metadata) %>% 
  .[!is.na(id_met) & pass_metQC==TRUE,id_met]

# Output lists
stats_short <- list()
stats_full <- list()

###############
## Load data ##
###############

# Load genomic contexts metadata
if (!is.null(opts$annos)) {
  anno_dt <- lapply(opts$annos, function(i) fread(sprintf("%s/%s.bed.gz",io$features,i))) %>%  
    rbindlist %>% setnames(c("chr","start","end","strand","id","anno")) %>%
    .[,strand:=NULL] %>%
    setkey(chr,start,end)
}

# Load methylation data and calculate hemimethylation 
for (i in opts$cells) {
  file <- sprintf("%s/%s.tsv.gz",io$data,i)
  if (file.exists(file)) {
    print(sprintf("Loading %s methylation...",i))

    # Load sample methylation data
    data <- fread(file, showProgress = F, select=c(1,2,5), 
                  colClasses=c("chr"="factor", "pos"="integer", "rate"="numeric")) %>%
      .[rate%in%c(0,1)] %>% .[!chr%in%c("MT")] %>% 
      setorder(chr,pos)
    total_number_cpgs <- nrow(data)
    
    print(table(data$chr))
    
    # Get genomic sequence
    seq <- unname( as.character(getSeq(Mmusculus, paste0("chr",data$chr), data$pos-1, data$pos+1)) )
    data[,c("base_up","base","base_down") := list(substr(seq,1,1),substr(seq,2,2),substr(seq,3,3))]
    
    print(table(data$base))
    
    data <- data %>%
      .[,strand:=ifelse(base=="C","+","-")] %>%   # Add strand information
      .[,pos:=ifelse(strand=="-",pos,pos+1)] %>%  # Convert all positions to the positive strand
      .[,N:=.N,by=c("chr","pos")] %>%
      .[N==2,.(rate=mean(rate)),c("chr","pos")]
    
    print(table(data$rate))

    # Compute genome-wide statistics
    stats_short[[i]] <- data.table(
      sample = i, 
      anno = "all", 
      number_cpgs = nrow(data),
      # fraction_cpgs = nrow(data)/total_number_cpgs
      hemimethylation = mean(!data$rate%in%c(0,1))
    )
    
    # Calculate statistics per genomic context
    data[,c("start","end") := list(pos,pos)] %>% setkey(chr,start,end)
    
    if (!is.null(opts$annos[1])) {
      data_k <- foverlaps(data, anno_dt, nomatch=0) %>%
        .[,c("start","end","i.start","i.end"):=NULL] %>%
        .[,sample:=as.factor(i)]
      
      stats_full[[i]] <- data_k %>%
        .[,.(number_cpgs=.N, hemimethylation=mean(!rate%in%c(0,1))), by=c("sample","id","anno")]
      
      stats_short[[i]] <- rbind(stats_short[[i]],
        data_k[,.(number_cpgs=.N, hemimethylation=mean(!rate%in%c(0,1))), by=c("sample","anno")]
      )
    }
  } else {
    print(sprintf("Sample %s not found for methylation",i))
  }
}

fwrite(rbindlist(stats_short), paste0(io$outdir,"/hemimethylation_short.txt.gz"), sep="\t")
fwrite(rbindlist(stats_full), paste0(io$outdir,"/hemimethylation_full.txt.gz"), sep="\t")
  