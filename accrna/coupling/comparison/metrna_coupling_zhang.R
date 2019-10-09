library(data.table)
library(purrr)
library(SingleCellExperiment)
library(scater)
source("/Users/ricard/NMT-seq_ESC/correlations/utils.R")

#########
## I/O ##
#########

io <- list()

# Zhang2017
io$zhang.rna <- "/Users/ricard/data/Zhang_2017/rna/GSE76505_RNA-FPKM_for_in_vivo_tissues.txt.gz"
io$zhang.sample_metadata <- "/Users/ricard/data/Zhang_2017/sample_metadata.txt"
io$zhang.met <- "/Users/ricard/data/Zhang_2017/methylation/GSE76505_STEM-seq_methylC_mCG/liftover/parsed"
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$annos_dir  <- "/Users/ricard/data/gastrulation/features/filt"

# Output directory
io$outdir <- "/Users/ricard/gastrulation/metrna/coupling/bulk_comparison/out"

#############
## Options ##
#############

opts <- list()

# Define which stage and lineages to look at in Zhang's data
opts$zhang.stage_lineage <- c(
  "E4.0_ICM",
  "E5.5_EPI",#"E5.5_PE", 
  "E6.5_EPI", "E6.5_PS", #"E6.5_VE", 
  "E7.5_Ectoderm", "E7.5_Mesoderm" # "E7.5_Endoderm"
)

# Genomic contexts
opts$annos <- c("genebody", "prom_2000_2000_noncgi", "CGI", "LINE", "LTR")

# Filtering parameters
# opts$min.CpGs <- 3       # Minimum number of CpGs per feature
opts$gene_window <- 1e4  # window length for the overlap between genes and features

# Multiple testing options
opts$threshold_fdr  <- 0.01   # pvalue threshold for significance in FDR

# Correlation type options
opts$method <- "pearson"      # correlation type
opts$weight <- FALSE          # weighted correlation? 

##########################
## Load sample metadata ##
##########################

# Zhang
sample_metadata <- fread(io$zhang.sample_metadata, stringsAsFactors=T) %>% 
  .[id_rna%in%opts$zhang.stage_lineage] %>% 
  .[,stage_lineage:=as.factor(paste(stage,lineage,sep="_"))]

###################
## Load RNA data ##
###################

# Zhang
rna_matrix <- read.table(io$zhang.rna, header=T) %>% tibble::column_to_rownames("Gene") %>%
  .[,colnames(.) %in% sample_metadata$stage_lineage]

###########################
## Load Methylation data ##
###########################

met_dt <- lapply(opts$annos, function(n) fread(sprintf("zcat < %s/%s.tsv.gz",io$zhang.met,n), stringsAsFactors=T)) %>% rbindlist
colnames(met_dt) <- c("sample","id","anno","rate")

######################################
## Load genomic context information ##
######################################

feature_metadata <- lapply(opts$annos, function(anno) 
  fread(sprintf("%s/%s.bed",io$annos_dir,anno), stringsAsFactors=T)[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

# Load gene metadata 
gene_metadata <- fread(io$gene_metadata,stringsAsFactors=T) %>% 
  setnames(c("ens_id","symbol"),c("id","gene")) %>% 
  .[,chr:=as.factor(substr(chr,4,4))]

####################
## Parse the data ##
####################

## RNA ##

# Log transofmration
rna_matrix.log <- log(rna_matrix+1)

# Convert to data.table
rna_dt <- rna_matrix.log %>% t %>% as.data.table(keep.rownames = "sample") %>% 
  melt(id.vars = "sample", value.name = "expr", variable.name = "gene")

## Methylation ##
met_dt <- met_dt[!is.na(rate)]

#############################################################
## Associate the non-genic contexts with overlapping genes ##
#############################################################

feature_metadata_filt <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% met_dt[anno==y,id]] ) %>%
  rbindlist

gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","gene")] %>% 
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)

met_list <- list()
for (ann in unique(met_dt$anno)){
  
  # Subset corresponding anno
  met_tmp <- met_dt[anno == ann, ]
  
  # Non gene-associated feature
  if (all(grepl("ENSMUSG", unique(met_tmp$id)) == FALSE)) {
    
    # Extract coordiantes for methylation sites and for genes
    feature_metadata_tmp <- feature_metadata_filt[anno==ann, c("chr","start","end","id")] %>% 
      .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
    
    # Do the overlap
    ov <- foverlaps(gene_metadata_filt, feature_metadata_tmp, nomatch=0) %>% .[,c("gene", "id")]# %>%
    # # If a feature overlap with more than one gene, concatenate all of them
    # .[,.(gene=paste(gene, collapse="_")),by="id"] %>%
    # .[,c("gene", "id")]
    
    # Merge with methylation data
    met_list[[ann]] <- merge(met_tmp, ov, by="id", allow.cartesian=T) 
  }
  # Gene-associated feature
  else if (all(grepl("ENSMUSG", unique(met_tmp$id)) == TRUE)) {
    met_list[[ann]] <- merge(met_tmp, gene_metadata[,c("id","gene")], by="id")
  }
}
met_dt <- rbindlist(met_list)
rm(met_list,met_tmp,ov,feature_metadata_tmp)

###################################################
## Merge DNA methylation and RNA expression data ##
###################################################

metrna_dt <- merge(
  met_dt, 
  rna_dt[,c("sample","gene","expr")], 
  by=c("sample","gene")
)

################################
## Merge with sample metadata ##
################################

metrna_dt <- metrna_dt %>% merge(sample_metadata[,by=c("stage","stage_lineage")], by="sample")

######################
## Compute coupling ##
######################

if (opts$weight == TRUE){
  if (opts$method != "pearson") { print("Weighted correlation only supported for pearson"); stop() }
  cor <- metrna_dt[, wtd.cor(rate, expr, N)[, c("correlation", "t.value", "p.value")], by = c("sample", "anno")]
  
  # Unweighted correlation
} else{
  cor <- metrna_dt[, .(V1 = unlist(cor.test(rate, expr, alternative = "two.sided", method = opts$method)[c("estimate", "statistic", "p.value")])), by = c("sample", "anno")]
}

# Compute adjusted p-values (both FDR and Bonferroni)
cor <- cor %>% .[,para := c("r", "t", "p")] %>% dcast(sample + anno ~ para, value.var = "V1") %>% 
  .[, c("padj_fdr", "padj_bonf") := list(p.adjust(p, method = "fdr"), p.adjust(p, method = "bonferroni")), by = anno] %>%
  .[, c("log_padj_fdr", "log_padj_bonf") := list(-log10(padj_fdr), -log10(padj_bonf))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% setorder(padj_fdr)

##########
## Plot ##
##########

# Divergent
opts$colors <- c(E4.0="#F8766D", E5.5="#7CAE00", E6.5="#00BFC4", E7.5="#C77CFF")


subset_annos <- c(
  "prom_2000_2000_noncgi"="Promoters",
  "genebody"="Gene\nbodies",
  # "H3K27ac_distal_E7.5_common_500"="Distal\nH3K27ac",
  # "H3K27ac_promoter_E7.5_common_500"="Proximal\nH3K27ac",
  # "H3K27ac_distal_E7.5_union_500"="E7.5 Distal\nH3K27ac",
  # "H3K27ac_promoter_E7.5_union_500"="E7.5 Proximal\nH3K27ac",
  # "H3K4me3_E7.5_common"="H3K4me3",
  # "H3K4me3_E7.5_union"="E7.5\nH3K4me3",
  # "ESC_p300"="p300 ChIP-seq (ESCs)",
  # "ESC_CTCF"="CTCF ChIP-seq (ESCs)",
  "CGI"="CpG\nislands",
  "LINE"="LINE",
  "LTR"="LTR"
)


tmp <- merge(cor, sample_metadata[,c("sample","stage","stage_lineage")], by="sample") %>%
  .[,.(r=mean(r)),by=c("stage","anno")] %>%
  .[anno%in%names(subset_annos)] %>%
  .[,anno:=stringr::str_replace_all(anno, subset_annos)] %>%
  .[,anno:=factor(anno,levels=subset_annos)]


# p <- ggplot(tmp, aes(x=anno, y=r, fill=stage)) +
p <- ggplot(tmp, aes(x=anno, y=r, fill=stage)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(-0.25,0.1)) +
  labs(title="", x="", y="Met/RNA correlation") +
  boxplot_theme() +
  # scale_fill_brewer(palette = "Accent") +
  theme(
    legend.position=c(0.3,1.05),
    axis.text.x = element_text(color="black", size=11)
  )
print(p)

pdf(paste0(io$outdir,"/metrna_coupling_barplots.pdf"), width=9, height=5, useDingbats = F)
print(p)
dev.off()