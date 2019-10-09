library(data.table)
library(purrr)
library(SingleCellExperiment)
library(scater)
source("/Users/ricard/NMT-seq_ESC/correlations/utils.R")

#########
## I/O ##
#########

io <- list()

# wu2017
io$wu.rna <- "/Users/ricard/data/gastrulation/public_data/Wu_2016/rna_counts.txt"
io$wu.sample_metadata <- "/Users/ricard/data/gastrulation/public_data/Wu_2016/Wu_2016_ATAC_meta_data.txt"
io$wu.acc <- "/Users/ricard/data/gastrulation/public_data/Wu_2016/Wu_2016_ATAC_parsed.tsv.gz"
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$annos_dir  <- "/Users/ricard/data/gastrulation/features/filt"

# Output directory
io$outdir <- "/Users/ricard/gastrulation/accrna/coupling/bulk_comparison/out"

#############
## Options ##
#############

opts <- list()

# Define which stage and lineages to look at in wu's data
# opts$wu.stage_lineage <- c(
#   "E4.0_ICM",
#   "E5.5_EPI",#"E5.5_PE", 
#   "E6.5_EPI", "E6.5_PS", #"E6.5_VE", 
#   "E7.5_Ectoderm", "E7.5_Mesoderm" # "E7.5_Endoderm"
# )
# early_2-cell
# 2-cell
# 8-cell
# ICM
# 8-cell
# mESC

# Genomic contexts
opts$annos <- c("genebody", "prom_2000_2000_noncgi", "CGI", "LINE", "LTR")

# Filtering parameters
opts$gene_window <- 1e4  # window length for the overlap between genes and features

# Multiple testing options
opts$threshold_fdr  <- 0.01   # pvalue threshold for significance in FDR

# Correlation type options
opts$method <- "pearson"      # correlation type
opts$weight <- FALSE          # weighted correlation? 

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$wu.sample_metadata)# %>% 
  # .[stage%in%opts$wu.stage_lineage] %>% 
  # .[,stage_lineage:=as.factor(paste(stage,lineage,sep="_"))]

###################
## Load RNA data ##
###################

rna_matrix <- read.table(io$wu.rna, header=T) %>% tibble::column_to_rownames("ens_id")
# %>% .[,colnames(.) %in% sample_metadata$stage_lineage]
earlycell2 <- rowMeans(rna_matrix[,c("early_2.cell_rep1_1_val_1","early_2.cell_rep2_1_val_1")]) %>% as.matrix()
cell2 <- rowMeans(rna_matrix[,c("X2.cell_rep1_1_val_1","X2.cell_rep1_1_val_1")]) %>% as.matrix()
cell8 <- rowMeans(rna_matrix[,c("X8.cell_rep1_1_val_1","X8.cell_rep1_1_val_1")]) %>% as.matrix()
icm <- rowMeans(rna_matrix[,c("ICM_rep1_1_val_1","ICM_rep2_1_val_1")]) %>% as.matrix()
rna_matrix <- cbind(earlycell2,cell2,cell8,icm)
colnames(rna_matrix) <- c("early2cell","2cell","8cell","icm")

###########################
## Load Accessibility data ##
###########################

acc_dt <- fread(sprintf("zcat < %s",io$wu.acc), select = c(1,2,3,6))

# IGNORED: c("2-cell_rep2_1_val_1","8-cell_rep2_2__1_val_1","ICM_rep2_1_trimmed") AND THE TRIMMED ONES
samples <- c(
  "early_2-cell_rep2_1_val_1"="early2cell",
  "2-cell_rep1_1_val_1"="2cell",
  "8-cell_rep1_2__1_val_1"="8cell",
  "ICM_rep1_1_val_1"="icm"
)
acc_dt <- acc_dt[sample%in%names(samples)] %>%
  .[,sample:=stringr::str_replace_all(sample, samples)]


######################################
## Load genomic context information ##
######################################

feature_metadata <- lapply(opts$annos, function(anno) 
  fread(sprintf("%s/%s.bed",io$annos_dir,anno))[,c(1,2,3,4,5,6)]) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

# Load gene metadata 
gene_metadata <- fread(io$gene_metadata) %>% 
  setnames(c("ens_id","symbol"),c("id","gene")) %>% 
  .[,chr:=substr(chr,4,4)]

####################
## Parse the data ##
####################

## RNA ##

# Log transofmration
rna_matrix.log <- log(rna_matrix+1)

# Convert to data.table
rna_dt <- rna_matrix.log %>% t %>% as.data.table(keep.rownames = "sample") %>% 
  melt(id.vars = "sample", value.name = "expr", variable.name = "gene")

## Accessibility ##

#############################################################
## Associate the non-genic contexts with overlapping genes ##
#############################################################

feature_metadata_filt <- feature_metadata %>% split(.$anno) %>% 
  map2(.,names(.), function(x,y) x[id %in% acc_dt[anno==y,id]] ) %>%
  rbindlist

gene_metadata_filt <- gene_metadata %>% .[,c("chr","start","end","id")] %>% 
  setnames("id","gene") %>%
  .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)

acc_list <- list()
for (ann in unique(acc_dt$anno)){
  
  # Subset corresponding anno
  acc_tmp <- acc_dt[anno == ann, ]
  
  # Non gene-associated feature
  if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == FALSE)) {
    
    # Extract coordiantes for Accessibility sites and for genes
    feature_metadata_tmp <- feature_metadata_filt[anno==ann, c("chr","start","end","id")] %>% 
      .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>% setkey(chr,start,end)
    
    # Do the overlap
    ov <- foverlaps(gene_metadata_filt, feature_metadata_tmp, nomatch=0) %>% .[,c("gene", "id")]# %>%
    # # If a feature overlap with more than one gene, concatenate all of them
    # .[,.(gene=paste(gene, collapse="_")),by="id"] %>%
    # .[,c("gene", "id")]
    
    # Merge with Accessibility data
    acc_list[[ann]] <- merge(acc_tmp, ov, by="id", allow.cartesian=T) %>%
      setcolorderc(c("sample","id","anno","fpkm"))
  }
  # Gene-associated feature
  else if (all(grepl("ENSMUSG", unique(acc_tmp$id)) == TRUE)) {
    acc_list[[ann]] <- acc_tmp[,gene:=id] %>%
      setcolorderc(c("sample","id","anno","fpkm"))
  }
}
acc_dt <- rbindlist(acc_list)
rm(acc_list,acc_tmp,ov,feature_metadata_tmp)

###################################################
## Merge Accessibility and RNA expression data ##
###################################################

accrna_dt <- merge(
  acc_dt, 
  rna_dt[,c("sample","gene","expr")], 
  by=c("sample","gene")
)

################################
## Merge with sample metadata ##
################################

# accrna_dt <- accrna_dt %>% merge(sample_metadata[,c("sample","stage")], by="sample")
accrna_dt[,stage:=sample]

######################
## Compute coupling ##
######################

cor <- accrna_dt[, .(V1 = unlist(cor.test(fpkm, expr, alternative = "two.sided", method = opts$method)[c("estimate", "statistic", "p.value")])), by = c("sample", "anno")]

# Compute adjusted p-values (both FDR and Bonferroni)
cor <- cor %>% .[,para := c("r", "t", "p")] %>% dcast(sample + anno ~ para, value.var = "V1") %>% 
  .[, c("padj_fdr", "padj_bonf") := list(p.adjust(p, method = "fdr"), p.adjust(p, method = "bonferroni")), by = anno] %>%
  .[, c("log_padj_fdr", "log_padj_bonf") := list(-log10(padj_fdr), -log10(padj_bonf))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% setorder(padj_fdr)

##########
## Plot ##
##########

# Divergent
# opts$colors <- c(E4.0="#F8766D", E5.5="#7CAE00", E6.5="#00BFC4", E7.5="#C77CFF")




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


# tmp <- merge(cor, sample_metadata[,c("sample","stage")], by="sample") %>%
tmp <- cor %>% copy %>% .[,stage:=sample] %>%
  .[,.(r=mean(r)),by=c("stage","anno")] %>%
  .[anno%in%names(subset_annos)] %>%
  .[,anno:=stringr::str_replace_all(anno, subset_annos)] %>%
  .[,anno:=factor(anno,levels=subset_annos)] %>%
  .[,stage:=factor(stage,levels=c("early2cell","2cell","8cell","icm"))]


# p <- ggplot(tmp, aes(x=anno, y=r, fill=stage)) +
p <- ggplot(tmp, aes(x=anno, y=r, fill=stage)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  # scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(-0.25,0.1)) +
  labs(title="", x="", y="Acc/RNA correlation") +
  boxplot_theme() +
  # scale_fill_brewer(palette = "Accent") +
  theme(
    legend.position=c(0.3,1.05),
    axis.text.x = element_text(color="black", size=11)
  )
print(p)

pdf(paste0(io$outdir,"/accrna_coupling_barplots.pdf"), width=9, height=5, useDingbats = F)
print(p)
dev.off()