library(GenomeGraphs) # this package has been depreciated..............

library(data.table)
library(purrr)
library(biomaRt)

# mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

gene.metadata <- fread("/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt")
# outdir <- "/Users/ricard/data/gastrulation_norsync_stuff/metaccrna/dynamics/pluripotency_markers/genomic_location"
outdir <- "/Users/ricard/data/Ecker_2017/mouse/scMET/revision1/running_window/output/selected_hits"

gene_name <- "Satb2"

tmp <- gene.metadata[symbol==gene_name]
chr <- stringr::str_replace(tmp$chr,"chr","")
start <- tmp$start
end <- tmp$end
strand <- tmp$strand

# window <- round( (end-start)/10  )
# window <- 25000
window <- 5000




# Query bioMart data base


# Create Gene object
# gene <- makeGene(id = "ENSMUSG00000029304", type="ensembl_gene_id", biomart = mart)

plusStrand <- makeGeneRegion(chromosome=chr, start=start-window, end=end+window, strand=strand, biomart = mart)
genomeAxis <- makeGenomeAxis(add53 = TRUE)


pdf(sprintf("%s/%s.pdf",outdir,gene_name), useDingbats = F, onefile = F, width=7, height=0.9)
gdPlot(list(plusStrand,genomeAxis))
dev.off()


