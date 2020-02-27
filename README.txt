#########################
## General description ##
#########################

This file contains the in vivo gastrulation data set (E4.5, E5.5, E6.5, E7.5) sequenced with scNMT-seq from Argelaguet et al. Nature. 2019.

#################################
## Sample metadata information ##
#################################

Cells are assigned different IDs for RNA expression (id_rna) DNA methylation (id_met) and chromatin accessibility (id_acc). Cells are identified in the different data modalities by the corresponding ID. The column "sample" contains an ID that can be used to match cells across data modalities.

Columns:
- sample: cell id used to anchor the three modalities
- id_rna: cell id for RNA data 
- id_met: cell id for DNA methylation data
- id_acc: cell id for chromatin accessibility data
- embryo: embryo ID
- plate: plate ID
- pass_rnaQC: does the cell pass QC for RNA expression? If NA, the RNA expression was not sequenced for that cell
- pass_metQC: does the cell pass QC for DNA methylation? If NA, the DNA methylation was not sequenced for that cell
- pass_accQC: does the cell pass QC for chromatin accessibility? If NA, the chromatin accessibility was not sequenced for that cell
- stage: embryonic stage
- lineage10x: level 1 (maximally resolved) lineage assignment by mapping the RNA expression profiles to the embryo atlas
- lineage10x_2: level 2 (aggregated) lineage assignment by mapping the RNA expression profiles to the embryo atlas


###################
## Files/Folders ##
###################

- features: gene and genomic annotations (descriptions below)
	genomic_contexts: genomic annotations
	genes: gene metadata
	motifs: motif information

- acc: data/results involvign only chromatin accessibility
    gpc_level: binary measurements at the GpC level (0=inaccessible, 1=accessible)
    feature_level: aggregated measurements over genomic contexts

- met: data/results involvign only DNA methylation
    cpg_level: binary measurements at the CpG level (0=unmethylated, 1=methylated)
    feature_level: aggregated measurements over genomic contexts

- rna: data/results involvign only RNA expression
	counts.txt.gz: count matrix
	SingleCellExperiment.rds: normalised SingleCellExperiment object

- accrna: data/results involving chromatin accessibility and RNA expression
	cor: correlation analysis of promoter accessibility and RNA expression

- metrna: data/results involving DNA mehylation and RNA expression
	cor: correlation analysis of promoter methylation and RNA expression

- metaccrna: data/results involving all three omics
	cor: scatterplots of Met/RNA vs Acc/RNA correlation
	differential: differential DNA methylation and chromatin accessibiliy analysis between germ layres
	mofa: MOFA analysis
	plot_dynamics_individual_examples: boxplots of the three omics per gene
	qc_stats: quality control and general statistics
	mesendoderm_commitment: study of the continuous epigenetic dynamics during mesoderm and endoderm commitment

- metacc: data/results involving DNA mehylation and chromatin accessibility
	boxplots:
	pseudobulked_profiles:

- H3K27ac: data/results involving only H3K27ac ChIP-seq data

- sample_metadata.txt: cell metadata (descriptions of the columns below)

########################################
## Description of genomic annotations ##
########################################

Genome assembly: mm10

- CGI: CpG islands (downloaded from UCSC)
- ESC_CTCF (CTCF binding sites in ESCs, downloaded from ENCODE)
- ESC_DHS (DNA-ase hypersensitive sites in ESCs, downloaded from ENCODE)
H3K27ac_distal_E7.5_Ect_intersect12: Ectoderm-specific distal H3K27ac ChIP-seq peaks (putative enhancers). Only peaks that intersect - between two replicates.
- H3K27ac_distal_E7.5_Ect_intersect12_500: same as above, but the peaks are extended by +-500bp to get better coverage
H3K27ac_distal_E7.5_End_intersect12: Endoderm-specific distal H3K27ac ChIP-seq peaks (putative enhancers). Only peaks that intersect - between two replicates
- H3K27ac_distal_E7.5_End_intersect12_500: same as above, but the peaks are extended by +-500bp to get better coverage
H3K27ac_distal_E7.5_Mes_intersect12: Mesoderm-specific distal H3K27ac ChIP-seq peaks (putative enhancers). Only peaks that intersect - between two replicates
- H3K27ac_distal_E7.5_Mes_intersect12_500: same as above, but the peaks are extended by +-500bp to get better coverage
- H3K4me3_E7.5_Ect: Ectoderm-specific H3K4me3 peaks (active transcription start sites)
- H3K4me3_E7.5_End: Endoderm-specific H3K4me3 peaks (active transcription start sites)
- H3K4me3_E7.5_Mes: Endoderm-specific H3K4me3 peaks (active transcription start sites)
- H3K4me3_E7.5_common: H3K4me3 peaks that are common between all three germ layers.
- H3K4me3_E7.5_union: all H3K4me3 peaks irrespective of whether they are unique or shared betweeen germ layers.
- LINE: (...) (repetitive elements)
- LTR: long terminal repeats (repetitive elements)
- exons: gene exons (only protein-coding genes, extracted from Ensembl v87)
- genebody: gene bodies (only protein-coding genes, extracted from Ensembl v87)
- introns: gene introns (only protein-coding genes, extracted from Ensembl v87)
- prom_2000_2000: gene promoters deifned as +-2kb around TSS (only protein-coding genes, extracted from Ensembl v87)
- prom_2000_2000_cgi: same as above, but promoters that overlap with CpG islands
- prom_2000_2000_noncgi: same as above, but promoters that do not overlap with CpG islands
- prom_200_200, prom_200_200_cgi, prom_200_200_noncgi: same as above, but promoters defined as +-200bp around TSS.

Note: in general we find that chromatin accessibility is very pronounced very close to the TSS, around +-200bp, and declines rapidly away from the TSS. In contrast, DNA methylation has broader signal up to +-2kb around the TSS (see Clark2019, Nature communications, Supplementary Figure 22)


##################################
## Data formats: RNA expression ##
##################################

RNA readouts are obtained following the standard Smart-seq2 protocol without UMIs.

We provide two files:
- counts.txt.gz: count matrix
- SingleCellExperiment.rds: Bioconductor SingleCellExperiment object


#################################################
## Data formats: methylation and accessibility ##
#################################################

DNA methylation and chromatin accessibility measurements are obtained using NOME-seq. 

NOME-seq uses a GpC methyltransferase (M.CviPI) to label open chromatin by methylating GpC sites. This enables the simultaneous profiling of endogenous DNA methylation (CpG context) and chromatin accessibility (GpC context) from the same experiment by a single round of bisulfite sequencing.

For every CpG site a DNA methylation rate is calculated as the fraction of methylated reads divided by all reads covering that position.
Similarly, Similarly, for every Gpc site an accessibility rate rate is calculated as the fraction of methylated reads divided by all reads covering that position (see folders met/cpg_level and acc/gpc_level). A similar procedure is used to calculate DNA methylation and accessibility rates at genomic features (see folders met/feature_level and acc/feature_level)

Importantly, this protocol has the drawback that one needs to filter out two trinucleotide contexts: 
- GCGs due to ambiguity between endogenous and enzymatic methylation (21% of all CpGs)
- CCG due to some off-target effects of the enzyme (27% of all CpGs). Note that we did not make a thorough exploration of how strong the off-target effect is. From Kelly2012 it seems rather weak but we have been conservative and removed all CCG sites.

You can read more about NOME-seq and NMT-seq in the following articles:
(bulk) https://genome.cshlp.org/content/22/12/2497.full.html
(single-cell) https://elifesciences.org/articles/23203
(single-cell) https://www.nature.com/articles/s41467-018-03149-4


#######################
## Code availability ##
#######################

https://github.com/rargelaguet/scnmt_gastrulation

#############
## Contact ##
#############

Ricard Argelaguet: ricard@ebi.ac.uk
Stephen Clark: stephen.clark@babraham.ac.uk 
