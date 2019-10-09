

#######################################
## Preprocessing of NMT-seq data ##
#######################################

Steps:
(1) Prepare .cov files for CG methylation and GC accessibility
	prepare_cov.sh
(2) Merge R1 and R2 files:
	merge.R
(3) (Optional) Filter non-CG or non-GC dinucleotides:  this is not required because this step is done by Bismark (report.txt -> .cov). If we find missmatches they are more likely to be true CG or GC dinucleotides that are in disagreement with the reference genome
	filter.R
(4) (Optional) Remove strnad-specific information: convert all CpG/GpC sites to the positive strand
	assign_strand.R
(5) Binarise CpG sites and calculate rates
	binarise.R
(6) Collapse single sites into features using known annotations:
	annotate.R
