These scripts are used to generate annotations to use for quantifying 
accessibility data. The idea is to use position weight matrices (PWMs) for as 
many transcription factors as we can find, scan the genome for instances of 
these PWMs then plot accessibility at these sites to see which PWMs have peaks
of accessibility (or other interesting patterns) which may indicate functional
relevance. 

1. generate_annotations_using_MotifDb.R
2. generate_profiles.R
3. annotation_curation.Rmd
4. generate_filtered_annotations_using_MotifDb.R
5. annotate_acc_data.R

These annotated data can then be used in downstream analysis such as:

6. acc_dimReduction.Rmd
7. local_cor.Rmd 
- this script plots the correlation coefficient between accessibility and 
expression of the nearest gene across loci in windows within the annotation in 
order to visuaslise where a functional region might be (and if this varies 
between stages or lineages).