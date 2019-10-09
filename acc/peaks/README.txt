pseudobulk.R
  - generate pseudobulk accessibility (or methylation) data from a specific
    stage and lineage

peak_caller.R 
  - simple script to quantify accessibility rates at overlapping 
    windows in pseudobulk data then pick the top x% to use as an annotation

fast_annotate.R
  - compute met/acc rates at specified genomic annotations for specified raw 
    met/acc files

unbiased_win_annotation.R
  - generate running window annotation and quantitate met/acc rates for
    specified raw met/acc files. processing is done in chunks due to the large
    size of the files otherwise.

unbiased_diff_acc.R
  - t-tests for differential acc (or met) at unbiased windows. Input is 
    quantified met/acc generated using unbiased_win_annotation.R. Uses 
    permutations to mitigate multiple-testing problem - as such can take a long 
    time to run.
    
