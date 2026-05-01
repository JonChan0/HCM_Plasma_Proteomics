The overall script to iterate over all the config files is finemap_iterator.sh.

This calls the finemap.smk which performs a set of operations to perform FINEMAP (i.e Bayesian fine-mapping) as well as calling finemap_analysis.R and ldstore2_prep.R scripts as helper scripts.

V2G acts downstream of the FINEMAP, taking the top-ranked causal configuration of SNP(s) at each genomic region and passing it to Open Target Genetics API to run its L2G scoring pipeline to ID the likely linked gene.
