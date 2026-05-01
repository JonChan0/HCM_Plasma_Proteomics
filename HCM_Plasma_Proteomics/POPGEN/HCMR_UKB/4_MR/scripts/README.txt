There are 2x sets of scripts

1) map_snptest_rsid.pl and rsid_mapper_pvalfilter which work on the meta-analysis results from the GWAS of HCM as per Tadros et al, 2023 
to map on rsIDs

2) Mendelian Randomisation scripts via mendelian_randomisation.smk which calls TwoSampleMR.R to perform MR

3) HCMR_UKB_MR_Comparison.Rmd compares the MR results for the same trait across the two different datasets UKB and HCMR.