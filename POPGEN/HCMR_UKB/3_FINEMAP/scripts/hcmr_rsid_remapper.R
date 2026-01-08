#Remap script from HCMR BGEN ID (i.e chr:position in b37) to rsID 
#Author: Jonathan Chan
#Date: 2024-07-03

###--------------------------------------------------------------------------------
library(rlang)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE) #Allows input of arguments from Rscript

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument (i.e the input vep.txt files) must be supplied", call.=FALSE)
}

pheno <- args[1]
output_finemap_basepath <- args[2]
original_manhattan_rsid2_tsv_basepath <- args[3]

#Local test code on RESCOMP
# pheno <- 'TnTStat'
# output_finemap_basepath <- '../output/hcmr/'
# original_manhattan_rsid2_tsv_basepath <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'

remapper <- function(pheno, output_finemap_basepath, original_manhattan_rsid2_tsv_basepath){
  
  input_overallconfig_tsv <- read_tsv(str_c(dirname(dirname(output_finemap_basepath)),'/',pheno,'_dataset_overallconfig.tsv'), col_types=c('nncc?????????'))
  
#Max number of SNPs in a single causal configuration for this phenotype's overallconfig.tsv
 max_num <- max(str_count(input_overallconfig_tsv$rsid, ','), na.rm=T) + 1
 
 if (is.infinite(max_num)){
   max_num <- 1 #i.e 1 SNP
 }
 
 if(max_num >1){
   split_config <- input_overallconfig_tsv %>%
     separate(rsid, into=str_c('snid',seq(1:max_num)),sep=', ')
 } else {
   split_config <- input_overallconfig_tsv %>%
     dplyr::rename(snid=rsid)
 }
  
#Remap each of the rsid columns to the original manhattan_rsid2.tsv
  rsid2_tsv <- read_tsv(str_c(original_manhattan_rsid2_tsv_basepath,pheno,'_manhattan_rsid2.tsv'),col_types=c('ccnnccnnnnnnnc'))
  
  remap_mini_function <- function(input_tb, rsid2_tsv, snid_name, rsid_name) {
    output_tb <- input_tb %>%
      left_join(select(rsid2_tsv, rsid, snid, snid_allele_A_allele_B), by = setNames("snid", snid_name)) %>%
      select(-!!sym(snid_name)) %>%
      rename(!!sym(rsid_name) := rsid, !!sym(snid_name) := snid_allele_A_allele_B)
    
    return(output_tb)
  }
  
  if(max_num >1){
      remapped <- split_config
      for (i in 1:max_num) {
        remapped <- remap_mini_function(remapped, rsid2_tsv, paste0("snid", i), paste0("rsid", i))
      }
  } else {
    remapped <- remap_mini_function(split_config, rsid2_tsv, 'snid','rsid')
  }

  # Combine all rsid columns back into a single column
  remapped_combined <- remapped %>%
    unite(col = "rsid", starts_with("rsid"), sep = ", ", na.rm = TRUE) %>%
    mutate(rsid = str_trim(rsid)) %>%  # Remove any leading/trailing whitespace
    unite(col = "snid", starts_with("snid"), sep = ", ", na.rm = TRUE) %>%
    mutate(snid = str_trim(snid)) %>%  # Remove any leading/trailing whitespace
    select(genomic_region, rank, rsid, snid, everything())  # Move rsid column to the front
  
  write_tsv(remapped_combined, str_c(dirname(dirname(output_finemap_basepath)),'/',pheno,'_dataset_overallconfig_remapped.tsv'))
  print(str_c('Remapped for ', pheno))

}

remapper(pheno, output_finemap_basepath, original_manhattan_rsid2_tsv_basepath)

#Run locally for all pp in RStudioServer
# pp <- c('NTproBNP','mmp1', 'st2', 'gal3','timp1','TnTStat', 'cicp')
# walk(pp,
#      ~remapper(
#        .,'../output/hcmr/','../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'))

# pp <- c('cicp')
# walk(pp,
#      ~remapper(
#        .,'../output/hcmr/','../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'))

