#FINEMAP2 Credible Set - Analysis Script
#Author: Jonathan Chan
#Date: 2024-04-08

###--------------------------------------------------------------------------------
library(tidyverse)

args = commandArgs(trailingOnly=TRUE) #Allows input of arguments from Rscript

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument (i.e the input vep.txt files) must be supplied", call.=FALSE)
}

pheno <- args[1]
output_finemap_basepath <- args[2]
original_manhattan_tsv <- args[3]
output_filepath <- args[4]

#Local test code
# pheno <- 'NTproBNP'
# output_finemap_basepath <- '../EDA_HCMR/popgen/4_finemap/output/NTproBNP/NTproBNP_dataset'

#Import--------------------------------------------------------------------------
#This imports in the .config files for each genomic region
#There are multiple .cred files corresponding to the credible set given k=1,2,3... (i.e number of credible SNPs)

#Importing in the .config file for each genomic region

config_importer <- function(pheno, input_basepath){
  config_files <- list.files(dirname(input_basepath))
  config_files <- config_files[str_detect(config_files,'config')]
  
  number_of_genomic_regions <- str_match(config_files,'dataset(\\d{1,3})\\.')[,2] %>%
    as.numeric() %>%
    max()
  number_of_genomic_regions <- number_of_genomic_regions +1 #Due to zero-based indexing of the filename
  
  print(str_c('There are ', number_of_genomic_regions, ' genomic regions of interest for this phenotype of ', pheno))
  config_file_list <- map(config_files, ~read_delim(str_c(dirname(input_basepath),'/',.x)))
  
  return(config_file_list)
  
}

configs <- config_importer(pheno, output_finemap_basepath)

#Remap to rsID and extract top 5-------------------------------------------------
#This remapping to rsID is only necessary for the HCMR data
#The UKB BGEN files use RSID already

finemap_snpfiles <- list.files(dirname(output_finemap_basepath))
finemap_snpfiles <- finemap_snpfiles[str_detect(finemap_snpfiles,'\\.snp')]

finemap_snpfile_tb <- map(finemap_snpfiles, ~read_delim(str_c(dirname(output_finemap_basepath),'/',.x),col_types = 'ncnnccnnnnnnnnnn')) %>%
  bind_rows()


config_remapper_extractor <- function(config_files, original_manhattan_tsv, finemap_snpfile_tb, top=5, remap=T){
  
  for (i in seq_along(config_files)){
    config_files[[i]] <- config_files[[i]] %>%
      filter(rank <= top) %>% #Extract top x
      mutate(genomic_region=i) %>%
      mutate(config_vector = str_split(config,','))
  }
  
  config_tb <- config_files %>% bind_rows()
  
  all_rsids_credible <- unlist(config_tb$config_vector) %>% unique() #This is derived from the .config file from FINEMAP
  
  finemap_snps_of_interest <- filter(finemap_snpfile_tb, rsid %in% all_rsids_credible) %>% #This is derived from the .snp file outputted from FINEMAP
    filter(log10bf>0) #Filter only for variants which are likely to be causal (this is supposed to clarify in the case where there are two different allelic changes for the same SNID)
  
  if(isTRUE(remap)){
    input_snp_tb <- read_tsv(original_manhattan_tsv, col_types = 'cccnccnnnnnnc')
    
    #Now remap = mapping back the original SNID to the rsID
    finemap_snps_of_interest <- finemap_snps_of_interest %>% 
      mutate(full_snid = str_c(rsid,'_',allele1,'_',allele2)) %>%
      dplyr::rename(snid=rsid) %>%
      left_join(select(input_snp_tb,full_snid,rsid),by='full_snid') #This gets all the SNPs of interest in their RSID form using the chr:position and AlleleA/B

    rsid <- config_tb$config_vector
    
    for (i in seq_along(config_tb$config_vector)){
      for (j in seq_along(config_tb$config_vector[[i]])){
        rsid[[i]][[j]] <- filter(finemap_snps_of_interest,snid == config_tb$config_vector[[i]][[j]])$rsid
      }
    }
    
    config_tb <- mutate(config_tb, rsid=rsid)
  } else{
    config_tb <- config_tb %>%
      rowwise() %>%
      mutate(rsid=str_c(config_vector, collapse=', '))
  }

  
  config_tb <- config_tb %>%
    select(genomic_region, rank, rsid, everything()) %>%
    rowwise() %>%
    mutate(rsid=str_c(rsid, collapse=', ')) %>%
    select(-config_vector)
  
  return(config_tb)
}

if(str_detect(original_manhattan_tsv,'ukb')){
  remap_indicator <- F
} else {
  remap_indicator <- T
}

output_config <- config_remapper_extractor(configs, original_manhattan_tsv, finemap_snpfile_tb, remap=remap_indicator)

write_tsv(output_config,str_c(dirname(dirname(output_finemap_basepath)),'/',pheno,'_dataset_overallconfig.tsv'))
print(str_c('Written out the overallconfig.tsv to',str_c(dirname(dirname(output_finemap_basepath)),'/',pheno,'_dataset_overallconfig.tsv')))
#End of script


