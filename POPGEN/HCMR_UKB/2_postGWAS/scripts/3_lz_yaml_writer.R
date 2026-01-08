#Script to output .yaml file for input into FINEMAP as config.yaml
#Author: Jonathan Chan
#Date: 2024-07-31

library(stringr)
library(readr)
library(haven)
library(dplyr)
library(purrr)
library(yaml)

args <- commandArgs(trailingOnly=TRUE) #Allows taking of arguments in bash command line #By default you should pass the path (relative to the script)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
chunk_file <- args

lz_params_tb <- args[1]
output_yaml_path <- args[2]
bgen_filepath <- args[3]
incl_filepath <- args[4]
sample_filepath <- args[5]
manhattan_rsid_folder <- args[6]
n_samples <- args[7]
basepath <- args[8]

pheno <- str_match(lz_params_tb, '/([^/]+)_lz_params.tsv$')[,2] #This assumes the format of output/phenoname_manhattan

yaml_writer <- function(lz_params_tb, pheno, output_yaml_path, bgen_filepath, incl_filepath, sample_filepath, manhattan_rsid_folder, n_samples,basepath, read=T){
  
  if(isTRUE(read)){
    lz_params_tb <- read_tsv(lz_params_tb)
    
    print('Imported lz_params_tb')
  }
  
  #Output a .yaml file equivalent for input as config.yaml for the FINEMAP pipeline
  
  output_yaml <- as.yaml(list(
    pheno = as.character(pheno),
    chr = str_c(lz_params_tb$chr, collapse=','),
    base_start = str_c(lz_params_tb$start, collapse=','),
    base_end = str_c(lz_params_tb$end, collapse = ','),
    bgen_basepath = as.character(bgen_filepath),
    incl_file = as.character(incl_filepath),
    sample_filepath=as.character(sample_filepath),
    manhattan_rsid_folder=as.character(manhattan_rsid_folder),
    n_samples=n_samples,
    basepath=as.character(basepath)
  ))
  
  print('Writing out the config.yaml')
  writeLines(output_yaml, str_c(output_yaml_path, pheno,'_config.yaml'))
}


yaml_writer(lz_params_tb, pheno, output_yaml_path, bgen_filepath, incl_filepath, sample_filepath, manhattan_rsid_folder,n_samples,basepath)

message("Script completed successfully.")