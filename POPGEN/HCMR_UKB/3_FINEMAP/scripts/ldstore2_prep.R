#LDSTORE2 - Prep script to prepare the files needed for LDSTORE2
#Author<- Jonathan Chan
#Date<- 2024-04-05

###--------------------------------------------------------------------------------
library(tidyverse)

args = commandArgs(trailingOnly=TRUE) #Allows input of arguments from Rscript

# test if there is at least one argument<- if not, return an error
if (length(args)==0) {
  stop("At least one argument (i.e the input vep.txt files) must be supplied", call.=FALSE)
}

pheno <- args[1]
chr <- as.numeric(unlist(str_split(args[2],',')))
base_start <- as.numeric(unlist(str_split(args[3],',')))
base_end <- as.numeric(unlist(str_split(args[4],',')))
input_snps <- args[5]
bgen_basepath <- args[6]
output_z_bcor_basepath <- args[7]
incl_file <- args[8]
output_finemap_basepath <- args[9]
master_output <- args[10]
sample_filepath <- args[11]
n_samples <- args[12]

if(str_detect(master_output, 'TRUE')){ #This is to make it play nice with Snakemake by separating the same script into two rules
  master_output <- T
  z_output <- F
} else{
  master_output <- F
  z_output<- T
}



#Local test code - UKB
# pheno <- 'EDN1'
# chr <- c(4,6,6)
# base_start <- c(186678787,11796255,31687637)
# base_end <- c(187678787,12796255,32687637)
# input_snps <- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/formatted/EDN1_manhattan_rsid.tsv'
# bgen_basepath <- '/well/ukbb-wtchg/v3/imputation/ukb_imp_chr'
# output_z_bcor_basepath <- '../pipeline_files/ukbpp_nonHarper_includeHCM/EDN1/EDN1_dataset'
# incl_file <- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/ukb/REGENIE/2_inclHCMcases/regenie_ukb_sampleinclusion_nonNAphenocovars.incl'
# output_finemap_basepath <- '../output/ukbpp_nonHarper_includeHCM/EDN1/EDN1_dataset'
# master_output <- TRUE
# sample_filepath <- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/ukb/ukb11223_imp_chr22_v3_s487324.sample'
# n_samples <- 48128

#Local test code - HCMR
# pheno<- 'cicp'
# chr<- '1,1,2,2,3,4,4,4,7,9,9,10,12,12,13,14,17,18'
# base_start<- '34659480,245380042,1866621,50821126,140813441,102399325,157430099,165025813,123091577,18142278,71033367,36215458,95482784,111358388,40678521,23980087,60220058,67666365'
# base_end<- '35659480,246380042,2866621,51821126,141813441,103399325,158430099,166025813,124091577,19142278,72033367,37215458,96482784,112358388,41678521,24980087,61220058,68666365'
# chr <- as.numeric(unlist(str_split(chr,',')))
# base_start <- as.numeric(unlist(str_split(base_start,',')))
# base_end <- as.numeric(unlist(str_split(base_end,',')))
# bgen_basepath<- '/well/PROCARDIS/agoel/hcm/wallthkmax/bgen/hcmr.p3.hrc.chr'
# incl_file<- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/hcmr/REGENIE/regenie_hcmr_sampleinclusion_nonNAphenocovars.incl'
# sample_filepath<- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/hcmr/hcmr_htscores_v2_essentialcols.sample'
# manhattan_rsid_folder<- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'
# n_samples<- 2465.0
# basepath<- 'hcmr'
#input_snps <- '../..//2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/cicp_manhattan_rsid2.tsv'
#output_z_bcor_basepath <- '../pipeline_files/hcmr/cicp/cicp_dataset'
# output_finemap_basepath <- '../output/hcmr/cicp/cicp_dataset'

#--------------------------------------------------------------------------------
#Rearrange the input_snps to the dataset.z file format here for FINEMAP and LDSTORE2 and filter to SNPs of interest
#1) Convert rsid column to bgenvarID
#3) Filter for variants in genomic regions of interest
#2) Write out as a space-delimited text file

#Accessory function
z_outputter <- function(input_snps, chr, base_start, base_end, output_z){ #Region number is 0-indexed so first region is 0

  
  input_snps<- input_snps%>%
    mutate(chromosome=as.numeric(chromosome))%>%
    filter(chromosome %in% chr) %>%
    filter(position %in% seq(base_start, base_end))
  
  if(str_detect(bgen_basepath,'hcmr')){
    input_snps <- input_snps %>%
      select(rsid=snid,chromosome,position,allele1=allele_A,allele2=allele_B,maf,beta,se=beta_se) 
    
  } else if (str_detect(bgen_basepath,'ukb')){
    input_snps <- input_snps %>%
      #mutate(full_snid=paste(snid, '_',allele_A,'_',allele_B,sep='')) %>%
      select(rsid,chromosome,position,allele1=allele_A,allele2=allele_B,maf,beta,se=beta_se) %>%# Use SNID
      mutate(chromosome = case_when(chromosome %in% seq(1,9) ~ paste('0',chromosome, sep=''),
                                    chromosome %in% seq(10,22) ~ as.character(chromosome))
             )#Add the 0 in front of chromosome as single digit for UKB BGEN files
  }
  
  output_snps <- input_snps %>%
    filter(!is.na(beta),!is.na(se))
  
  write_delim(output_snps,output_z,delim=' ')
  return(output_snps)
}

#Start of main code
#Either output the z files
output_z_names <- str_c(output_z_bcor_basepath, seq_along(chr)-1,'.z')
if(isTRUE(z_output)){
  input_snp_tb <- read_tsv(input_snps, col_types = 'cccnccnnnnnnc')
  
  pwalk(list(chr, base_start, base_end, output_z_names),~z_outputter(..1,..2,..3,output_z=..4,input_snps=input_snp_tb))
}

#Or output the master files7
if(isTRUE(master_output)){
  #--------------------------------------------------------------------------------
  #Define the master file here for LDSTORE2
  
  if(str_detect(bgen_basepath,'hcmr')){
    z <- output_z_names
    append <- '.bgen'
    n_samples <- rep(n_samples, length(z))
    
  } else if (str_detect(bgen_basepath,'ukb')){
    z <- output_z_names
    append <- '_v3.bgen'
    sample <- sample_filepath
    n_samples <- rep(n_samples, length(z)) #This is without HCM cases
  }
  
  bgen <- str_c(bgen_basepath,chr,append)
  
  if(str_detect(bgen_basepath,'hcmr')){
    sample <- str_replace(bgen,'.bgen$','.sample')
  }
  
  bgi <- str_c(bgen,'.bgi')
  bcor <- str_c(output_z_bcor_basepath, seq_along(chr)-1,'.bcor')
  ld <- str_c(output_z_bcor_basepath, seq_along(chr)-1,'.ld')
  

  incl_file <- rep(incl_file, length(z))
  
  output_master_ldstore2 <- tribble(
    ~z, ~bgen, ~bgi, ~bcor, ~ld, ~n_samples,~sample,~incl,
    z,bgen,bgi,bcor,ld, n_samples, sample,incl_file
  ) %>% unnest(cols = c(z, bgen, bgi, bcor, ld, n_samples, sample,incl))
  
  write_delim(output_master_ldstore2,str_c(output_z_bcor_basepath,'_master.txt'),delim=';' )
  
  #--------------------------------------------------------------------------------
  #Define the master file here for FINEMAP
  
  z <- output_z_names
  bcor <- str_c(output_z_bcor_basepath, seq_along(chr)-1,'.bcor')
  snp <- str_c(output_finemap_basepath, seq_along(chr)-1,'.snp')
  config<- str_c(output_finemap_basepath, seq_along(chr)-1,'.config')
  cred<- str_c(output_finemap_basepath, seq_along(chr)-1,'.cred')
  log <- str_c(output_finemap_basepath, seq_along(chr)-1,'.log')
  
  output_master_finemap <- tribble(
    ~z, ~bcor, ~snp, ~config, ~cred, ~n_samples,~log,
    z,bcor,snp,config,cred,n_samples,log
  ) %>% unnest(cols = c(z,bcor,snp,config,cred,n_samples,log))
  
  write_delim(output_master_finemap,str_c(output_z_bcor_basepath,'_master2.txt'),delim=';' )
}



