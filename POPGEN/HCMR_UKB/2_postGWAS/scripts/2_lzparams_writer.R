## Script to plot a Manhattan plot 
## Author: Jonathan Chan
## Date: 2024-05-06

#----------------------------------------------------------------------
# Setup
library(ggplot2)
library(stringr)
library(readr)
library(haven)
library(dplyr)
library(ggrepel)
library(purrr)
library(yaml)

args <- commandArgs(trailingOnly=TRUE) #Allows taking of arguments in bash command line #By default you should pass the path (relative to the script)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
chunk_file <- args

input <- args[1]
input_lz_tsv <- args[2]
output_path <- args[3]
pvalcutoff <- args[4]
pvalue_format <- args[5]
online_formatting <- args[6]

theme_set(
  theme_classic()
)

if(online_formatting == 'True'){
  online_formatting <- T
} else{
  online_formatting <- F
}

#Local test code
# input <- 'popgen/2_GWAS/output/cicp_bap_manhattan_rsid.tsv'
# input_lz_tsv <- 'popgen/2_GWAS/output/cicp_lz.tsv'
# output_path <- 'popgen/2_GWAS/output/pipeline_files/'
# output_path <- 'popgen/2_GWAS/output/plots/manhattan/'

#-------------------------------------------------------------------------------
#Import
lz_params_writer <- function(pheno, input, input_lz_tsv, output_path, pvalcutoff, pvalue_format, online_formatting, min_points=2){
  
  import <- read_tsv(input, progress=T)
  
  if(pvalue_format == 'log'){
    import <- import %>%
      mutate(pval = 10^-pval)
  } 
  
  if(nrow(filter(import, pval<=5e-8)<20)){
    #Label the SNP if the pvalue < 5e^-8 and less than 20 GWS SNPs
    import <- import %>%
      mutate(Label = ifelse(pval<= 5e-8, rsid, ''))
  } else { #If more than 20 GWS SNPs just label the most significant SNP
    import <- import %>%
      mutate(Label = ifelse(pval==min(import$pval,na.rm=T), rsid, ''))
  }
  
  #Filter out NA values for pval due to failing QC filters in REGENIE
  
  import <- filter(import, !is.na(pval))
  
  #--------------------------------------------------------------------------------
  # Identification of peaks for LocusZoom plotting
  #ID peaks that are meet certain criteria for LocusZoom plotting
  p_cutoff <- pvalcutoff
  span <- 1e+6 #this is the span on each side of the main SNP
  min_points <- min_points
  
  peak_snps <- import %>%
    # Filter out SNPs with p-value > 1e-5
    filter(pval <= p_cutoff) %>%
    # Sort by chromosome, position, and p-value
    arrange(chromosome, position, pval) %>%
    # Group by chromosome
    group_by(chromosome) %>%
    unique()
  
  group <- 1
  
  #This applies the grouping criterion to ID different groups in the SNPs
  print('Applying grouping criteria to SNPs')
  for (i in seq(nrow(peak_snps))){
    if(i ==1){
      peak_snps$group[[i]] <- group #Assign group 1 for the first SNP
    } else{
      current_pos <- peak_snps$position[[i]]
      current_chr <- peak_snps$chromosome[[i]]
      prev_pos <- peak_snps$position[[i-1]]
      prev_chr <- peak_snps$chromosome[[i-1]]
      
      if((current_pos-prev_pos) <= span && current_chr == prev_chr){ #If same chromosome and proximal, then same group
        peak_snps$group[[i]] <- group
        
      } else if (current_pos-prev_pos > span && current_chr == prev_chr){ #If same  chromosome and span distal, then different group
        peak_snps$group[[i]] <- group+1
        group <- group+1
      } else if (current_pos-prev_pos <= span && current_chr != prev_chr){ #If different chromosome and span proximal, then different group
        peak_snps$group[[i]] <- group+1
        group <- group+1
      } else if (current_pos-prev_pos > span && current_chr != prev_chr){#If different chromosome and span distal, then different group
        peak_snps$group[[i]] <- group+1
        group <- group+1
      }
    }
  }
  
  print('Identifying the most significant peak in each group')
  #ID the most significant peak in each group
  peak_snps <- peak_snps %>%
    ungroup() %>%
    group_by(group) %>%
    # Filter groups with at least 2 SNPs
    filter(n() >= min_points) %>%
    # Select the SNP with the minimum p-value in each group
    filter(pval == min(pval)) %>%
    mutate(several_SNPs = ifelse(n() > 1, T,F)) %>% #When there are identical p-values for several SNPs in a group, take the most centrally located one
    ungroup()
  
  peak_snps_multiple <- filter(peak_snps, several_SNPs==TRUE) %>%
    group_by(group) %>%
    mutate(mean_position = mean(position)) %>%
    mutate(abs_distance = abs(position-mean_position))%>%
    filter(abs_distance == min(abs_distance)) %>%
    group_by(group) %>%
    filter(rank(position)==1) %>% #This corrects for cases where there are multiple SNPs with the smallest difference to average by taking the one with the earliest position of those tiebreaker SNPs
    select(-mean_position, -abs_distance)
  
  peak_snps <- filter(peak_snps, several_SNPs ==F) %>%
    bind_rows(peak_snps_multiple) %>%
    arrange(chromosome,position,pval)
  
  rm(peak_snps_multiple)
  
  #Output a tibble which contains the parameter to run LocusZoom standalone to plot LZ plots output
  lz_params_tb <- peak_snps %>%
    mutate(start=position-5e+5, end=position+5e+5,
           metal=input_lz_tsv)
  
  if(isTRUE(online_formatting)){
    lz_params_tb <- lz_params_tb %>%
      select(rsid, chr=chromosome, start, end, allele_A, allele_B, metal, pos=position, pval, beta, beta_se, maf) 
  } else{
    lz_params_tb <- lz_params_tb %>%
      select(rsid, chr=chromosome, start, end, metal, pos=position, pval, beta, beta_se, maf) 
      
  }
  
  lz_params_tb <- lz_params_tb %>%
    mutate(chr=as.integer(chr))
  
  #If >10 individual LZ SNPs, filter to the top 10. The script includes all GWS loci as well!
  if(nrow(lz_params_tb >10)){
    lz_params_tb <- lz_params_tb %>%
      mutate(rank = rank(pval)) %>%
      filter(rank <= 10 | pval <= 5e-8)
  } else if (nrow(lz_params_tb)==0){
    print('No SNPs in the lz_params_tb output so no file written')
  }
  
  if(isTRUE(online_formatting)){ #For online submission to LocusZoom you need to include all the SNPs which are also near your peak SNP 
    lz_params_tb <- lz_params_tb %>%
      select(-metal)
    
    proximal_SNPs <- vector('list', nrow(lz_params_tb))
    
    for(i in seq(1, nrow(lz_params_tb))){
      proximal_SNPs[[i]] <- import %>%
        filter(chromosome == lz_params_tb$chr[i],
               position >= lz_params_tb$pos[i]-span,
               position <= lz_params_tb$pos[i]+span) %>%
        select(rsid, chr=chromosome, allele_A, allele_B, pos=position, pval, beta, beta_se, maf) %>%
        mutate(chr=as.integer(chr))
    }
    proximal_SNPs_tb <- bind_rows(proximal_SNPs)
    rm(proximal_SNPs)
    
    #Join the lz_params_tb to the proximal_SNPs_tb
    lz_params_tb <- lz_params_tb %>%
      bind_rows(proximal_SNPs_tb) %>%
      arrange(chr, pos)
    
  }
  
  print('Writing the lz_params.tsv')
  
  write_tsv(lz_params_tb,str_c(output_path,pheno,'_lz_params.tsv'))
  
  return(lz_params_tb)
}

safe_lzparams_writer <- safely(lz_params_writer) #In case of error due to lack of variants passing the filter
lz_param_tb <- safe_lzparams_writer(input, input_lz_tsv, output_path, pvalcutoff, pvalue_format, min_points)

# message("Script completed successfully.")

#Local Run code via RStudioServer

#ForUKB
# phenos <- c('NTproBNP','MMP1','LGALS3','IL1RL1','TIMP1')
# input <- str_c('popgen/2_gwas/output/gwas/ukb/REGENIE/step2/formatted/',phenos,'_manhattan_rsid.tsv')
# input_lz_tsv <- str_c('/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/ukb/REGENIE/step2/formatted/',phenos,'_lz.tsv')
# output_path <- 'popgen/2_gwas/output/pipeline_files/nonHarperUKB_selectpp/'
# pvalcutoff <- 0.00000005

#walk2(input, input_lz_tsv, ~lz_params_writer(.x,.y,output_path,pvalcutoff))


#ForUKB
# phenos <- c('MAMDC2','LTBP2','FABP3','SHISA5')
# input <- str_c('./popgen/2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part2/formatted/',phenos,'_manhattan_rsid.tsv')
# input_lz_tsv <- str_c('/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part2/formatted/',phenos,'_lz.tsv')
# output_path <- 'popgen/2_gwas/output/pipeline_files/nonHarperUKB_selectpp/2_inclHCMcases/lz_params/'
# pvalcutoff <- 0.00000005
# 
# pwalk(list(phenos,input, input_lz_tsv), ~lz_params_writer(..1,..2,..3,output_path,pvalcutoff))

#ForUKB
# phenos <- c('ANGPT2', 'STC2')
# input <- str_c('./popgen/2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part3/formatted/',phenos,'_manhattan_rsid.tsv')
# input_lz_tsv <- str_c('/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part3/formatted/',phenos,'_lz.tsv')
# output_path <- 'popgen/2_gwas/output/pipeline_files/nonHarperUKB_selectpp/2_inclHCMcases/lz_params/'
# pvalcutoff <- 0.00000005
# 
# pwalk(list(phenos,input, input_lz_tsv), ~lz_params_writer(..1,..2,..3,output_path,pvalcutoff))

#ForHCMR
# phenos <- c('NTproBNP','TnTStat','gal3','st2','timp1','mmp1','cicp')
# input <- str_c('popgen/2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/',phenos,'_manhattan_rsid2.tsv')
# input_lz_tsv <- str_c('/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/',phenos,'_lz.tsv')
# output_path <- 'popgen/2_gwas/output/pipeline_files/hcmr/lz_params/'
# pvalcutoff <- 0.00001
# output_yaml_path <- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/4_finemap/input/configs/hcmr/'
# bgen_filepath <- '/well/PROCARDIS/agoel/hcm/wallthkmax/bgen/hcmr.p3.hrc.chr'
# incl_filepath<- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/hcmr/REGENIE/regenie_hcmr_sampleinclusion_nonNAphenocovars.incl'
# sample_filepath<- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/data/hcmr/hcmr_htscores_v2_essentialcols.sample'
# n_samples <- 2465 #Number of samples in the incl filepath
# basepath <- 'hcmr'
# manhattan_rsid_folder <- '/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'
# 
# lz_param_tb_paths <- str_c('./popgen/2_gwas/output/pipeline_files/hcmr/lz_params/',phenos,'_lz_params.tsv')
# 
# lz_param_tbs <- pmap(list(phenos,input, input_lz_tsv), ~lz_params_writer(..1,..2,..3,output_path,pvalcutoff))
# walk2(lz_param_tb_paths,phenos, ~yaml_writer(.x,.y, output_yaml_path, bgen_filepath, incl_filepath, sample_filepath, manhattan_rsid_folder,n_samples,basepath, read=T))

#For AoUS
# phenos <- c('hcm')
# input <- str_c('popgen/2_gwas/output/gwas/aous/formatted/',phenos,'_manhattan_rsid.tsv')
# input_lz_tsv <- str_c('/well/PROCARDIS/jchan/hcmr_ukbb/popgen/2_gwas/output/gwas/aous/formatted/',phenos,'_lz.tsv')
# output_path <- 'popgen/2_gwas/output/pipeline_files/aous/lz_params/'
# pvalcutoff <- 10^-5
# pvalue_format <- 'log'
# min_points <- 2
# online_formatting <- F #Needed because offline LZ doesn't work for hg38 rsIDs
# 
# pwalk(list(phenos,input, input_lz_tsv), ~lz_params_writer(..1,..2,..3,output_path,pvalcutoff, pvalue_format,online_formatting, min_points))
