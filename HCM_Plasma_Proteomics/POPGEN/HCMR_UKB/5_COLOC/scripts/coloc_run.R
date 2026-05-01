#Script to run colocalisation between trait A and trait B via coloc package
# Author: Jonathan Chan
# Date: 2024-10-03

## Inputs:
#GWAS Summary statistics of trait A
#GWAS Summary statistics of trait B
#Genomic region of interest

library(tidyverse)
library(coloc)
library(tictoc)
# library(TwoSampleMR) #Imported for using to harmonise SNPs

args <- commandArgs(trailingOnly=T)

if(length(args)==0){
  print('No arguments defined')
  stop()
}

pathA <- args[1]
pathB <- args[2]

chr <- args[3]
lead_snp_position <- as.numeric(args[4])
rsid <- args[5]
locusname <- args[6]

traittype_A <- args[7]
traittype_B <- args[8]
traitsd_A <- as.numeric(args[9])
traitsd_B <- as.numeric(args[10])

output_path <- args[11]
traitname_A <-args[12]
traitname_B <- args[13]
N_A <- as.numeric(args[14])
N_B <- as.numeric(args[15])
ld_path <- args[16]
susie <- args[17]

if(susie=='TRUE'){
  susie <- T
} else{
  susie <- F
}


tic()
##Test code
# pathA <- '2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/NTproBNP_manhattan_rsid2.tsv'
# pathB <- '5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_rsid.tsv'

# chr <- 7
# lead_snp_position <- 128438284 # FLNC hit from HCM meta analysis = top hit in the per-SNP MR
# rsid <- 'rs66520020' # FLNC hit from HCM meta analysis
# locusname <- 'CCDC136|FLNC'

# chr <- 1
# lead_snp_position <- 16340879 # FLNC hit from HCM meta analysis = top hit in the per-SNP MR
# rsid <- 'rs1048302' # FLNC hit from HCM meta analysis
# locusname <- 'HSPB7'
# 
# traittype_A <- 'quant'
# traittype_B <- 'cc'
# traitsd_A <- 1 #Due to rank-based inverse normalisation in REGENIE for quantitative traits
# traitsd_B <- 5900/74286 #This is standard deviation for trait if quantitative and proportion of samples which are cases if cc
# output_path <- '6_coloc/output/'
# traitname_A <-'NTproBNP'
# traitname_B <- 'HCM'
# N_A <- 2465
# N_B <- 74286
# 
# ld_path <- '/well/PROCARDIS/agoel/ukbb_full/ld/chr__.ukbbv3b.eur.cont.ld'


##Import-------------------------------------------


##Filter---------------------------------------------
#This filters the summary statistics only for SNPs in the genomic region of interest.

gr_filter <- function(summstat,chr_of_interest, region_start, region_end, chr_colname, position_colname ,snps_only=F){
  
  if('rs_number' %in% colnames(summstat)){ #Output from GWAMA
    summstat<- summstat %>% 
      mutate(chromosome = as.numeric(str_match(rs_number, 'chr(\\d{1,2}):')[,2]),
             position = as.numeric(str_match(rs_number, ':(\\d+)_')[,2]))
  }
  out <- summstat %>%
    filter(eval(parse(text=chr_colname))==chr_of_interest) %>%
    filter(eval(parse(text=position_colname)) >= region_start, eval(parse(text=position_colname)) <= region_end)
  
  if('reference_allele' %in% colnames(out)){
    A_colname <- 'other_allele'
    B_colname <- 'reference_allele'
  } else if ('allele_A' %in% colnames(out)){
    A_colname <- 'allele_A'
    B_colname <- 'allele_B'
  }
  
  if(isTRUE(snps_only)){
    print(str_c('Number of indels = ', nrow(filter(str_length(eval(parse(text=A_colname)))!=1, str_length(eval(parse(text=B_colname)))!=1))))
    
    out <- out %>%
      filter(str_length(eval(parse(text=A_colname)))==1, str_length(eval(parse(text=B_colname)))==1)
  }

  return(out)
}

##LD Cleaning/Transformation-------------------------
#This uses bash to parse the LD file to filter only for the genomic region of interest, writing it out to a temp.txt which is imported
#It also checks that the effects are harmonised across the summary statistics and the LD information

ld_filter_importer <- function(ld_path, ld_region_start, ld_region_end, rsid){
  #This uses position to filter for only rows in the LD file corresponding to pairs of SNPs within the region of interest
  #It assumes that there is a column called pos1 and pos2 corresopnding to the SNPs in question and are the 3rd and 6th column
  
  bash_command <- str_glue("awk '{{if ($3 >= {region_start} && $3 <= {region_end} || $6 >= {region_start} && $6 <= {region_end} ) print}}' {ld_path}")
  system2('bash', input=bash_command, stdout='temp.txt')
  
  ld_info <- read.table('temp.txt')
  colnames(ld_info) <- c('snp1','rs1','pos1','snp2','rs2','pos2','r2','d')
  ld_info <<- ld_info #Set as global variable for testing purposes
  system2('bash', input='rm temp.txt')
  
  shared_snps <- intersect(unique(ld_info$snp1), unique(ld_info$snp2))
  
  ld_matrix <- select(ld_info, snp1, snp2, r2) %>%
    unique() %>%
    filter(snp1 %in% shared_snps | snp2 %in% shared_snps) %>%
    arrange(snp1, snp2) 
  
  #Check for those rows which have multiple values of r-sq for the same 2-SNP pair
  multi_rsq <- ld_matrix %>%
    group_by(snp1, snp2) %>% 
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% 
    dplyr::filter(n > 1L) %>%
    unite('snp1_snp2',snp1, snp2)
  
  multi_rsq_vals <- ld_matrix %>% 
    mutate(snp1_snp2 = str_c(snp1, snp2, sep='_')) %>%
    filter(snp1_snp2 %in% multi_rsq$snp1_snp2) %>%
    group_by(snp1_snp2) %>%
    mutate(r2=max(r2, na.rm=T)) %>%
    ungroup() %>%
    unique() %>%
    select(-snp1_snp2)
  
  ld_matrix <- ld_matrix %>%
    mutate(snp1_snp2 = str_c(snp1, snp2, sep='_')) %>%
    filter(!snp1_snp2 %in% multi_rsq$snp1_snp2) %>%
    select(-snp1_snp2) %>%
    bind_rows(multi_rsq_vals)
  
  #Double check that the lead SNP of interest is present within the LD matrix
  if(rsid %in% filter(ld_info, snp1 %in% ld_matrix$snp1)$rs1 | rsid %in% filter(ld_info, snp2 %in% ld_matrix$snp2)$rs2){
    print('Lead SNP found within UKB LD matrix')
  } else{
    print('Lead SNP NOT found within UKB LD matrix')
    stop()
  }
  
  return(ld_matrix)
}

##SNP Harmonisation Across Datasets-------------------------------------------
#This harmonsies the SNPs across the two filtered summary statistics, creating identical column names and making sure the effect alleles are the same (i.e what allele does beta correspond to).

summstat_colname_harmoniser <- function(summstat, traitname){ #The columns of interest are rsid, snptestid, chromosome, position, allele_A, allele_B, beta, beta_se
  
  if (traitname  == 'HCM'){ #i.e summary statistics outputted from GWAMA
    out <- summstat %>%
      mutate(rs_alleleA = str_match(rs_number, '\\d+_([A-Za-z]+)_')[,2],
             rs_alleleB = str_match(rs_number, '\\d+_[A-Za-z]+_([A-Za-z]+$)')[,2]) %>%
      mutate(rs_alleleA2 = case_when(other_allele == rs_alleleA ~ rs_alleleA,
                                    other_allele != rs_alleleA & other_allele == 'I' ~ ifelse(str_length(rs_alleleA) > str_length(rs_alleleB), rs_alleleA, rs_alleleB),
                                    other_allele != rs_alleleA & other_allele == 'D' ~ ifelse(str_length(rs_alleleA) < str_length(rs_alleleB), rs_alleleA, rs_alleleB),
                                    rs_alleleA == reference_allele ~ other_allele)) %>% #In case it's just swapped reference and other allele in the rs_number and is also a SNP
      mutate(rs_alleleB2 = case_when(reference_allele == rs_alleleB ~ rs_alleleB,
                                     reference_allele != rs_alleleB & reference_allele == 'I' ~ ifelse(str_length(rs_alleleB) > str_length(rs_alleleA), rs_alleleB, rs_alleleA),
                                     reference_allele != rs_alleleB & reference_allele == 'D' ~ ifelse(str_length(rs_alleleB) < str_length(rs_alleleA), rs_alleleB, rs_alleleA),
                                     rs_alleleB == other_allele ~ reference_allele))%>% #In case it's just swapped reference and other allele in the rs_number and is also a SNP
      select(rsid, chromosome, position, allele_A=rs_alleleA2, allele_B=rs_alleleB2, beta, beta_se=se, eaf, pval=`p-value`) %>% #No need to inverse beta values because I map everything to the reference_allele and other_allele
      mutate(rsid=ifelse(!str_detect(rsid, 'rs'), str_c(as.character(chromosome),':',as.character(position), '_',allele_A, '_', allele_B), rsid)) #For harmonisation purposes, make the rsID correctly in the order of alleleA_alleleB
    
  } else { #i.e summary statistics outputted from REGENIE
    out <- summstat %>%
      select(rsid, chromosome, position, allele_A, allele_B, beta, beta_se, eaf, pval)
  }
  
  return(out)
}


summstats_harmoniser <- function(summstatA, summstatB){ #This by default harmonises to summstatA
  
  merged_summstats <- left_join(summstatA, summstatB, by=c('rsid')) %>% #By default this harmonises the SNPs in summstatB to match that of summstatA
    filter(!is.na(chromosome.x), !is.na(chromosome.y)) %>%
    filter(position.x==position.y) %>%
    filter(allele_A.x != allele_A.y | allele_B.x != allele_B.y) %>% #Find the mismatches
    mutate(beta.y = ifelse(`allele_A.y` == `allele_B.x` & `allele_B.y` == `allele_A.x`, -beta.y, beta.y),
           allele_A.y2 = ifelse(`allele_A.y` == `allele_B.x` & `allele_B.y` == `allele_A.x`, `allele_A.x`, `allele_A.y`),
           allele_B.y2 = ifelse(`allele_A.y` == `allele_B.x` & `allele_B.y` == `allele_A.x`, `allele_B.x`, `allele_B.y`)
           ) #If in summstatB, alleles A and B are directly swapped relative to summstatA
    
  new_summstatB <- summstatB %>%
    filter(!rsid %in% merged_summstats$rsid) #i.e filter out those rsIDs which don't have mismatch problems with summstatA
  
  merged_summstats <- merged_summstats %>%
    filter(`allele_A.x` == `allele_A.y2`) %>%#Filter for only harmonised SNPs
    select(rsid, chromosome=chromosome.y, position=position.y, allele_A=`allele_A.y2`, allele_B = `allele_B.y2`, beta=beta.y, beta_se=beta_se.y, eaf=eaf.y, pval=pval.y)
    
  new_summstatB <- new_summstatB %>%
    bind_rows(merged_summstats)
  
  #Double check there are no longer harmonisation issues
  merged_summstats <- left_join(summstatA, new_summstatB, by=c('rsid')) %>%
    filter(!is.na(chromosome.x), !is.na(chromosome.y)) %>%
    filter(position.x==position.y) %>%
    filter(allele_A.x != allele_A.y | allele_B.x != allele_B.y) #Find the mismatches
  
  print('Printing the mismatches')
  print(merged_summstats) #Print the mismatches. If they are simply multiallelic, that is fine because you output summstatB separately
  
  return(list(summstatA, new_summstatB))
  
}

##Transform-----------------------------------
#This transforms the summary statistics each into coloc dataset objects.
#It also harmonises the ld_info to make sure that the SNP effects are in the same orientation.
#Coloc dataset object = list of
  ## beta = per-SNP beta
  ## varbeta = per-SNP standard error^2 = variance
  ## snp = ID for each SNP (preferably rsID) but uses snptestID due to multi-allelic SNPs under the same rsID
  ## position = position on the chromosome of interest
  ## type = either quant or cc for trait type
  ## sdY = standard deviation of trait (required for estimation of prior on true effect size for interpretation of scale of beta)
  ## LD = for SuSiE, you need LD information in that SNP region = UKB LD information for EUR ancestry used

summstats_to_coloc_dataset <- function(summstat, ld_matrix, plot_output_path, trait_N, snp_colname='snptestid',beta_colname='beta', betase_colname='beta_se',trait_type='qt', trait_sd=1,filtersnps_to_LDpanel=T, susie=T){
  
  summstat <- summstat %>% #This snptestID should correctly reflect the allele_A and allele_B columns given harmonisation has previously occurred in summstat_colname_harmoniser 
    mutate(snptestid = str_c(chromosome,':',position,'_',allele_A,'_', allele_B))
  #If still duplicates e.g due to the reference/alt allele referring to the same across two 'rsIDs' but not being meta-analysed together due to different rsID names e.g 7:128092841_G_GAGAA 7:128092841_GAGAA_G differnet names but reference_allele = 'I' so not meta-analysed correctly
  #In that case, take the row with smaller standard error
  summstat <- summstat %>%
    group_by(snptestid) %>%
    filter(beta_se == min(beta_se, na.rm=T)) %>%
    ungroup()

  if(isTRUE(susie)){
    
    #Harmonise the summary statistics to have the same snptestID as the LD panel and reverse beta sign if necessary
    ld_check <- c(unique(ld_matrix$snp1), unique(ld_matrix$snp2))
    ld_check_tb <- tibble(
      'snptestid'=ld_check,
      'position' = as.numeric(str_match(ld_check, ':(\\d+)_')[,2]),
      'allele_A' = str_match(ld_check, ':\\d+_([A-Za-z]+)_')[,2],
      'allele_B' = str_match(ld_check, ':\\d+_[A-Za-z]+_([A-Za-z]+)$')[,2],
    )

    
    mismatches <- left_join(summstat, ld_check_tb, by=c('position')) %>%
      filter(!is.na(snptestid.y)) %>% #Filter for only rows which are present in the summary statistic and LD panel
      mutate(snp_flipped = ifelse(`allele_A.x` == `allele_B.y` & `allele_B.x` == `allele_A.y`, T, F)) %>%
      filter(snp_flipped==T) %>%
      mutate(beta=ifelse(`allele_A.x` == `allele_B.y` & `allele_B.x` == `allele_A.y`, -beta, beta),
             allele_A.x2 = ifelse(`allele_A.x` == `allele_B.y` & `allele_B.x` == `allele_A.y`, allele_A.y, allele_A.x), #Flip the SNPs in the summary statistic (not the LD panel) to align them to the LD panel direction
             allele_B.x2 = ifelse(`allele_A.x` == `allele_B.y` & `allele_B.x` == `allele_A.y`, allele_B.y, allele_B.x)) 
    
    print('Checking for if any mismatches did not flip properly in the alignment of SNPs to the LD panel')
    print(nrow(filter(mismatches, `allele_B.y` != `allele_B.x2` | `allele_A.y` != `allele_A.x2`))) #Check if any mismatches didn't flip correctly
    
    summstat <- summstat %>%
      filter(!snptestid %in% mismatches$snptestid.x)
    
    mismatches <- mismatches %>%
      select(rsid, chromosome, position, allele_A=allele_A.x2, allele_B=allele_B.x2, beta, beta_se, eaf, pval)
    
    summstat <- summstat %>%
      bind_rows(mismatches) %>%
      mutate(snptestid = str_c(chromosome,':',position,'_',allele_A,'_', allele_B)) %>% #Recompute snptestid after alignment to LD panel SNPs
      unique()
    
    #Check for duplicates again
    print(length(unique(summstat$snptestid))==length(summstat$snptestid)) #TRUE if no duplicates
    
    #Filter the LD panel to the same as the summary statistics SNPs
    ld_matrix <- ld_matrix %>%
      filter(snp1 %in% summstat[[snp_colname]], snp2 %in% summstat[[snp_colname]]) %>%
      pivot_wider(names_from = snp2, values_from = r2)  %>%
      column_to_rownames("snp1") %>%
      as.matrix()
    
    #Make matrix square
    common_snps <- intersect(rownames(ld_matrix ), colnames(ld_matrix ))
    square_matrix <- ld_matrix[common_snps, common_snps]
    
    #Replace NUll or NA values in the square matrix with 0 because the initial LD matrix is sparse so NULL value if r2 < 0.2
    square_matrix[is.na(square_matrix)] <- 0
    
    #Filter SNPs in summmary statistics to those found in the LD panel
    if(isTRUE(filtersnps_to_LDpanel)){
      summstat <- summstat %>%
        filter(eval(parse(text=snp_colname)) %in% c(colnames(square_matrix), rownames(square_matrix)))
    }
    
    summstat <- summstat %>% unique() #Ensure no duplicated SNPs
    
    #Create coloc dataset object
    coloc_dataset_obj <- list(
      beta=summstat[[beta_colname]],
      varbeta = summstat[[betase_colname]]^2,
      snp=summstat[[snp_colname]],
      position=summstat[['position']],
      type=trait_type,
      sdY=trait_sd,
      N=trait_N,
      LD=square_matrix #The LD information should correspond to the correct allele order because snptestid is used for SNP names in both the LD matrix and the summary statistics
    )
    
  } else{ #i.e if not SuSiE - no need for LD or filtering of summary statistics to match LD panel
    #Create coloc dataset object
    coloc_dataset_obj <- list(
      beta=summstat[[beta_colname]],
      varbeta = summstat[[betase_colname]]^2,
      snp=summstat[[snp_colname]],
      position=summstat[['position']],
      type=trait_type,
      sdY=trait_sd,
      N=trait_N
    )
  }

  if(is_null(check_dataset(coloc_dataset_obj))){
    png(str_c(plot_output_path, '_qc_locusplot.png'))
    plot <- plot_dataset(coloc_dataset_obj)
    dev.off()
    
    if(isTRUE(susie)){
      png(str_c(plot_output_path, '_qc_alignment.png'))
      check_alignment(coloc_dataset_obj) # It compares the product of the Z scores for all SNP pairs against their correlation, and plots a histogram of this ratio for SNPs in some degree of LD.
      #We do not expect to see strongly negative products of Z scores for strongly positive correlations or vice versa. 
      #Therefore this plot should be skewed to have more positive values than negative.
      dev.off()
    }

    return(coloc_dataset_obj)
  } else{
    print('Error with check_dataset() function')
  }
}

##Run SuSiE---------------------------------------
#This runs fine-mapping via SuSiE for each trait individually at that region of interest and save output to .rds file.

susie_function <- function(input_coloc_obj, output_path, max_iterations=100){
  
  susie_obj <- runsusie(input_coloc_obj, maxit=max_iterations, coverage=0.1)
  
  # saveRDS(susie_obj, str_c(output_path,'_susie_obj.rds'))
  
  return(susie_obj)
}

##Main-------------------------------------

main <- function(pathA, pathB, chr, lead_snp_position,traittype_A,traittype_B, traitsd_A, traitsd_B, output_path, traitname_A, traitname_B, N_A, N_B, locusname, rsid, susie=T){
  
  #Define genomic region
  region_start <- lead_snp_position - 500000
  region_end <- lead_snp_position + 500000
  
  ld_region_start <- lead_snp_position - 700000
  ld_region_end <- lead_snp_position + 700000
  
  #Import summary statistics
  summstats <- map(list(pathA, pathB), ~read_tsv(.))
  
  #Filter for summary stats in genomic region
  filtered_summstats <- map(summstats, ~gr_filter(.,chr, region_start, region_end, 'chromosome','position'))
  walk(filtered_summstats, ~print(nrow(.))) #Number of SNPs after filtering
  rm(summstats)
  
  #Harmonise the SNPs across the two datasets
  harmonised_summstats <- map2(filtered_summstats, c(traitname_A, traitname_B), ~summstat_colname_harmoniser(.x, .y ))
  harmonised_summstats <- summstats_harmoniser(harmonised_summstats[[1]], harmonised_summstats[[2]])
  
  if(isTRUE(susie)){
    #Import in LD information from UKB
    ld_path <- str_replace(ld_path, '__',as.character(chr)) #Choose the correct chromosome's LD file
    ld_matrix <- ld_filter_importer(ld_path, ld_region_start, ld_region_end, rsid)
    
    #Convert summary statistic to coloc dataset object + plot some QC plots
    if(!dir.exists(str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname))){dir.create(str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname))}
    
    plot_output_path <- str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname, '/susie_')
    coloc_datasets <- pmap(list(harmonised_summstats, c(N_A,N_B), c(traittype_A, traittype_B), c(traitsd_A, traitsd_B), c('hcmr','hcmMA')), 
                           ~summstats_to_coloc_dataset(..1, trait_N=..2,  trait_type=..3, trait_sd=..4, plot_output_path=str_c(plot_output_path,..5), ld_matrix=ld_matrix))
    
    #Run SuSiE for fine-mapping within each trait 
    susie_objs <- map2(coloc_datasets, list('hcmr','hcmMA'), ~susie_function(.x, str_c(susie_output_path,.y), max_iterations=1000))
    
    #Run coloc-SuSiE
    coloc_results <- coloc.susie(susie_objs[[1]], susie_objs[[2]])
    
  } else{
    #Convert summary statistic to coloc dataset object + plot some QC plots
    if(!dir.exists(str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname))){dir.create(str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname))}
    plot_output_path <- str_c(output_path,'/plots/', traitname_A, '_',traitname_B,'/', locusname, '/abf_')
    coloc_datasets <- pmap(list(harmonised_summstats, c(N_A,N_B), c(traittype_A, traittype_B), c(traitsd_A, traitsd_B), c('hcmr','hcmMA')), 
                           ~summstats_to_coloc_dataset(..1, trait_N=..2, trait_type=..3, trait_sd=..4, plot_output_path=str_c(plot_output_path,..5), ld_matrix=ld_matrix, susie=F))
    
    #Run coloc-abf i.e enumeration method
    coloc_results <- coloc.abf(coloc_datasets[[1]], coloc_datasets[[2]], p12=5e-6)
    
    #Run sensitivity analyses showing how prior and posterior pribabilities change as a function of p12 (i.e P(H4) that shared causal variants across traits)
    png(str_c(plot_output_path, 'results_sensitivity.png'))
    plot <- sensitivity(coloc_results,"H4 > 0.4")
    dev.off()
  }
  
  print(coloc_results)
  
  #Output results as .tsv format

  write_tsv(coloc_results[['results']], str_c(plot_output_path, 'fullresults.tsv'))
  sink(str_c(plot_output_path,rsid,'_',chr,'_',lead_snp_position, '_summaryresults.txt'))
  print(coloc_results$summary)
  sink()
  
  print('Colocalisation finished')
  
  
}

main(pathA, pathB, chr, lead_snp_position,traittype_A,traittype_B, traitsd_A, traitsd_B, output_path, traitname_A, traitname_B, N_A, N_B, locusname, rsid, susie=susie)

toc()
