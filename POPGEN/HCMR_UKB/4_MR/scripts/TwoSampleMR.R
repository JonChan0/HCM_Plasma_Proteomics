#Script to run the Mendelian Randomisation itself via TwoSampleMR and MendelianRandomization
#Author: Jonathan Chan
#Date: 2024-04-30

###--------------------------------------------------------------------------------
library(tidyverse)
library(TwoSampleMR)
library(tictoc)
library(ggmanh)
#source('popgen/5_MR/scripts/sparsePCA_MVMR.R')
# source('./sparsePCA_MVMR.R')

args = commandArgs(trailingOnly=TRUE) #Allows input of arguments from Rscript

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}

exposure_input_path <- args[1]
outcome_input_path <- args[2]
selected_exposure_snps_string <- args[3]
output_path <- args[4]
exposure_name <- args[5]
outcome_name <- args[6]
ld_clumping <- args[7]
ld_eur_bed_file <- args[8]
plink_binary_path <- args[9]
multivariate <- args[10]

if(ld_clumping=='TRUE'){
  ld_clumping <- TRUE
} else if(ld_clumping == 'FALSE' | !exists('ld_clumping')){
  ld_clumping <- FALSE
}

if(multivariate=='TRUE'){
  multivariate <- TRUE
} else if(multivariate == 'FALSE' | !exists('multivariate') ){
  multivariate <- FALSE
}

if(selected_exposure_snps_string !='FALSE'){
  #Assume that the selected exposure SNPs come in as a single string
  selected_exposure_snps <- str_split_1(selected_exposure_snps_string,',')
  
  print(selected_exposure_snps)
}


#Local test code for bidirectional MR for HCMR
#Direction1
# exposure_input_path <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/NTproBNP_manhattan_rsid2_logpval.tsv'
# outcome_input_path <- '../../5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.nohcmr_rsid.tsv'
# selected_exposure_snps <- c('rs198388')
# output_path <- '../..//5_MR/output/hcmr/'
# exposure_name <- 'NTproBNP'
# outcome_name <- 'hcm'
# ld_clumping <-  FALSE
# multivariate <- F

#Direction 2 #Uses 19 GWS SNPs in HCMR-less meta-analysis
# exposure_input_path <- '../../5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.nohcmr_filtered.tsv'
# outcome_input_path <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/NTproBNP_manhattan_rsid2_logpval.tsv'
# selected_exposure_snps <- c("rs1048302","rs3845778","rs4894803","rs2191446","rs11748963","rs3176326","rs12210733","rs66520020","rs2177843","rs11196085","rs17617337","rs182427065","rs41306688","rs8006225","rs8033459","rs28768976","rs7210446","rs2644262","rs5760054")
# output_path <- '../../5_MR/output/hcmr/'
# exposure_name <- 'hcm'
# outcome_name <- 'NTproBNP'
# ld_clumping <- FALSE
# multivariate<- F

#Rescomp test code for bidirectional MR for UKB PP 
# Direction 1
# exposure_input_path <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/formatted/APOM_manhattan_rsid_logpval.tsv'
# outcome_input_path <- '../input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_rsid.tsv'
# selected_exposure_snps <- c('rs805270')
# output_path <- '../../5_MR/output/ukb_noHarper_inclHCM/'
# exposure_name <- 'APOM'
# outcome_name <- 'hcm'
# ld_clumping <-  FALSE
# multivariate <- F

#Direction 2 #Using the full MTAG SNPs
# exposure_input_path <- '../../5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_filtered.tsv'
# outcome_input_path <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/formatted/NTproBNP_manhattan_rsid_logpval.tsv'
# selected_exposure_snps <- c('rs2503715','rs11121483','rs1048302','rs6699769','rs2810883','rs74139614','rs547048','rs11687178','rs3845778','rs1009358','rs2540277','rs6747402','rs295114','rs1699340','rs7612736','rs6796333','rs79502300','rs4894803','rs11730129','rs75096272','rs2191446','rs11748963','rs6914805','rs3176326','rs816379','rs12210733','rs1159974','rs2184370','rs66520020','rs2945236','rs12216858','rs7824244','rs7835298','rs35006907','rs13297385','rs2645210','rs2182400','rs2177843','rs11196085','rs17617337','rs12270374','rs182427065','rs11571171','rs10841524','rs58747679','rs7487962','rs116904997','rs9515201','rs41306688','rs113907726','rs8019721','rs140763915','rs8006225','rs2002854','rs8033459','rs4984741','rs11860560','rs28575249','rs28768976','rs17817857','rs7210446','rs497999','rs2644262','rs6566955','rs12460541','rs62222424','rs5760054','rs139920')
# output_path <- '../../5_MR/output/ukb_noHarper_exclHCM/'
# exposure_name <- 'hcm'
# outcome_name <- 'NTproBNP'
# ld_clumping <- FALSE
# multivariate<- F

#32 GWS SNPs for HCM for HCM -> TNNI3
# exposure_input_path <- '../input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_rsid.tsv'
# outcome_input_path <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/formatted/TNNI3_manhattan_rsid_logpval.tsv'
# selected_exposure_snps <- c("rs1048302","rs11687178","rs3845778","rs2540277","rs6747402","rs7612736","rs4894803","rs2191446","rs11748963","rs3176326","rs12210733","rs66520020","rs7824244","rs35006907","rs2645210","rs2177843","rs11196085","rs17617337","rs12270374","rs182427065","rs7487962","rs41306688","rs113907726","rs8006225","rs8033459","rs28768976","rs7210446","rs2644262","rs6566955","rs12460541","rs62222424","rs5760054")
# output_path <- '../../5_MR/output/ukb_noHarper_inclHCM/'
# exposure_name <- 'hcm'
# outcome_name <- 'TNNI3'
# ld_clumping <-  FALSE
# multivariate <- F

#Triadin locus in HCM -> HRC levels in UKB
# exposure_input_path <- '../../5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_filtered.tsv'
# outcome_input_path <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/formatted/HRC_manhattan_rsid_logpval.tsv'
# selected_exposure_snps <- c('rs9320939')
# output_path <- '../../5_MR/output/ukb_noHarper_inclHCM/'
# exposure_name <- 'hcmTRDN'
# outcome_name <- 'HRC'
# ld_clumping <- FALSE
# multivariate<- F

#Rescomp test code for bidirectional MR for AoUS
# Direction 1
# exposure_input_path <- '../../2_gwas/output/gwas/aous/formatted/NTproBNP_manhattan_rsid_logpval.tsv'
# outcome_input_path <- '../input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_rsid.tsv'
# selected_exposure_snps <- c('rs12406383')
# output_path <- '../../5_MR/output/aous/'
# exposure_name <- 'NTproBNP'
# outcome_name <- 'hcm'
# ld_clumping <-  FALSE
# multivariate <- F

#Direction 2
# exposure_input_path <- '../../5_MR/input/gwas_summary_statistics/HCM/hcm_meta.230523.fix.gwama.noukb_filtered.tsv'
# outcome_input_path <- '../../2_gwas/output/gwas/aous/formatted/NTproBNP_manhattan_rsid_logpval.tsv'
# selected_exposure_snps <- c("rs1048302","rs11687178","rs3845778","rs2540277","rs6747402","rs7612736","rs4894803","rs2191446","rs11748963","rs3176326","rs12210733","rs66520020","rs7824244","rs35006907","rs2645210","rs2177843","rs11196085","rs17617337","rs12270374","rs182427065","rs7487962","rs41306688","rs113907726","rs8006225","rs8033459","rs28768976","rs7210446","rs2644262","rs6566955","rs12460541","rs62222424","rs5760054")
# output_path <- '../../5_MR/output/aous/'
# exposure_name <- 'hcm'
# outcome_name <- 'NTproBNP'
# ld_clumping <- FALSE
# multivariate<- F

#Test code for unidirectional MR for CMR
# exposure_input_path <- 'popgen/2_gwas/output/cicp_manhattan_rsid.tsv'
# outcome_input_path <- '../GeneticBasis_Cardiomyopathies_Review/data/tadros21_gwas/LVM_metal.TBL.gz'
# selected_exposure_snps <- c('rs2109080','rs61594226')
# output_path <- 'popgen/5_MR/output/pp_x_cmr_crosspheno/'
# exposure_name <- 'cicp'
# outcome_name <- 'LVM'
# ld_clumping <-  FALSE

#Local test code for multivariate MR
# exposure_input_path <- 'popgen/5_MR/input/gwas_summary_statistics/tadros21_cmr_gwas/allcmr_noextra_tadros21_gws_MRpackage_loci.tsv'
# outcome_input_path <- 'popgen/2_gwas/output/NTproBNP_manhattan_rsid.tsv'
# output_path <- 'popgen/5_MR/output/pp_x_cmr_crosspheno/multivariate/'
# exposure_name <- 'allcmr_noextra'
# outcome_name <- 'NTproBNP'
# multivariate<- T

#Local test code for univariate MR -> NTproBNP
# exposure_input_path <- 'popgen/5_MR/input/gwas_summary_statistics/tadros21_cmr_gwas/meanLVWT_tadros21_gws_loci.tsv'
# outcome_input_path <- 'popgen/2_gwas/output/NTproBNP_manhattan_rsid.tsv'
# output_path <- 'popgen/5_MR/output/pp_x_cmr_crosspheno/'
# exposure_name <- 'meanLVWT_only'
# selected_exposure_snps <- c('rs12906223','rs4820654')
# outcome_name <- 'NTproBNP'
# multivariate<- F
# ld_clumping <-  FALSE

#---------------------------------------------------------------------------------
#Load the exposure GWAS summary statistics + Extract the instruments
tic()
print('Importing the exposure GWAS summary statistics')

## In the case of non-multvariate MR
if(isFALSE(multivariate)){
  if(str_detect(exposure_input_path, '_manhattan_rsid\\.tsv|_manhattan_rsid2\\.tsv')){
    #phenotype_name <- str_match(exposure_input_path,'/([^_/]+)_manhattan_rsid.tsv')[,2]
    
    full_exposure_summstats <- read_tsv(exposure_input_path, col_types=c('ccnnccnnnnnnc')) %>%
      filter(pval < 10^-5) %>% #Use arbitrary threshold for GWAS significance
      mutate(Phenotype=exposure_name)
    
    # #Plot a manhattan but only taking SNPs iin the same chromosome which have pval < 10^-5 to ensure that they have peaks and are not flat or singletons 
    # leadSNP_ranges <- full_exposure_summstats %>%
    #   filter(rsid %in% selected_exposure_snps) %>%
    #   select(rsid,chromosome, position) 
    # # %>%
    # #   mutate(upper_range = position + 500000,
    # #          lower_range = ifelse(position -500000 <0, 0,position-5000000))
    # 
    # plot_snps <- full_exposure_summstats %>%
    #   filter(chromosome %in% leadSNP_ranges$chromosome)
    # 
    # g <- manhattan_plot(x = mutate(plot_snps, position=as.numeric(position)), pval.colname = "pval", chr.colname = "chromosome", pos.colname = "position", 
    #                     plot.title = str_c('Manhattan plot for selected SNPs with pvalue < 10^-5 for phenotype ',exposure_name), y.label = "-log10(pval)",
    #                     rescale=F, label.colname='rsid')
    
    #This formats it for TwoSampleMR
    exposure_summstats <- format_data(full_exposure_summstats, type='exposure',
                                      snps=selected_exposure_snps,
                                      snp_col='rsid',beta_col='beta',se_col='beta_se',eaf_col='eaf',effect_allele_col='allele_B',other_allele_col='allele_A',
                                      pval_col='pval',chr_col='chromosome',pos_col='position', samplesize_col='all_total')
    
    
  } 
  else if(str_detect(exposure_input_path, '_manhattan_rsid_logpval\\.tsv|_manhattan_rsid2_logpval\\.tsv')){
      
      full_exposure_summstats <- read_tsv(exposure_input_path, col_types=c('ccnnccnnnnnnc')) %>%
        filter(pval >= 5) %>% #Use arbitrary threshold for GWAS significance
        mutate(Phenotype=exposure_name)
      
      #This formats it for TwoSampleMR
      exposure_summstats <- format_data(full_exposure_summstats, type='exposure',
                                        snps=selected_exposure_snps,
                                        snp_col='rsid',beta_col='beta',se_col='beta_se',eaf_col='eaf',effect_allele_col='allele_B',other_allele_col='allele_A',
                                        pval_col='pval',chr_col='chromosome',pos_col='position', samplesize_col='all_total', log_pval = T)
      
  } else if(str_detect(exposure_input_path, 'gwama')){
    
    full_exposure_summstats <- read_tsv(exposure_input_path, col_types=c('ccccnnnnnnnnnnnnncc')) %>%
      filter(`p-value` < 10^-5) %>% #Use arbitrary threshold for GWAS significance
      mutate(Phenotype=exposure_name)
    
    if(!all(str_detect(colnames(full_exposure_summstats),'chr'))){
      full_exposure_summstats <- full_exposure_summstats %>%
        mutate(chr = str_match(snptestid,'(\\d{1,2}):')[,2],
               pos = str_match(snptestid, '\\d{1,2}:(\\d+)_')[,2])
    }
    
    #Plot a manhattan but only taking SNPs iin the same chromosome which have pval < 10^-5 to ensure that they have peaks and are not flat or singletons 
    
    # leadSNP_ranges <- full_exposure_summstats %>%
    #   filter(rsid %in% selected_exposure_snps) %>%
    #   select(rsid,chr, pos) 
    # # %>%
    # #   mutate(upper_range = position + 500000,
    # #          lower_range = ifelse(position -500000 <0, 0,position-5000000))
    # 
    # plot_snps <- full_exposure_summstats %>% #This is from the 
    #   filter(chr %in% leadSNP_ranges$chr)
    # 
    # g <- manhattan_plot(x = mutate(plot_snps, pos=as.numeric(pos)), pval.colname = "p-value", chr.colname = "chr", pos.colname = "pos", 
    #                     plot.title = str_c('Manhattan plot for selected SNPs with pvalue < 10^-5 for phenotype ',exposure_name), y.label = "-log10(pval)",
    #                     rescale=F, label.colname='rsid')
    
    #This formats it for TwoSampleMR
    exposure_summstats <- format_data(full_exposure_summstats, type='exposure', 
                                      snps=selected_exposure_snps,
                                      snp_col='rsid',beta_col='beta',se_col='se',eaf_col='eaf',effect_allele_col='reference_allele',other_allele_col='other_allele',
                                      pval_col='p-value', samplesize_col='n_samples')
    
  } else if(str_detect(exposure_input_path, 'tadros21')){ #This is for univariate Tadros21 CMR -> ... for sensitivity analsyes
    full_exposure_summstats <-read_tsv(exposure_input_path) %>%
      dplyr::rename(SNP=rsid) %>%
      filter(pval < 10^-5) #Assume no need to specify which SNPs to use as exposure and just use all the significant ones because these are GWS from Tadros21 CMR GWAS
    
    #This formats it for TwoSampleMR
    exposure_summstats <- format_data(full_exposure_summstats, type='exposure',
                                      snps=selected_exposure_snps,
                                      snp_col='SNP',beta_col='beta',se_col='se',eaf_col='EAF',effect_allele_col='EA',other_allele_col='NEA',
                                      pval_col='pval',chr_col='chr',pos_col='pos', samplesize_col='samplesize')
    
  }
  toc()
  
  tic()
  #LD clumping to only return independent SNPs
  if (isTRUE(ld_clumping)){
    print('LD clumping the associated SNPs with exposure')
    exposure_clumped <- clump_data(exposure_summstats,
                                   bfile=ld_eur_bed_file,
                                   plink_bin=plink_binary_path) #Uses 1KGP EUR LD Panel
    
    #Plot a Manhattan to determine validity of instruments after LD clumping
    g2 <- manhattan_plot(x = mutate(exposure_clumped, pos.exposure=as.numeric(pos.exposure)), pval.colname = "pval.exposure", chr.colname = "chr.exposure", pos.colname = "pos.exposure", 
                         plot.title = str_c('Manhattan plot for LD clumped SNPs with pvalue < 10^-5 for phenotype ',exposure_name), y.label = "-log10(pval)",
                         rescale=F, label.colname='SNP')
    #print(g)
    ggsave(str_c(output_path,'plots/',exposure_name,'_to_',outcome_name,'/',exposure_name,'_ldclumped_manhattan.png'),g2,dpi=600)
    rm(g2)
    
    if(nrow(exposure_clumped) ==0){
      print('LD clumping removed all SNPs or you ran out of API calls')
      stop()
    }
    
  } else{
    exposure_clumped <- exposure_summstats
  }
  
  toc()
  

  
#------------------------------------------------------------------------------------
#If there is multivariate MR
} else if (isTRUE(multivariate)){
  if(str_detect(exposure_input_path,'tadros21')){
    
    #TwoSampleMR version
      # exposure_input_paths <- str_c(exposure_input_path,list.files(exposure_input_path, pattern='.tsv'))
      # 
      # #This allows the use of SNPs even if they're not associated with ALL exposures (that is the default for TwoSampleMR) 
      # exposure_clumped <- mv_extract_exposures_local_adjusted(exposure_input_paths, sep='\t',
      #                                                  snp_col='rsid',beta_col='beta',se_col='se',eaf_col='EAF',effect_allele_col='EA',other_allele_col='NEA',
      #                                                  pval_col='pval', samplesize_col='samplesize',
      #                                                  plink_bin=plink_binary_path, bfile=ld_eur_bed_file
      #                                                  ) #It automatically performs LD clumping
    
    
    #MendelianRandomization Version
    exposure_clumped <-read_tsv(exposure_input_path) %>%
      dplyr::rename(SNP=rsid)
  }
  
  
}

#-------------------------------------------------------------------------------------
#Continue regardless of multivariate or not
print(head(exposure_clumped))

# Filter for the selected input instruments of choice (after fine-mapping for HCMR GWAS)
# if(!identical(selected_exposure_snps,'')){
#   exposure_instruments <- filter(exposure_clumped, SNP %in% selected_exposure_snps)
# } else{
#   print('No exposure instruments defined!')
# }

exposure_instruments <- exposure_clumped

# if(nrow(exposure_instruments ==0)){
#   print('You have no instruments for the exposure')
#   stop()
# }

#print(head(exposure_instruments))

#--------------------------------------------------------------------------------
#Load the outcome GWAS summary statistics
tic()
print('Importing the outcome GWAS summary statistics')
if(str_detect(outcome_input_path, '_manhattan_rsid\\.tsv|_manhattan_rsid2\\.tsv')){
  test <- read_tsv(outcome_input_path) %>%
    mutate(Phenotype=outcome_name)
  
  outcome_summstats <- format_data(test,snps=exposure_instruments$SNP, 
                                         type='outcome',
                                          snp_col='rsid',beta_col='beta',se_col='beta_se',eaf_col='eaf',effect_allele_col='allele_B',other_allele_col='allele_A',
                                           pval_col='pval',chr_col='chromosome',pos_col='position', samplesize_col='all_total')

} else if(str_detect(outcome_input_path, '_manhattan_rsid_logpval\\.tsv|_manhattan_rsid2_logpval\\.tsv')){
  test <- read_tsv(outcome_input_path) %>%
    mutate(Phenotype=outcome_name)
  
  outcome_summstats <- format_data(test,snps=exposure_instruments$SNP, 
                                   type='outcome',
                                   snp_col='rsid',beta_col='beta',se_col='beta_se',eaf_col='eaf',effect_allele_col='allele_B',other_allele_col='allele_A',
                                   pval_col='pval',chr_col='chromosome',pos_col='position', samplesize_col='all_total', log_pval = T)

  
} else if (str_detect(outcome_input_path, 'gwama')){
  test <- read_tsv(outcome_input_path) %>%
    mutate(Phenotype=outcome_name)
  
  outcome_summstats <- format_data(test,snps=exposure_instruments$SNP, #Only extract the exposure instruments
                              type='outcome',
                                          snp_col='rsid',beta_col='beta',se_col='se',eaf_col='eaf',effect_allele_col='reference_allele',other_allele_col='other_allele',
                                         pval_col='p-value', samplesize_col='n_samples')

} else if (str_detect(outcome_input_path, '_metal.TBL.gz')){
  test <- read_tsv(outcome_input_path) %>%
    mutate(Phenotype=outcome_name) %>%
    mutate(samplesize=19260) #For the Tadros 2021 CMR GWAS
  
  outcome_summstats <- format_data(test,snps=exposure_instruments$SNP, #Only extract the exposure instruments
                                   type='outcome',
                                   snp_col='MarkerName',beta_col='Effect',se_col='StdErr',eaf_col='Freq1',effect_allele_col='Allele2',other_allele_col='Allele1',
                                   pval_col='P-value')

}

toc()

if(nrow(outcome_summstats) != nrow(exposure_instruments)){
  print('Failed to identify the some of the exposure instruments in the outcome GWAS so ignoring those instruments')
  missing_snps <- exposure_instruments$SNP[!exposure_instruments$SNP %in% outcome_summstats$SNP]
  print(str_c('The following SNPs are not present in the outcome summary statistics: ',missing_snps))
}
rm(test)
#print(head(outcome_summstats))

#-------------------------------------------------------------------------------
#Harmonise the data
if(isFALSE(multivariate)){
  print('Harmonising exposure and outcome data')
  harmonised_data <- harmonise_data(exposure_instruments, outcome_summstats,
                                    action=1) #Assume that all alleles are presented on the forward strand
  rm(exposure_clumped, exposure_instruments, exposure_summstats, outcome_summstats)
  
  #Compute variance for binary traits given error in report Calculating approximate SNP-exposure and/or SNP-outcome correlations, assuming all are quantitative traits. Please pre-calculate r.exposure and/or r.outcome using get_r_from_lor() for any binary traits
  
  # if (exposure_name %in% c('hcm')){
  #   if(max(harmonised_data$samplesize.exposure) == 34914){ #This is without HCMR and with UKB but with everything else so
  #     #NL = 999 + 2117 (cases + controls) As per Tadros Supp.
  #     #BRRD = 239 + 7203
  #     #CAN = 1035 + 13889
  #     #ITALY = 277 + 1293
  #     #RBH = 448 + 1219 
  #     #GEL = 470 + 2355
  #     #UKB = 674 + 2695
  #     ncases <- 4142
  #     ncontrols <-30772
  #     prevalence <- 1/543 #As per McGurk et al, 2023
  #   } else if(max(harmonised_data$samplesize.exposure) == 74259){ #This is with HCMR and without UKB
  #     ncases <- 2431+999+239+1035+277+448+470
  #     ncontrols <-40283+2117+7203+13889+1293+1219+2355
  #     prevalence <- 1/543 #As per McGurk et al, 2023
  #   }
    
  #   harmonised_data <- harmonised_data %>%
  #     mutate(ncase.exposure = ncases,ncontrol.exposure=ncontrols) %>%
  #     mutate(r.exposure=get_r_from_lor(lor=beta.exposure,af=eaf.exposure, ncase.exposure, ncontrol.exposure,
  #                                      prevalence = rep(prevalence, nrow(harmonised_data))))
    
  # }else if (outcome_name %in% c('hcm')) {
  #   if(max(harmonised_data$samplesize.outcome) == 34914){ #This is without HCMR and with UKB but with everything else so
  #     ncases <- 4142
  #     ncontrols <-30772
  #     prevalence <- 1/543 #As per McGurk et al, 2023
  #   } else if(max(harmonised_data$samplesize.outcome) == 74259|max(harmonised_data$samplesize.outcome) == 72592){ #This is with HCMR and without UKB
  #     ncases <- 2431+999+239+1035+277+448+470
  #     ncontrols <-40283+2117+7203+13889+1293+1219+2355
  #     prevalence <- 1/543 #As per McGurk et al, 2023
  # }
  
  # harmonised_data <- harmonised_data %>%
  #   mutate(ncase.outcome = ncases,ncontrol.outcome=ncontrols) %>%
  #   mutate(r.outcome=get_r_from_lor(lor=beta.outcome,af=eaf.outcome, ncase.outcome, ncontrol.outcome,
  #                                   prevalence = rep(prevalence, nrow(harmonised_data))))
  # }
}

#---------------------------------------------------------------------------------
#If multivariate MR, harmonise the exposure and outcome data here i.e double check that the effect alleles refer to the same
if(isTRUE(multivariate)){
  #Exposure side
  exposure_instruments <- filter(exposure_instruments, SNP %in% outcome_summstats$SNP)
  
  #Outcome side - if any of the effect alleles don't match, reverse the sign of the beta of the outcome to convert it SAME EA as the EXPOSURES!
  outcome_summstats_harmonised <- left_join(outcome_summstats, select(exposure_instruments, SNP, EA)) %>%
    mutate(EA_nonharmonised = ifelse(effect_allele.outcome == EA, F, T)) %>%
    rowwise() %>%
    mutate(beta.outcome = ifelse(isTRUE(EA_nonharmonised), -beta.outcome, beta.outcome),
           effect_allele.outcome =ifelse(isTRUE(EA_nonharmonised), other_allele.outcome, effect_allele.outcome),
           other_allele.outcome =ifelse(isTRUE(EA_nonharmonised), effect_allele.outcome, other_allele.outcome)) %>%
    select(-EA_nonharmonised,-EA)
  
  #Harmonising together into single tb
  harmonised_data <- left_join(exposure_instruments, outcome_summstats_harmonised, by='SNP')
  if (identical(harmonised_data$EA, harmonised_data$effect_allele.outcome)){print('SNPs harmonised')}
  
  #Extract the relevant matrices for MendelianRanomization Multivariate MR
  exposure_beta_cols <- harmonised_data[str_detect(colnames(harmonised_data),'_beta')]
  exposure_se_cols <- harmonised_data[str_detect(colnames(harmonised_data),'_se')]
  exposure_names <- str_match(colnames(harmonised_data), '(.+)_beta')[,2][which(!is.na(str_match(colnames(harmonised_data), '(.+)_beta')[,2]))]
  
  outcome_beta_col <- harmonised_data$beta.outcome
  outcome_se_col <- harmonised_data$se.outcome
  
  mvmr_input <- MendelianRandomization::mr_mvinput(bx=as.matrix(exposure_beta_cols), bxse=as.matrix(exposure_se_cols),
                                                   by=outcome_beta_col, byse=outcome_se_col, exposure=exposure_names, outcome=outcome_name, 
                                                   snps=harmonised_data$SNP, effect_allele=harmonised_data$EA, other_allele = harmonised_data$NEA)
    
  rm(exposure_clumped, exposure_instruments, outcome_summstats, outcome_summstats_harmonised)
}

#--------------------------------------------------------------------------------
#Run MR
print('Running MR')

if(isFALSE(multivariate)){
  #Mainline MR analyses
  ##Select IVW with MRE if n_instruments >5, otherwise fixed-effects as per Burgess et al, 2023
  if(nrow(harmonised_data)==1){
    mainline_mr_method_list <- c('mr_wald_ratio')
  } else {
    mainline_mr_method_list <- c('mr_ivw')
  }
  mainline_mr_results <- mr(harmonised_data, method_list=mainline_mr_method_list) #Runs MR with a bunch of different methods- call mr_method_list() to see\

  #Sensitivity analyses
  sensitivity_mr_method_list <- c('mr_ivw','mr_weighted_median','mr_weighted_mode','mr_egger_regression')
  sensitivity_mr_results <- mr(harmonised_data, method_list=sensitivity_mr_method_list)
  pleiotropy_results <- mr_pleiotropy_test(harmonised_data) #MR Egger test for directional pleiotropy via evaluating MR Egger intercept difference from 0
  res_single <- mr_singlesnp(harmonised_data) #Single-SNP analysis
  res_loo <- mr_leaveoneout(harmonised_data) #Leave-one-out analysis

  #Plotting plots
  print('Plotting output plots')
  theme_set(theme_classic())
  scatter <- mr_scatter_plot(mainline_mr_results, harmonised_data)
  scatter2 <- mr_scatter_plot(sensitivity_mr_results, harmonised_data)
  # print(scatter[[1]])
  forest <- mr_forest_plot(res_single)
  #print(forest[[1]])
  loo_forest <- mr_leaveoneout_plot(res_loo)
  #print(loo_forest[[1]])

  #Add an output step to output odds-ratio with 95% confidence interval (instead of log-odds with standard error)
  safe_generate_odds_ratios <- safely(generate_odds_ratios)
  or_result <- safe_generate_odds_ratios(sensitivity_mr_results)

  # Output the harmonised instrument details (Beta and SE for exposure/outcome)
  print("Writing individual instrument details (beta/se) to TSV file.")
  instrument_details_df <- harmonised_data %>%
    dplyr::select(
      SNP,
      effect_allele.exposure, # Keep alleles for reference
      other_allele.exposure,
      eaf.exposure,          # Keep EAF for context
      beta.exposure,
      se.exposure,
      pval.exposure,         # Keep p-values
      beta.outcome,
      se.outcome,
      pval.outcome,
      mr_keep              # Indicates if SNP was kept after harmonisation checks
    )
  
  instrument_details_filename <- file.path(output_path, paste0(exposure_name, '_to_', outcome_name, '_instrument_details.tsv'))
  write_tsv(instrument_details_df, instrument_details_filename)
  print(paste("Instrument details saved to:", instrument_details_filename))
  
  #MendelianRandomization Addendum - computing an approximation of the first-stage F-statistic from genetic variants -> exposure
  # mr_object <- MendelianRandomization::mr_input(bx=harmonised_data$beta.exposure, bxse=harmonised_data$se.exposure,
  #                                               by=harmonised_data$beta.outcome, byse= harmonised_data$se.outcome,
  #                                               exposure=exposure_name,
  #                                               outcome=outcome_name,
  #                                               snps=harmonised_data$SNP)
  # MRpackage_ivw_results <- MendelianRandomization::mr_ivw(mr_object, model='default')

  #Output plots and data --------------------------------------------------------------------------------
  #Write out a table of results
  if (isFALSE(dir.exists(str_c(output_path,exposure_name,'_to_',outcome_name,'/')))){
    dir.create(str_c(output_path,exposure_name,'_to_',outcome_name,'/')) #Create the directory if it doesn't exist already
  }

  #This outputs the print output as a txt file
  # sink(file=str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_MRpackage_IVW_results.txt'))
  # print(MRpackage_ivw_results)
  # sink(file=NULL)

  #Save 
  ggsave(str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_scatter_mainline.png'),scatter[[1]],dpi=600)
  ggsave(str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_scatter_sensitivity.png'),scatter2[[1]],dpi=600)
  ggsave(str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_forest.png'),forest[[1]],dpi=600)
  ggsave(str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_looforest.png'),loo_forest[[1]],dpi=600)

  write.table(mainline_mr_results,str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_results_mainline.tsv'),sep='\t')
  write.table(sensitivity_mr_results,str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_results_sensitivity.tsv'),sep='\t')
  if (!is.null(or_result$result)){ 
    write.table(or_result$result,str_c(output_path,exposure_name,'_to_',outcome_name,'/',exposure_name,'_mr_results_sensitivity_oddsratio.tsv'),sep='\t')
  } else {
    print('Odds ratio function failed')
  }

  #Print out a HTML report
  print('Constructing HTML Report')
  safe_mrreport <- safely(mr_report)
  safe_mrreport(harmonised_data, output_path=str_c(output_path,exposure_name,'_to_',outcome_name)) #Generate html report
#--------------------------------------------------------------------------------
#Run multivariate MR
} else if(isTRUE(multivariate)){
  
  if(str_detect(exposure_input_path,'tadros21')){ exposure_gwas_samplesize=19260 } #This defines a vector for the sample size of the GWAS used for each of the exposures
  
  #ivw_results_intercept0 <- MendelianRandomization::mr_mvivw(mvmr_input, model='default',correl=F,nx=exposure_gwas_samplesize, distribution='normal') #Default is fixed effect model with <=3 genetic variants and random effects if more variants
  
  if(str_detect(exposure_input_path,'contractility')){ n_pcs = 3 } else if(str_detect(exposure_input_path,'allcmr_noextra')){ n_pcs = 5 } else if(str_detect(exposure_input_path,'allcmr?')){ n_pcs = 8 }
  
  #This is multivariate MR with sparse PCA beforehand as per Karageorgiou et al, 2023
  spca_mvmr <- sca_MR(mvmr_input$betaX, mvmr_input$betaY, mvmr_input$betaYse,nfold=5,N_PC=n_pcs, spars_length=n_pcs)
  
  #Write out a table of results
  if (isFALSE(dir.exists(str_c(output_path,'plots/',exposure_name,'_to_',outcome_name,'/')))){
    dir.create(str_c(output_path,'plots/',exposure_name,'_to_',outcome_name,'/')) #Create the directory if it doesn't exist already
  }
  
  #This outputs the print output as a txt file
  sink(file=str_c(output_path,'plots/',exposure_name,'_to_',outcome_name,'/',exposure_name,'_sPCA_MVMR_IVW_results.txt'))
  print(spca_mvmr)
  sink(file=NULL)
  
  saveRDS(spca_mvmr,file=str_c(output_path,'plots/',exposure_name,'_to_',outcome_name,'/',exposure_name,'_sPCA_MVMR_IVW_object.rds'))
  
}

