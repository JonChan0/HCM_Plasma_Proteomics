#Script to summarise for the top-ranked causal SNPs their GWAS summary statistics from REGENIE; FINEMAP causal configuration statistics; OpenTargetGenetics V2G statistics.
#Author: Jonathan Chan
#Date: 2024-07-04

args <- commandArgs(trailingOnly=T)

if(length(args) < 1){print('No arguments defined')}

library(tidyverse)
library(ggrepel)
library(liftOver)
theme_set(theme_classic())


pheno <- args[1]
basepath <- args[2]
summstat_folder <- args[3]
chainfilepath <- args[4]

#Local test code
# pheno <- 'cicp'
# basepath <- 'hcmr/'
# summstat_folder <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'
# chainfilepath <- '/well/PROCARDIS/jchan/bin/liftover/hg38ToHg19.over.chain'

#Import/Format-----------------------------------------------------------------------------

  ## FINEMAP Import/Format-------
finemap_importer_formatter <- function(pheno, basepath){
  
  #Import in the FINEMAP overallconfig.tsv which contains the posterior probabilities of x being a causal SNP; and Bayes Factors
  if(str_detect(basepath,'ukb')){
    finemap_configtsv_path <- str_c('../output/',basepath, '/', pheno, '/',pheno,'_dataset_overallconfig.tsv')
  } else if (str_detect(basepath,'hcmr')){
    finemap_configtsv_path <- str_c('../output/',basepath,'/', pheno, '/',pheno,'_dataset_overallconfig_remapped.tsv')
  }
  
  finemap_tsv <- read_tsv(finemap_configtsv_path, col_types=c('nncc',rep('n',9)))
  
  toprank <- finemap_tsv %>%
    filter(rank==1) %>% #Only grab top configuration
    dplyr::select(genomic_region,rsid, config,prob) %>%
    mutate(rsid = ifelse(str_detect(rsid,','), str_split(rsid,', '), rsid)) #This is for multiple SNPs in a single configuration
  
  toprank <- toprank %>%
    mutate(rsid = map(rsid, enframe))%>%
    unnest(rsid) %>%
    dplyr::rename(rsid=value) %>%
    mutate(snptest_rsid=rsid) %>% #This is the maintain the SNPtest RSID for use in filtering against manhattan_rsid.tsv
    mutate(rsid=ifelse(!str_detect(rsid,'rs'),str_replace(rsid,':','_'),rsid)) %>% #This is in case the variant does not have RSID (and obeys OTG notation)
    dplyr::rename(SNP=rsid)
  
  toprank <- toprank %>%
    mutate(config = ifelse(str_detect(config,','), str_split(config,','), config))
  
  for (i in seq(nrow(toprank))){
    if(class(toprank$config[[i]])=='character'){
      index <- toprank$name[[i]]
      toprank$config[[i]]<- toprank$config[[i]][[index]]
    } else{
      toprank$config[[i]] <- toprank$config[[i]]
    }
  }
  
  toprank <- mutate(toprank, config=as.character(config))
  
  return(toprank)
}
  
  ## Summstats Import--------
  
summstat_importer_formatter <- function(input_path, pheno, rsids_of_interest, pattern='_manhattan_rsid.tsv'){
  
    tsvpath <- str_c(input_path,'/',list.files(input_path, pattern=str_c(pheno,pattern)))
    
    summstats<- read_tsv(tsvpath) %>%
      filter(rsid %in% rsids_of_interest) %>% #Filters for only SNPs which are in the finemap top-ranked configurations
      dplyr::rename(SNP=rsid) %>%
      mutate(pheno=pheno)%>%
      dplyr::select(-snid)
    
    if(str_detect(input_path, 'hcmr')){
      summstats<- summstats%>% dplyr::select(-snid_allele_A_allele_B)
    }
    
    write_tsv(summstats, str_c('../output/',basepath,pheno,'_FINEMAP_toprank_GWASsummstats.tsv')) #Outputs a file with the GWAS summary statistics for the SNPs which are in FINEMAP top-ranked causal configs
    return(summstats)
  }
  
  ## V2G OTG Import-----------

v2g_importer_formatter <- function(pheno, basepath){
  
  #Import in the V2G summary tables for each phenotype which includes the L2G scores for the top-ranked configuration at each SNP
  otargen_l2g_summary_path <- str_c('../output/',basepath, pheno,'/v2g/',pheno,'_v2g_summary_l2g.tsv')
  
  v2g_tsv <- read_tsv(otargen_l2g_summary_path,show_col_types = FALSE) %>%
    filter(!is.na(overallScore))
  return(v2g_tsv)
  
  }

#Merger---------------------------------------------------------------------------------

summstats_finemap_v2g_merger <- function(summstats2, finemap_tsv, v2g_tsv){
  
  #Harmonise the summstats2 chr notation to use _ instead of : for the first one detected
  summstats2 <- mutate(summstats2,SNP=ifelse(!str_detect(SNP,'rs'),str_replace(SNP,':','_'),SNP))
  
  merged <- full_join(summstats2, finemap_tsv, by='SNP') #Merges the summary statistics and the FINEMAP statistics
  
  #Merges the V2G statistics with the  others
  merged <- full_join(merged, v2g_tsv, by=c('genomic_region','SNP'))
  
  return(merged)
}


#Liftover--------------------------------------------------------------------------
#Need to map the nonRSID (b38 coordinates from the L2G back to the b37 coordinates for FINEMAP/GWAS summstats)
#Likely using genomic region and alleles

b38_b37_remapper <- function(input_tb, pheno_name, chainfilepath){
  #filter for only rows which have NA at snptest_rsid so are from the V2G whereby their nonrsIDs were lifted over to b38 to pass to OTG
  b38 <- input_tb %>%
    filter(is.na(snptest_rsid))
  
   #Need to liftOver the chromosomal positions of the b38 SNPS which aren't rsIDs
  b38 <- b38 %>%
    mutate(chr = str_match(variant, '^(\\d{1,2})')[,2]) %>%
    mutate(start=str_match(variant,'^[^_]+_(\\d+)_')[,2]) %>%
    mutate(string=str_c('chr',chr,':',start,'-',start))
  
  ## Convert to GRanges object and LiftOver
  granges <- as(b38$string, "GRanges")
  b37 <- liftOver(granges, import.chain(chainfilepath))
  
    #Remap back to the original tibble
  output_tb <- dplyr::select(b38, c(genomic_region, gene.symbol:gene_rank)) %>%
    bind_cols(cbind(unlist(as.character(seqnames(b37))),unlist(start(b37))))
  colnames(output_tb) <- c('genomic_region','gene.symbol','b38','overallScore','gene.id','rank', 'chr','pos')
  
   output_tb<-  output_tb %>%
    mutate(b37 = str_c(str_extract(chr,'\\d+'),pos,str_match(b38,'_([A-Za-z]+)_[A-Za-z]+$')[,2],str_match(b38,'_[A-Za-z]+_([A-Za-z]+)$')[,2],sep='_')) %>%
    dplyr::select(-chr, -pos)
  
  out <- left_join(filter(input_tb, !is.na(snptest_rsid)),output_tb, by=c('SNP'='b37','genomic_region')) %>%
    mutate(gene.symbol=coalesce(gene.symbol.x, gene.symbol.y),
           overallScore=coalesce(overallScore.x, overallScore.y),
           gene.id = coalesce(gene.id.x, gene.id.y))%>%
    dplyr::select(-contains('.x'), -variant, -rank, -contains('.y')) %>%
    arrange(genomic_region)
  
  write_tsv(out, str_c('../output/',basepath,pheno_name,'_FINEMAP_toprank_allstats.tsv')) #Write out all the statistics for that phenotype after lifting back
  
  return(out)
}

#Plotter-----------------------------------------------------------------------

scatter_plotter <- function(input_tb, output_basepath, log10pval_or_pval='pval'){
  
  input_tb <- filter(input_tb, !is.na(pval)) %>%#Removes all rows without a summary statistic pvalue
    mutate(pval=as.numeric(pval))
  
  pheno <- unique(input_tb$pheno)
  output_path <- str_c(output_basepath, pheno,'/', collapse='')
  
  input_tb <- mutate(input_tb, genomic_region = factor(genomic_region))
  
  if(log10pval_or_pval=='pval'){
    if(length(unique(input_tb$genomic_region))<=8 && length(unique(input_tb$SNP))<= 12){
      
      scatter <- ggplot(input_tb, aes(x=overallScore, y=-log10(pval)))+
        geom_point(aes(size=prob, col=SNP, shape=genomic_region), alpha=0.5)+
        scale_colour_brewer(name='SNP',palette='Dark2')+
        scale_shape_discrete(name='Genomic Region')+
        guides(col=guide_legend(ncol=2), shape=guide_legend(ncol=2))
      
    } else {
      scatter <- ggplot(input_tb, aes(x=overallScore, y=-log10(pval)))+
        geom_point(aes(size=prob), alpha=0.5)
    }
    scatter <- scatter+
      ylab('-log10(GWAS p-value)') 
  } else if (log10pval_or_pval=='log10pval'){
    
    if(length(unique(input_tb$genomic_region))<=8 && length(unique(input_tb$SNP))<= 12){
      
      scatter <- ggplot(input_tb, aes(x=overallScore, y=pval))+
        geom_point(aes(size=prob, col=SNP, shape=genomic_region), alpha=0.5)+
        scale_colour_brewer(name='SNP',palette='Dark2')+
        scale_shape_discrete(name='Genomic Region')+
        guides(col=guide_legend(ncol=2), shape=guide_legend(ncol=2))
      
    } else {
      scatter <- ggplot(input_tb, aes(x=overallScore, y=pval))+
        geom_point(aes(size=prob), alpha=0.5)
    }}
  
  
  scatter <- scatter +
    geom_text_repel(aes(label=gene.symbol),force=3, size=3)+
    scale_size_continuous(name=str_wrap('FINEMAP Posterior Probability',width=20))+
    xlab('OpenTarget Genetics L2G Score')+
    ylab('-log10(GWAS p-value)') +
    labs(title=str_wrap(str_c('GWAS; FINEMAP; Open Target Genetics L2G Summary for Top-Ranked Finemapped SNPs for ', pheno),width=70),
         caption=str_wrap('Each point refers to one potential prioritised gene per SNP and each shape represents one genomic region',width=70))
  
  
  print(scatter)
  
  ggsave(str_c(output_path,pheno,'_gwas_finemap_v2g_summary.png'),scatter, dpi=600)
}

#Summariser-------------------------------------------------------------------------
#I also output summary TSV files for each genomic region with its top-ranked linked gene.

perpheno_summariser <- function(summ, pheno_name, basepath){
  
  summ2 <- summ %>%
    #filter(pval >= 10^-(1.7e-11)) %>% #For -log10pval form of FDR used by Sun et al, 2023
    group_by(genomic_region, name, pheno) %>%
    filter(overallScore == max(overallScore)) %>%
    ungroup() %>%
    group_by(genomic_region, pheno) %>%
    # filter(pval == min(pval)) %>%
    ungroup() %>%
    dplyr::select(genomic_region,chromosome, position, SNP, gene.symbol, allele_A, allele_B, maf, pheno, beta, beta_se, pval, overallScore) %>%
    arrange(desc(pval))
  
  write_tsv(summ2, str_c('../output/',basepath, pheno_name,'_summstats_topV2G.tsv'))
  
}

#Main function------------------------------------------------------------------
main <- function(pheno, basepath, summstat_folder, chainfilepath){
  
  finemap_tsv <- finemap_importer_formatter(pheno, basepath)
  
  if(str_detect(basepath,'ukb')){
    summstat_pattern <- '_manhattan_rsid_logpval.tsv'
  } else if (str_detect(basepath,'hcmr')){
    summstat_pattern <- '_manhattan_rsid2_logpval.tsv'
  }
  
  summstats <-summstat_importer_formatter(summstat_folder, pheno, finemap_tsv$snptest_rsid,pattern=summstat_pattern)
  
  v2g_tsv <- v2g_importer_formatter(pheno, basepath)
  
  #Double check that the number of SNPs is equal in both the summstats and the finemap_tsv
  summstats2 <- unique(summstats)
  print(str_c('There is a difference in number of SNPs between the summary statistics and the FINEMAP config.tsv of ',nrow(summstats2)-nrow(finemap_tsv)))
  
  #Merge
  merged <- summstats_finemap_v2g_merger(summstats2, finemap_tsv, v2g_tsv)
  merged <- b38_b37_remapper(merged,pheno,chainfilepath)
  
  #Plot
  scatter_plotter(merged, str_c('../output/',basepath),log10pval_or_pval='log10pval')   #Lack of dot means V2G data unavailable for that variant
  
  #Summarise
  perpheno_summariser(merged, pheno, basepath)
}

#Carry out the function
main(pheno, basepath, summstat_folder, chainfilepath)

#Run locally on RStudioServer for all pp
#For HCMR
# phenos <- c('NTproBNP', 'TnTStat','cicp','gal3','mmp1','st2') #For HCMR
# basepath <- 'hcmr/'
# summstat_folder <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'
# chainfilepath <- '/well/PROCARDIS/jchan/bin/liftover/hg38ToHg19.over.chain'

phenos <- c('hcmr2_ecvfwhole','hcmr2_lge_total') #For HCMR
basepath <- 'hcmr/'
summstat_folder <- '../../2_gwas/output/gwas/hcmr/REGENIE/step2/formatted/'
chainfilepath <- '/well/PROCARDIS/jchan/bin/liftover/hg38ToHg19.over.chain'

walk(phenos, ~main(., basepath, summstat_folder, chainfilepath))

#For UKB
# phenos <- c('NTproBNP','MMP1', 'IL1RL1', 'LGALS3',
#         'NPPB','ACE2','EDN1','HRC',
#         'APOM','F7'
#         #,TNNI3
#         )

# phenos <- c('SHISA5','MAMDC2','FABP3','LTBP2')

# phenos <- c('ANGPT2','STC2')
# basepath <- 'ukbpp_nonHarper_includeHCM/'
# # summstat_folder <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/formatted/'
# # summstat_folder <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part2/formatted/'
# summstat_folder <- '../../2_gwas/output/gwas/ukb/REGENIE/step2/2_inclHCMcases/part3/formatted/'
# chainfilepath <- '/well/PROCARDIS/jchan/bin/liftover/hg38ToHg19.over.chain'
# 
# walk(phenos, ~main(., basepath, summstat_folder, chainfilepath))

