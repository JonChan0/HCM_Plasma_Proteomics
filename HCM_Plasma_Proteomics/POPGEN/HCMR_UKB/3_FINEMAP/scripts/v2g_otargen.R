#Script to run V2G via Open Target Genetics downstream of the FINEMAP step on the top-ranked SNP causal configurations at each genomic region.
#Author: Jonathan Chan
#Date: 2024-05-29

library(tidyverse)
library(otargen)
library(fs)
library(rmarkdown)
library(liftOver)

args = commandArgs(trailingOnly=TRUE) #Allows input of arguments from Rscript
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument (i.e the input vep.txt files) must be supplied", call.=FALSE)
}

input_file <- args[1]
pheno <- args[2]
output_basepath <- args[3]
l2g_score_cutoff <- args[4]
Rmd_template_path <- args[5]
liftover_chain_filepath <- args[6]

#Local test code for HCMR
# pheno <- 'NTproBNP'
# input_file <- str_c('popgen/4_finemap/output/hcmr/',pheno,'/',pheno,'_dataset_overallconfig_remapped.tsv')
# output_basepath <- str_c('popgen/4_finemap/output/hcmr/',pheno,'/v2g/')
# l2g_score_cutoff <-0.1
# Rmd_template_path <- 'popgen/4_finemap/scripts/v2g_otargen_tables_template.Rmd'
# liftover_chain_filepath <- 'popgen/4_finemap/scripts/hg19ToHg38.over.chain'


#Make the entire script a function
v2g_otargen_overall <- function(input_file, pheno, output_basepath, l2g_score_cutoff, Rmd_template_path, liftover_chain_filepath){
  
  print(str_c('Running OpenTarget Genetics L2G for ', pheno))
    
  if (!dir_exists(output_basepath)) {
    dir_create(output_basepath, recurse = TRUE)
  }
  
  #Import-----------------------------------------------------------------------------
  finemapped_loci <- read_tsv(input_file, col_types=c('nnccnnnnnncnn'))
  
  #Filter----------------------------------------------------------------------------
  #Filter for the top configuration at each genomic locus to create a tibble of all unique variants which are deemed highest-rank for causality
  
  toprank <- finemapped_loci %>%
    filter(rank==1) %>% #Only grab top configuration
    dplyr::select(genomic_region,rsid) %>%
    mutate(rsid = ifelse(str_detect(rsid,','), str_split(rsid,', '), rsid))
  
  toprank <- toprank %>%
    mutate(rsid = map(rsid, enframe))%>%
    unnest(rsid) %>%
    dplyr::rename(rsid=value) %>%
    mutate(rsid=ifelse(!str_detect(rsid,'rs'),str_replace(rsid,':','_'),rsid)) #This is in case the variant does not have RSID
  
  #Need to liftOver the chromosomal positions of the b37 SNPS which aren't rsIDs
  nonrsID <- filter(toprank, !str_detect(rsid,'rs')) %>%
    mutate(chr = str_match(rsid, '^(\\d{1,2})')[,2]) %>%
    mutate(start=str_match(rsid,'^[^_]+_(\\d+)_')[,2]) %>% #The _ is for the presence of the allele codes afterwards
    mutate(string=str_c('chr',chr,':',start,'-',start))
  
  ## Convert to GRanges object and LiftOver
  nonrsID_granges <- as(nonrsID$string, "GRanges")
  b38_nonrsID <- liftOver(nonrsID_granges, import.chain(liftover_chain_filepath))
  
  #Remap back to the toprank tibble
  v3_nonrsID_tb <- dplyr::select(nonrsID, -c('chr','start','string')) %>%
    bind_cols(cbind(unlist(as.character(seqnames(b38_nonrsID))),unlist(start(b38_nonrsID))))
  colnames(v3_nonrsID_tb) <- c('genomic_region','name','rsid','chr','pos')
  v3_nonrsID_tb <- v3_nonrsID_tb %>%
    mutate(b38 = str_c(str_extract(chr,'\\d+'),pos,str_match(rsid,'_([A-Za-z]+)_[A-Za-z]+$')[,2],str_match(rsid,'_[A-Za-z]+_([A-Za-z]+)$')[,2],sep='_')) %>%
    dplyr::select(-chr, -pos, -rsid) %>%
    dplyr::rename(rsid=b38)
  
  toprank <- bind_rows(filter(toprank, !rsid %in% nonrsID$rsid),v3_nonrsID_tb) %>%
    arrange(genomic_region, name)
  
  rm(b38_nonrsID, nonrsID, nonrsID_granges,v3_nonrsID_tb)
  
  #V2G Call via Open Target Genetics-----------------------------------------------
  safe_v2g <- safely(genesForVariant) #Use purrr:safely to wrap my function and get it to output either the result or the error message
                              #This is especially useful for occasions where OTG doesn't have the variant of interest
  v2g_raw <- map(toprank$rsid,~safe_v2g(.))
  names(v2g_raw) <- toprank$rsid
  
  #For those that fail, try again with the snID instead-----------------------
  failed_v2g <- names(v2g_raw)[map_lgl(v2g_raw, ~is.null(.$result))]
  
  if(length(failed_v2g) >= 1 && str_detect(input_file,'hcmr')){
    # Get snIDs for failed rsIDs (this is in b37)
    failed_snids <- finemapped_loci %>%
      filter(str_detect(rsid, paste(failed_v2g, collapse = "|"))) %>%
      dplyr::select(genomic_region, rank,rsid, snid) %>%
      separate_rows(rsid, sep = ", ") %>%
      filter(rsid %in% failed_v2g)
    
    #Liftover to b38
    snids_to_lift <- failed_snids %>%
      mutate(chr = str_match(snid, '^(\\d{1,2})')[,2],
             start = str_match(snid, ':(\\d+)_')[,2],
             string = str_c('chr', chr, ':', start, '-', start))
    
    # Convert to GRanges object and LiftOver
    snids_granges <- as(snids_to_lift$string, "GRanges")
    b38_snids <- liftOver(snids_granges, import.chain(liftover_chain_filepath))
    
    #Remap back to the failed_snids tibble
    snids_to_lift <- dplyr::select(snids_to_lift, -c('chr','start','string')) %>%
      bind_cols(cbind(unlist(as.character(seqnames(b38_snids))),unlist(start(b38_snids))))
    colnames(snids_to_lift) <- c('genomic_region','rank','rsid','snid','chr','pos')
    
    failed_snids_b38 <- snids_to_lift %>%
      mutate(b38 = str_c(str_extract(chr,'\\d+'),pos,str_match(snid,'_([A-Za-z]+)_[A-Za-z]+$')[,2],str_match(snid,'_[A-Za-z]+_([A-Za-z]+)$')[,2],sep='_')) %>%
      dplyr::select(-chr, -pos, -snid)
    
    rm(failed_snids, failed_v2g, snids_to_lift, snids_granges, b38_snids)
    
    v2g_raw2 <- map(failed_snids_b38$b38,~safe_v2g(.))
    names(v2g_raw2) <- failed_snids_b38$rsid #Label with the old one to enable remapping to manhattan_rsid2.tsv
    
    v2g_raw <- v2g_raw[map_lgl(v2g_raw, ~!is.null(.$result))] #Filter for those which contains non-null results
    v2g_raw <- c(v2g_raw, v2g_raw2) #Append the new data for the failed results
    rm(v2g_raw2)
    
    #N.B This code keeps the old rsID which failed as the 'rsid' and uses another column called b38 which corresponds to the b38 call which was actually made to OTG
    #This enables matching to precursor files in the pipeline
    toprank <- toprank %>%
      filter(!rsid %in% failed_snids_b38$rsid) %>%
      bind_rows(failed_snids_b38 %>% dplyr::rename(name=rank)) %>%
      arrange(genomic_region, name)
    
  }

  #Reformat V2G for output----------------------------------------------------------
  #This outputs 4x dataframes: 1) Overall summary of V2G using the LG2 score for each lead variant = one tibble
  #For each lead variant, a knitted HTML containing the for the prioritised genes ()
  #2) Colocalisation scores with QTLs/chromatin data for each lead variant
  #3) QTL associations for each lead variant;
  #4) Chromatin interactions for each lead variant
  
  #N.B The L2G model produced a well-calibrated score, ranging from 0 to 1, which reflects the approximate fraction of gold standard positive genes among all genes above a given threshold (Fig. 2). 
  #Extract L2G score for each causal variant
  topv2g_extractor <- function(v2g_raw_listelement, SNP_name,l2g_score_cutoff){
    
    if(class(v2g_raw_listelement$result) == 'NULL' || length(v2g_raw_listelement$result)==0){
      print(str_c('No OpenTarget Genetics data for ', SNP_name))
      return(v2g_raw_listelement$result)
      
    } else{
      out <- v2g_raw_listelement$result$v2g %>%
        mutate(gene_rank = base::rank(plyr::desc(overallScore))) %>% #Rank refers to the per-SNP rank of associated genes
        filter(overallScore > l2g_score_cutoff | gene_rank==min(gene_rank,na.rm=T)) %>% #Take the genes with a score > score or top ranked gene for that variant (if it is a tie, that is fine)
        mutate(SNP = SNP_name) %>%
        dplyr::select(SNP, everything())
      
      return(out)
    }

  }
  
  overall_v2g_out <- map2(v2g_raw,names(v2g_raw), ~topv2g_extractor(.x,.y, l2g_score_cutoff)) %>%
    bind_rows() %>%
    full_join(toprank, by=c('SNP'='rsid')) %>%
    dplyr::select(genomic_region, within_region_SNPnumber=name, SNP,everything()) %>%
    arrange(genomic_region, desc(overallScore))
  
  #Graphing--------------------------------------------------------------------
  #Construct a graph of overall L2G scores for each lead variant
  theme_set(theme_classic())
  
  plot_tb <- overall_v2g_out %>%
    mutate(genomic_region=as.integer(genomic_region), SNP=factor(SNP)) %>%
    group_by(genomic_region,SNP) %>%
    mutate(gene.symbol=factor(gene.symbol)) %>%
    mutate(gene.symbol=reorder(gene.symbol, overallScore))
  
  if(length(unique(plot_tb$gene.symbol))<=12){
    
    l2g_summary_plot <- ggplot(plot_tb, 
                               aes(x=overallScore, y=SNP,group=gene.symbol))+
      geom_col(aes(fill=gene.symbol),
               position='dodge')+
      scale_fill_brewer(name='Gene',palette='Set3')
  } else{
    l2g_summary_plot <- ggplot(plot_tb, 
                               aes(x=overallScore, y=SNP,group=gene.symbol))+
      geom_col(position='dodge', col='white')
  }
  
  l2g_summary_plot <- l2g_summary_plot+
    geom_text(aes(label=gene.symbol), position=position_dodge(width=0.9), hjust=-0.2, size=3)+
    xlim(c(0,1))+
    #theme(axis.text.x=element_text(angle = 90))+
    xlab('Overall L2G Score') +
    ylab('SNP')+
    facet_wrap(~genomic_region, scale='free_y', dir='v',ncol=1)+
    labs(title=str_wrap(str_c('L2G Scores from Open Target Genetics for V2G from FINEMAP-derived causal SNPs for ', pheno),width=40))
  
  safe_ggsave <- safely(ggsave)
  safe_ggsave(str_c(output_basepath,pheno,'_l2g_summary_plot.png'),l2g_summary_plot,dpi=600, height = max(plot_tb$genomic_region) * 2.5)
  
  #Now knit a HTML looking for each of the lead variants, all the genes which pass the l2g_score_cutoff - the evidence for their association with that locus
  #Create a temporary .Rmd script
  
  #This grabs a list of dataframes for each lead variant where the data is only for the prioritised genes which pass the l2g_score_cutoff
  table_grabber <- function(v2g_raw_listelement, SNP_name,l2g_score_cutoff){
    
    if(class(v2g_raw_listelement$v2g) == 'NULL'){
      print(str_c('No OpenTarget Genetics tabular evidence data for ', SNP_name))
      return(list(v2g_raw_listelement$v2g))
      
    } else{
      SNPs <- v2g_raw_listelement$v2g %>%
        filter(overallScore > l2g_score_cutoff) #Take the genes with a score > score
      
      out <- list(
        tssd=filter(v2g_raw_listelement$tssd,gene.symbol %in% SNPs$gene.symbol),
        qtls=filter(v2g_raw_listelement$qtls,gene.symbol %in% SNPs$gene.symbol),
        chromatin=filter(v2g_raw_listelement$chromatin,gene.symbol %in% SNPs$gene.symbol)
      )
      return(out)
    }

  }
  
  output_support_tables <- map2(v2g_raw,names(v2g_raw), ~table_grabber(.x,.y,l2g_score_cutoff))
  names(output_support_tables) <- str_c(toprank$genomic_region, 'R_',names(output_support_tables))
  
  #Write out all the supporting tables for each causal variant
  walk2(output_support_tables, names(output_support_tables), 
        ~rmarkdown::render(Rmd_template_path, output_dir=output_basepath,output_file = str_c(.y, '_supporttables.html'), #This path is relative to the .Rmd file
                    params = list(leadvar = .y, tssd = .x$tssd, qtls = .x$qtls, chromatin=.x$chromatin),
                    clean=T))
  
  #Export as TSV----------------------------------------------------------
  
  write_tsv(overall_v2g_out, str_c(output_basepath,pheno,'_v2g_summary_l2g.tsv'))

}



#This is for running it as a normal Rscript with input arguments by CLI
v2g_otargen_overall(input_file, pheno, output_basepath, l2g_score_cutoff, Rmd_template_path, liftover_chain_filepath)


#Code to run locally for multiple traits


#For HCMR
# pp <- c('NTproBNP','mmp1', 'st2', 'gal3','TnTStat', 'cicp')
# walk(pp,
#      ~v2g_otargen_overall(
#        str_c('popgen/4_finemap/output/hcmr/',.,'/',.,'_dataset_overallconfig_remapped.tsv'),
#        .,
#        str_c('popgen/4_finemap/output/hcmr/',.,'/v2g/'),
#        0.1,
#        'popgen/4_finemap/scripts/v2g_otargen_tables_template.Rmd',
#        'popgen/4_finemap/scripts/hg19ToHg38.over.chain'
#      ))

#For HCMR
# pp <- c('hcmr2_ecvfwhole','hcmr2_lge_total')
# walk(pp,
#      ~v2g_otargen_overall(
#        str_c('popgen/4_finemap/output/hcmr/',.,'/',.,'_dataset_overallconfig_remapped.tsv'),
#        .,
#        str_c('popgen/4_finemap/output/hcmr/',.,'/v2g/'),
#        0.1,
#        'popgen/4_finemap/scripts/v2g_otargen_tables_template.Rmd',
#        'popgen/4_finemap/scripts/hg19ToHg38.over.chain'
#      ))


# FOr UKB - incl HCM
# pp <- c('NTproBNP','MMP1', 'IL1RL1', 'LGALS3',
#         'NPPB','ACE2','EDN1','HRC',
#         'APOM','F7'
#         #,TNNI3
#         )

# pp <- c('SHISA5','MAMDC2','LTBP2','FABP3')

# pp <- c('ANGPT2','STC2')
# 
# walk(pp,
#      ~v2g_otargen_overall(
#        str_c('popgen/4_finemap/output/ukbpp_nonHarper_includeHCM/','/',.,'_dataset_overallconfig.tsv'),
#        .,
#        str_c('popgen/4_finemap/output/ukbpp_nonHarper_includeHCM/',.,'/v2g/'),
#        0.1,
#        'popgen/4_finemap/scripts/v2g_otargen_tables_template.Rmd',
#        'popgen/4_finemap/scripts/hg19ToHg38.over.chain'
#      ))

# FOr UKB - excl HCM
# pp <- c('NTproBNP','MMP1', 'IL1RL1', 'LGALS3')
# walk(pp,
#      ~v2g_otargen_overall(
#        str_c('popgen/4_finemap/output/ukbpp_nonHarper_excludeHCM/',.,'/',.,'_dataset_overallconfig.tsv'),
#        .,
#        str_c('popgen/4_finemap/output/ukbpp_nonHarper_excludeHCM/',.,'/v2g/'),
#        0.1,
#        'popgen/4_finemap/scripts/v2g_otargen_tables_template.Rmd',
#        'popgen/4_finemap/scripts/hg19ToHg38.over.chain'
#      ))

