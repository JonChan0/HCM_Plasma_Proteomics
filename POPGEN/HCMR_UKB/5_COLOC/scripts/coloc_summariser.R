#Script to merge and plots results for multiple coloc results between two traits of interest
# Author: Jonathan Chan
# Date: 2024-10-11

## Inputs:
#Summary_results of interest
#Locusname
#Lead SNP

library(tidyverse)
theme_set(theme_classic())

args <- commandArgs(trailingOnly = T)

if(length(args)==0){
  print('No arguments defined')
}

results_folderpath <- args[1]

#Test local input
# results_folderpath<- '../output/plots/NTproBNP_HCM/'
# results_folderpath<- '../output/plots/TnTStat_HCM/'

#Import----------------------------

importer_merger <- function(path, output_path, traitA, traitB, results_prefix='abf'){

  loci <- dir(path)[!str_detect(dir(path), '\\.png|\\.tsv')]
  rsids <- vector('character', length(loci))
  summary_path <- vector('character', length(loci))
  
  for (i in seq_along(loci)){
    files <- dir(str_c(path, '/',loci[[i]],'/'))
    summary_path[[i]] <- str_c(path, '/', loci[[i]], '/',files[str_detect(files,'summaryresults.txt')])
    rsids[[i]] <- str_match(summary_path[[i]], str_c(results_prefix, '_(rs\\d+)_'))[,2]
  }
  
  results <- map(summary_path, ~read.table(., header=T)) %>% bind_rows() %>%
    mutate(locus = loci, rsid = rsids)
  
  write_tsv(results, str_c(output_path, '/',traitA, '_', traitB,'_summary_results.tsv'))
  
  return(results)
}



#Plotter---------------------------

plotter <- function(results, output_path, traitA, traitB){
  
  results <- mutate(results, name = str_c(locus, ': ', rsid))
  
  pp4 <- ggplot(data= mutate(results,name= reorder(factor(name),`PP.H4.abf`))) +
    geom_vline(xintercept=0.4, alpha=0.5, linetype='dashed')+
    geom_point(aes(x=`PP.H4.abf`, y=name))+
    geom_segment(aes(x=0, xend=`PP.H4.abf`, y=name, yend=name))+
    
    ylab('Genomic Region: Lead rsID in Trait B')+
    xlab('Posterior Probability for H4: Shared Causal Variant')+
    labs(title=str_wrap(str_c('Posterior probabilities for shared causal variant at genomic regions between traits ', traitA, 'and ', traitB)))+
    guides(colour=guide_legend('Colour', position='bottom'))
  
  print(pp4)
  ggsave(str_c(output_path, '/pp4_summary.png'),pp4,dpi=600)
  
  pp3 <- ggplot(data= mutate(results,name= reorder(factor(name),`PP.H3.abf`))) +
    geom_vline(xintercept=0.4, alpha=0.5, linetype='dashed')+
    geom_point(aes(x=`PP.H3.abf`, y=name))+
    geom_segment(aes(x=0, xend=`PP.H3.abf`, y=name, yend=name))+
    
    ylab('Genomic Region: Lead rsID in Trait B')+
    xlab('Posterior Probability for H3: Distinct Causal Variants @ Shared Locus')+
    labs(title=str_wrap(str_c('Posterior probabilities for distinct causal variants at shared genomic regions between traits ', traitA, 'and ', traitB)))+
    guides(colour=guide_legend('Colour', position='bottom'))
  
  print(pp3)
  ggsave(str_c(output_path, '/pp3_summary.png'),pp3,dpi=600)
  
  results2 <- results %>% 
    mutate(name= reorder(factor(name),`PP.H4.abf`)) %>%
    pivot_longer(contains('PP'), names_to='Hypothesis',values_to='PP') %>%
    mutate(Hypothesis = str_match(Hypothesis, '\\.(H\\d)\\.')[,2])
  
  pp_stacked <- ggplot(data= results2) +
    geom_col(aes(y=name, x=PP, fill=Hypothesis)) +
    geom_vline(xintercept=0.4, alpha=0.5, linetype='dashed')+
    ylab('Genomic Region: Lead rsID in Trait B')+
    xlab('Posterior Probability')+
    labs(title=str_wrap(str_c('Posterior probabilities for each hypothesis at genomic regions between traits ', traitA, ' and ', traitB)))
  
  print(pp_stacked)
  ggsave(str_c(output_path, '/pp_summary.png'),pp_stacked,dpi=600)
  
}

main <- function(results_folderpath){
  
  traitA <- str_match(results_folderpath, '/([^_/]+)_[^_/]+/$')[,2] #Assume traitA_traitB notation
  traitB <- str_match(results_folderpath, '/[^_/]+_([^_/]+)/$')[,2]
  
  results <- importer_merger(results_folderpath,output_path=results_folderpath,traitA, traitB)
  
  plotter(results, results_folderpath, traitA, traitB) #Use the same results folderpath as the output path
  
}

main(results_folderpath)
