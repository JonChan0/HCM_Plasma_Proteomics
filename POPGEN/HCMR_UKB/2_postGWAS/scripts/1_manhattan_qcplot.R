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
library(ggmanh)

args <- commandArgs(trailingOnly=TRUE) #Allows taking of arguments in bash command line #By default you should pass the path (relative to the script)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
chunk_file <- args

input <- args[1]
output_plot_path <- args[2]
pvalue_format <- args[3]

theme_set(
  theme_classic()
)

#Local test code
#input <- 'popgen/2_GWAS/output/cicp_bap_manhattan_rsid.tsv'
# output_path <- 'popgen/2_GWAS/output/pipeline_files/'
# output_plot_path <- 'popgen/2_GWAS/output/plots/manhattan/'

#ForUKB
# input <- 'popgen/2_gwas/output/ukb/REGENIE/step2/formatted/NTproBNP_manhattan_rsid.tsv'
# output_path <- 'popgen/2_gwas/output/pipeline_files/nonHarperUKB_selectpp/'
# output_plot_path <- 'popgen/2_gwas/output/plots/nonHarperUKB_selectpp/'
# pvalcutoff <- 0.00000005

#-------------------------------------------------------------------------------
#Import

pheno <- str_match(input, '/([^/]+)_manhattan_rsid')[,2] #This assumes the format of output/phenoname_manhattan
import <- read_tsv(input, progress=T)

#Reformat the log-pvalues to raw pvalue if only in log form

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

#----------------------------------------------------------------------
#Output Manhattan plot

manhattan <- manhattan_plot(import, chr.colname = "chromosome", pos.colname = "position", 
                            plot.title = str_c('Manhattan plot for phenotype ',pheno), y.label = "-log10(pval)",
                            rescale=T,
                            #, label.colname='Label'
)

#print(manhattan)
ggsave(str_c(output_plot_path,pheno,'_manhattan.png'),manhattan,dpi=600,width=12, height=6)
rm(manhattan)

#----------------------------------------------------------------------------------
#Additional QC e.g QQ plots; p-value histograms

#Compute the genomic inflation factor
chisq <- qchisq(import$pval,1, lower.tail=F)
#qchisq(assoc.df$P,1,lower.tail=FALSE) can convert p-value (even < 5.5e-17) to chisq
# while qchisq(1-assoc.df$P,1) fails to convert small p-value (Inf in this case) 
#As per https://bioinformaticsngs.wordpress.com/2016/03/08/genomic-inflation-factor-calculation/
lambdagc <- round((median(chisq)/qchisq(0.5,1)), digits=3)


#Compute the QQ
qq <- qqunif(import$pval)+
  labs(title=str_wrap(str_c('QQplot for ', pheno, ' GWAS results')),
       caption=str_c('Genomic Inflation Factor = ',as.character(lambdagc)))
ggsave(str_c(output_plot_path,pheno,'_qq.png'),qq,dpi=600,width=9, height=6)

#P-value histogram
pvalhist <- ggplot(import)+
  geom_histogram(aes(x=pval), bins = 100)+
  scale_x_continuous(limits = c(0,1))+
  theme_classic()+
  ylab('Frequency')+
  xlab('p-value')+
  labs(title=str_wrap(str_c('p-value histogram for ', pheno, ' GWAS results')),
       caption=str_c('Genomic Inflation Factor: ', as.character(lambdagc)))

ggsave(str_c(output_plot_path,pheno,'_pvalhist.png'),pvalhist,dpi=600,width=9, height=6)

rm(qq, pvalhist)

message("Script completed successfully.")




