#!/bin/bash
#SBATCH -p long
#SBATCH -A procardis.prj
#SBATCH --constraint="skl-compat"
#SBATCH -c 4

date

module purge

#Required for RSid mapping script
module load Perl/5.28.1-GCCcore-8.2.0 
export PERL5LIB=/well/PROCARDIS/bin/perl_modules/lib/perl5/site_perl/5.28.1/:$PERL5LIB 

#Define the file to map to RSID
basepath='../input/gwas_summary_statistics/HCM/'
name="$basepath"hcm_meta.230523.fix.gwama.nohcmr

input_path="$name".txt.out
output_path="$name"_rsid.tsv
filtered_path="$name"_filtered.tsv

#Executing code to map snptestID to RSID 
echo 'Appending the snptestID'
awk '{$18=substr($1,4)} 1' $input_path | awk 'NR>1' | awk 'BEGIN{print "rs_number reference_allele other_allele eaf beta se beta_95L beta_95U z p-value _-log10_p-value q_statistic q_p-value i2 n_studies n_samples effects snptestid"}1'> $basepath/test.tsv
tr ' ' '\t' < $basepath/test.tsv > $basepath/test2.tsv
rm $basepath/test.tsv

echo "Adding on the rsid column using the snptestID"
perl "map_snptest_rsid.pl" $basepath/test2.tsv $output_path 17
rm $basepath/test2.tsv

echo 'Filtering for variants with pval < 10^-5'
awk '{if ($12 >=5) print $0}' $output_path > $filtered_path
