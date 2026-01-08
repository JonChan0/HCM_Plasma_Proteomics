#example
#./run.bash locusname 22 start end CHD_SAIGE_30_1_19.metal
###module use /mgmt/modules/eb/modules/all/
# module purge
# module use -a /apps/eb/2022b/skylake/modules/all
# module load python/2.7.11
# module load htslib/1.8-gcc5.4.0
# module load R/3.2.2
#module load R/3.5.1-foss-2018b

region=$1
chr=$2
st=$3
en=$4
metal=$5
echo $region $metal $chr
#If want to use Phase1 V3 VCF files
#--ld-vcf ../vcf/phase1v3/chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz \
#If want to use Phase3 VCF files
#--ld-vcf /well/PROCARDIS/bin/locuszoom/vcf/phase1v3/chr$chr\.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz \

#If you want to use UKBB imputed .vcf for LD calculation
# --ld-vcf /well/PROCARDIS/bin/locuszoom/vcf/ukbb_v3_maf_imp_$chr\_eur2.vcf.gz \

/well/PROCARDIS/bin/locuszoom/bin/locuszoom --metal $metal \
                 --ld-vcf /well/PROCARDIS/bin/locuszoom/vcf/phase1v3/chr$chr\.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.vcf.gz \
                 title=$region recombColor=orange recombOver=TRUE signifLine="5,7.3" signifLineColor="blue,red" \
                 --build hg19 \
                 --cache None \
                 --plotonly \
                 --chr $chr \
                 --start $st \
                 --end $en  \
                 --snpset NULL \
                 --prefix $region \
                 >& $region.$chr.$st.log
echo "Done!"

#The command below works
#If you want to do any edits, do it on the (uncommented) command above
#../bin/locuszoom --metal $metal --refsnp $snp --ld-vcf ../vcf/ukbb_v3_maf_imp_$chr\_eur2.vcf.gz recombColor=orange recombOver=TRUE signifLine="5,7.3" signifLineColor="blue,red" --build hg19 --cache None --plotonly --snpset NULL

