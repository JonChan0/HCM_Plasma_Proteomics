#!/bin/bash
#SBATCH -p long
#SBATCH -A procardis.prj
#SBATCH --constraint="skl-compat"
#SBATCH -c 3
date

module load PLINK/1.9b_6.21-x86_64

#Define your input details here
chr=8
snp=rs2701918

plink \
	--r2 \
	--vcf /well/PROCARDIS/bin/locuszoom/vcf/ukbb_v3_maf_imp_$chr\_eur2.vcf.gz \
	--ld-window-r2 0.8 \
	--ld-snp $snp

date
