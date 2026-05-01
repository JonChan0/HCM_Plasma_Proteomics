'''
Snakemake script to run REGENIE for GWAS in UKB and HCMR
'''

#configfile: "config.yaml"#Define this when you run the snakemake.

rule all:
    input:
        config['step1_output_path']+config['output_prefix']+"_pred.list", #The output from the first step of REGENIE
        expand("{path}{output_prefix}_{chr}.regenie",path=config['step2_output_path'],output_prefix=config['output_prefix'],chr=config['chromosomes']), #The output from the second step of REGENIE
        expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']),
        expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid_logpval.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']), #The output from the collation of all chromosome's REGENIE results and split into per-phenotype manhattan_rsid.tsv
        expand("{path}formatted/{output_prefix}_{pheno}_lz.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']), #The per-phenotype lz.tsv
        # expand("{path}formatted/{pheno}_1Eneg5.tsv",path=config['step2_output_path'],pheno=config['phenotypes']), #The per-phenotype tsv with pvalues < 1e-5
        expand("{path}formatted/{output_prefix}_{pheno}_ldsc.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']) #The per-phenotype tsv for LDSC

if config["run_hardcalled_prep"]:
    rule hardcalled_snps_prep:
        input:
            hardcalled_snps_basepath=config['hardcalled_genotypes_basepath']
        output:
            hardcalled_snps_allchr_bedfile=config['hardcalled_genotypes_allchr_merge_path']+'.bed',
            hardcalled_snps_allchr_bimfile=config['hardcalled_genotypes_allchr_merge_path']+'.bim',
            hardcalled_snps_allchr_famfile=config['hardcalled_genotypes_allchr_merge_path']+'.fam'
        params:
            hardcalled_snps_fam_path=config['hardcalled_genotypes_fam'],
            hardcalled_snps_allchr_merge_path=config['hardcalled_genotypes_allchr_merge_path'],
            mergelist_txt=config['mergelist_file'],
            starter_bed=config['hardcalled_merged_starter_bed'],
            starter_bim=config['hardcalled_merged_starter_bim']
        resources:
            mem_mb=64000
        shell:'''
            module purge
            module load PLINK/1.9b_6.21-x86_64

            plink \
            --bed {input.hardcalled_snps_basepath}/{params.starter_bed} \
            --bim {input.hardcalled_snps_basepath}/{params.starter_bim} \
            --fam {params.hardcalled_snps_fam_path} \
            --merge-list {params.mergelist_txt} \
            --make-bed --out {params.hardcalled_snps_allchr_merge_path}
        '''

    rule hardcalled_snps_qc:
        input: 
            rules.hardcalled_snps_prep.output.hardcalled_snps_allchr_bedfile,
            rules.hardcalled_snps_prep.output.hardcalled_snps_allchr_bimfile,
            rules.hardcalled_snps_prep.output.hardcalled_snps_allchr_famfile
        output:
            hardcalled_snps_allchr_qc_bedfile=config['hardcalled_genotypes_allchr_qc_path']+'.bed',
            hardcalled_snps_allchr_qc_bimfile=config['hardcalled_genotypes_allchr_qc_path']+'.bim',
            hardcalled_snps_allchr_qc_famfile=config['hardcalled_genotypes_allchr_qc_path']+'.fam'
        params:
            hardcalled_snps_allchr_merge_path=config['hardcalled_genotypes_allchr_merge_path'],
            hardcalled_snps_allchr_qc_path=config['hardcalled_genotypes_allchr_qc_path'],
            ukb_redacted_ids=config['ukb_redactedIDs']
        resources:
            mem_mb=32000
        shell:'''
            module purge
            module load PLINK/1.9b_6.21-x86_64

            plink \
            --bfile {params.hardcalled_snps_allchr_merge_path} \
            --maf 0.01 \
            --mac 100 \
            --geno 0.01 \
            --hwe 1e-15 \
            --mind 0.1 \
            --remove {params.ukb_redacted_ids} \
            --make-bed --out {params.hardcalled_snps_allchr_qc_path}
        '''

    rule hardcalled_snps_ldprune:
        input:
            rules.hardcalled_snps_qc.output.hardcalled_snps_allchr_qc_bedfile,
            rules.hardcalled_snps_qc.output.hardcalled_snps_allchr_qc_bimfile,
            rules.hardcalled_snps_qc.output.hardcalled_snps_allchr_qc_famfile
        output:
            hardcalled_snps_allchr_ldpruned_bedfile=config['hardcalled_genotypes_allchr_ldpruned2_path']+'.bed',
            hardcalled_snps_allchr_ldpruned_bimfile=config['hardcalled_genotypes_allchr_ldpruned2_path']+'.bim',
            hardcalled_snps_allchr_ldpruned_famfile=config['hardcalled_genotypes_allchr_ldpruned2_path']+'.fam'
        params:
            hardcalled_snps_allchr_qc_path=config['hardcalled_genotypes_allchr_qc_path'],
            hardcalled_snps_allchr_ldpruned1_path=config['hardcalled_genotypes_allchr_ldpruned1_path'],
            hardcalled_snps_allchr_ldpruned2_path=config['hardcalled_genotypes_allchr_ldpruned2_path']
        resources:
            mem_mb=32000
        shell:'''
            module purge
            module load PLINK/1.9b_6.21-x86_64

            plink \
            --bfile {params.hardcalled_snps_allchr_qc_path} \
            --indep-pairwise 1000 100 0.8 \
            --out {params.hardcalled_snps_allchr_ldpruned1_path}

            plink \
            --bfile {params.hardcalled_snps_allchr_qc_path} \
            --extract {params.hardcalled_snps_allchr_ldpruned1_path}.prune.in \
            --make-bed --out {params.hardcalled_snps_allchr_ldpruned2_path}
        '''

rule regenie_step1:
    input:
        config['hardcalled_genotypes_allchr_ldpruned2_path']+'.bed',
        ldpruned_hardcalled_snps=rules.hardcalled_snps_ldprune.output.hardcalled_snps_allchr_ldpruned_bedfile if config["run_hardcalled_prep"] else config['hardcalled_genotypes_allchr_ldpruned2_path']+'.bed', #Induce a dependence on the previous step
        phenoFile=config['pheno_file'],
        covarFile=config['covar_file'],
        sampleInclusionFile=config['sampleInclusion_file']
    output:
        map_file=config['step1_output_path']+config['output_prefix']+"_pred.list"
    resources:
        mem_mb=32000
    threads: 8
    conda:
        'regenie3.4.1_env'
    params:
        hardcalled_snps_allchr_ldpruned2_path=config['hardcalled_genotypes_allchr_ldpruned2_path'],
        nThreads=8,
        output_prefix=config['step1_output_path']+config['output_prefix'],
        sampleInclusion_file=config['sampleInclusion_file'], 
        optional_params=config['step1_optional_params']
    shell:'''
        regenie \
        --step 1 \
        --bed {params.hardcalled_snps_allchr_ldpruned2_path} \
        --bsize 1000 \
        --keep {params.sampleInclusion_file} \
        --phenoFile {input.phenoFile} \
        --covarFile {input.covarFile} \
        --loocv \
        --apply-rint \
        --threads {params.nThreads} \
        --out {params.output_prefix} \
        {params.optional_params}
    '''

rule regenie_step2_singleAssoc:
    input:
        phenoFile=config['pheno_file'],
        covarFile=config['covar_file'],
        step1_map_file=rules.regenie_step1.output.map_file,
        bgen_file=config['imputed_genotypes_bgens']+"{chr}"+config['imputed_genotypes_bgen_suffix'],
        sample_file=config['imputed_genotypes_sample']
    output:
        config['step2_output_path']+config['output_prefix']+"_{chr}.regenie"
    resources:
        mem_mb=32000
    threads: 8
    conda:
        'regenie3.4.1_env'
    params:
        nThreads=8,
        output_prefix=config['step2_output_path']+config['output_prefix']+"_{chr}",
        sampleInclusion_file=config['sampleInclusion_file']
    shell:'''
        regenie \
        --step 2 \
        --bgen {input.bgen_file} \
        --ref-first \
        --sample {input.sample_file} \
        --covarFile {input.covarFile} \
        --phenoFile {input.phenoFile} \
        --keep {params.sampleInclusion_file} \
        --minINFO 0.7 \
        --minMAC 50 \
        --bsize 1000 \
        --apply-rint \
        --pred {input.step1_map_file} \
        --out {params.output_prefix} \
        --no-split \
        --threads {params.nThreads}
    '''

rule chr_merger_pheno_extractor:
    input:
        expand("{path}{output_prefix}_{chr}.regenie",path=config['step2_output_path'],output_prefix=config['output_prefix'],chr=config['chromosomes']) #The output from the second step of REGENIE
    output:
        expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix']),
        expand("{path}formatted/{output_prefix}_{pheno}_manhattan_rsid_logpval.tsv",path=config['step2_output_path'],pheno=config['phenotypes'],output_prefix=config['output_prefix'])
    params:
        regenie_path = config['step2_output_path'],
        output_prefix=config['output_prefix']
    resources:
        mem_mb=64000
    shell:'''
        module purge
        module load R/4.3.2-gfbf-2023a

        Rscript REGENIE_collater_formatter.R {params.regenie_path} {params.output_prefix}
    '''

if config['rsid_mapping']:
    rule snid_to_rsid_mapper:
        input:
            main=config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid.tsv",
            logpval=config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid_logpval.tsv"
        output:
            main=config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid2.tsv",
            temp_snid = temp(config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_snid.tsv"),
            temp_snid_logpval = temp(config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_snid_logpval.tsv"),
            logpval=config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid2_logpval.tsv"
        params:
            rsid_mapping_script = config['rsid_mapping_script']
        resources:
            mem_mb=64000
        shell:'''
            awk '{{$14=$2"_"$5"_"$6}} 1' {input.main} | awk '{{ $1=""; sub(/^ /, ""); print }}' | awk 'BEGIN {{FS=" "; OFS="\t" }} {{$1=$1; print}}' > {output.temp_snid}

            module load Perl/5.36.1-GCCcore-12.3.0
            export PERL5LIB=/well/PROCARDIS/bin/perl_modules/lib/perl5/site_perl/5.28.1/:$PERL5LIB 
            perl {params.rsid_mapping_script} {output.temp_snid} {output.main} 12

            awk '{{$14=$2"_"$5"_"$6}} 1' {input.logpval} | awk '{{ $1=""; sub(/^ /, ""); print }}' | awk 'BEGIN {{FS=" "; OFS="\t" }} {{$1=$1; print}}' > {output.temp_snid_logpval}
            perl {params.rsid_mapping_script} {output.temp_snid_logpval} {output.logpval} 12
        '''
#Note that the column number of perl is 0-based index so 12 refers to 13th column (14-1 when you remove the first column)

rule lz_tsv_generator:
    input: 
        config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid2.tsv" if config['rsid_mapping'] else config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid.tsv"
    output:
        config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_lz.tsv"
    resources:
        mem_mb=8000
    shell:'''
        awk 'BEGIN {{ print "MarkerName\tP-value" }} NR > 1 {{ if (sub(/_.*/, "", $1))$1 = "chr"$1; print $1 "\t" $11 }}' {input} > {output}
    '''

# rule Eneg5_filter:
#     input:
#         config['step2_output_path']+"formatted/{pheno}_manhattan_rsid2_logpval.tsv" if config['rsid_mapping'] else config['step2_output_path']+"formatted/{pheno}_manhattan_rsid_logpval.tsv"
#     output:
#         temp=temporary(config['step2_output_path']+"formatted/{pheno}_gws_temp.tsv"),
#         gws_out=config['step2_output_path']+"formatted/{pheno}_1Eneg5.tsv"
#     resources:
#         mem_mb=8000
#     shell:'''
#         awk '{{if($11 >= 5) print}}' {input} > {output.temp}
#         (head -n 1 {input} && cat {output.temp}) > {output.gws_out}
#     '''

rule ldsc_file_generator:
    input: 
        config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid2.tsv" if config['rsid_mapping'] else config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_manhattan_rsid.tsv" #Effect allele = column 6 (Allele_B) which needs to map onto A1 (for LDSC.tsv)
    output: config['step2_output_path']+"formatted/"+config['output_prefix']+"_{pheno}_ldsc.tsv"
    shell:'''
	    awk 'BEGIN {{ print "SNP\tCHR\tPOS\tA1\tA2\tINFO\tN\tEAF\tP\tBETA\tBETA_SE" }} NR > 1 {{$1 = ($1 ~ /_/) ? $13 : $1; print $1,"\t",$3,"\t",$4,"\t",$6,"\t",$5,"\t",$7,"\t",$8,"\t",$9,"\t",$11,"\t",$12,"\t",$13 }}' {input} > {output}
    '''