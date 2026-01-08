''' Snakemake script to run FINEMAP pipeline for all traits in output/?_gwas folder
'''

from os import listdir
from os.path import isfile, join
from re import findall

phenotype = config['pheno']
genomic_region = list(range(len(config['chr'].split(','))))

#Define some additional parameters using the phenotype variable
basepath=config['basepath']
output_z_bcor_basepath= "../pipeline_files/"+basepath+"/{pheno}/{pheno}_dataset".format(pheno=phenotype)
output_finemap_basepath= "../output/"+basepath+"/{pheno}/{pheno}_dataset".format(pheno=phenotype)
output_v2g_basepath = "../output/"+basepath+"/{pheno}/v2g/".format(pheno=phenotype)

#If HCMR, need to use manhattan_rsid2.tsv
if 'hcmr' in basepath:
    manhattan_rsid_suffix = '_manhattan_rsid2.tsv'
else:
    manhattan_rsid_suffix = '_manhattan_rsid.tsv'

common_inputs=[
        expand('../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.z',phenotype=phenotype,genomic_region=genomic_region),
        expand('../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset_master.txt',phenotype=phenotype),
        expand('../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset_master2.txt',phenotype=phenotype),
        expand('../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.bcor', phenotype=phenotype, genomic_region=genomic_region),
        expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.snp', phenotype=phenotype, genomic_region=genomic_region),
        expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.config', phenotype=phenotype, genomic_region=genomic_region)
        # ,expand('../output/'+basepath+'/{phenotype}/v2g/{phenotype}_v2g_summary_l2g.tsv', phenotype=phenotype)
        # expand('../output/'+basepath+'/{phenotype}_summstats_topV2G.tsv', phenotype=phenotype)
]

else_inputs = [
    expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset_overallconfig.tsv', phenotype=phenotype) 
]

hcmr_only_inputs=[
    expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset_overallconfig_remapped.tsv', phenotype=phenotype)
]

#If HCMR, defome the additional inputs
if 'hcmr' in basepath:
    all_inputs = common_inputs + hcmr_only_inputs
else:
    all_inputs = common_inputs + else_inputs

rule all:
    input: all_inputs

rule ldstore2_finemap_prep: #Prep the configuration files z file for both LDSTORE2 and FINEMAP
    input: config['manhattan_rsid_folder']+'{phenotype}'+manhattan_rsid_suffix
    output: 
        z='../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.z'
    conda: 
        'gms'
    resources:
        mem_mb=32000
    params:
        pheno=config['pheno'],
        chr=config['chr'],
        base_start=config['base_start'],
        base_end=config['base_end'],
        bgen_filepath=config['bgen_basepath'],
        output_z_bcor_basepath=output_z_bcor_basepath,
        incl_file=config['incl_file'],
        output_finemap_basepath=output_finemap_basepath,
        master_True = 'FALSE',
        sample_filepath = config['sample_filepath'],
        n_samples=config['n_samples']
    shell:'''
        Rscript ldstore2_prep.R {params.pheno} {params.chr} {params.base_start} {params.base_end} {input} {params.bgen_filepath} {params.output_z_bcor_basepath} {params.incl_file} {params.output_finemap_basepath} {params.master_True} {params.sample_filepath} {params.n_samples}
        '''

rule ldstore2_finemap_prep2:#Prep the configuration file of master file for both LDSTORE2 and FINEMAP
    input: rules.ldstore2_finemap_prep.input
    output: 
        ldstore2_master='../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset_master.txt',
        finemap_master='../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset_master2.txt'
    conda: 
        'gms'
    params:
        pheno=config['pheno'],
        chr=config['chr'],
        base_start=config['base_start'],
        base_end=config['base_end'],
        bgen_filepath=config['bgen_basepath'],
        output_z_bcor_basepath=output_z_bcor_basepath,
        incl_file=config['incl_file'],
        output_finemap_basepath=output_finemap_basepath,
        master_True = 'TRUE',
        sample_filepath = config['sample_filepath'],
        n_samples=config['n_samples']
    shell:'''
        Rscript ldstore2_prep.R {params.pheno} {params.chr} {params.base_start} {params.base_end} {input} {params.bgen_filepath} {params.output_z_bcor_basepath} {params.incl_file} {params.output_finemap_basepath} {params.master_True} {params.sample_filepath} {params.n_samples}
        '''


rule ldstore2: #Run LDSTORE2 to compute the SNP correlations for SNPs of interest in the z file i.e per genomic region you are fine-mapping. This computes the LD in the same sample as per their recommended guidelines.
    input: 
        master=rules.ldstore2_finemap_prep2.output.ldstore2_master,
        z=rules.ldstore2_finemap_prep.output.z
    resources:
        mem_mb=32000
    output:
        bcor='../pipeline_files/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.bcor'
    shell:'''
        /well/PROCARDIS/jchan/bin/FINEMAP/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 --write-bcor --read-only-bgen --in-files {input.master}
        '''

rule finemap: #Run finemap
    input: 
        master=rules.ldstore2_finemap_prep2.output.finemap_master,
        z=rules.ldstore2_finemap_prep.output.z,
        bcor=rules.ldstore2.output.bcor
    resources:
        mem_mb=32000
    output:
        snp='../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.snp',
        config='../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.config'
    shell:'''
        /well/PROCARDIS/jchan/bin/FINEMAP/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 --sss --log --in-files {input.master} 
        '''

rule finemap_credible_set_analysis: #This analyses the phenotype overall (all its genomic regions of interest) and remaps the ID to rsID
    input: 
        main=rules.ldstore2_finemap_prep.input,
        snp=expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.snp', phenotype=phenotype, genomic_region=genomic_region),
        config=expand('../output/'+basepath+'/{phenotype}/{phenotype}_dataset{genomic_region}.config', phenotype=phenotype, genomic_region=genomic_region)
    output: 
        main = '../output/'+basepath+'/{phenotype}/{phenotype}_dataset_overallconfig_remapped.tsv' if 'hcmr' in basepath else '../output/'+basepath+'/{phenotype}/{phenotype}_dataset_overallconfig.tsv' #This remaps the bgen chr:pos ID in the HCMR BGEN files and thus in the overallconfig.tsv files to rsID or at least snid (chr:pos_alleleA_alleleB)
    conda: 
        'gms'
    resources:
        mem_mb=32000
    params:
        pheno=config['pheno'],
        output_finemap_basepath=output_finemap_basepath,
        basepath=basepath,
        manhattan_rsid2_basepath= config['manhattan_rsid_folder']
    shell:'''
        Rscript finemap_analysis.R {params.pheno} {params.output_finemap_basepath} {input.main} {output.main}

        if [[ {params.basepath} =~ 'hcmr' ]]; then
            Rscript hcmr_rsid_remapper.R {params.pheno} {params.output_finemap_basepath} {params.manhattan_rsid2_basepath}
        fi
    '''

# rule v2g_otargen: #This doesn't actually run in rescomp because it doesn't have internet access so can't call API
#     input: 
#         main= rules.finemap_credible_set_analysis.output.main
#     output:
#         main='../output/'+basepath+'/{phenotype}/v2g/{phenotype}_v2g_summary_l2g.tsv'
#     conda:
#         'R_4.3.3'
#     params:
#         pheno=config['pheno'],
#         output_v2g_basepath=output_v2g_basepath,
#         l2g_score_cutoff=0.1,
#         Rmd_template_path='v2g_otargen_tables_template.Rmd'
#     shell:
#         '''
#         Rscript v2g_otargen.R {input.main} {params.pheno} {params.output_v2g_basepath} {params.l2g_score_cutoff} {params.Rmd_template_path}
#         '''

#Final step is running the summstat_finemap_v2g_summariser.R which enables plotting of overalls summaries from GWAS; FINEMAP and V2G
# rule summstat_finemap_v2g_summariser:
#     input:
#         main='../output/'+basepath+'/{phenotype}/v2g/{phenotype}_v2g_summary_l2g.tsv'
#     output:
#         main='../output/'+basepath+'/{phenotype}_summstats_topV2G.tsv'
#     conda:
#         'gms'
#     resources: 
#         mem_mb=32000
#     params:
#         pheno=config['pheno'],
#         basepath=basepath,
#         summstats_basepath= config['manhattan_rsid_folder'],
#         chainfilepath = '/well/PROCARDIS/jchan/bin/liftover/hg38ToHg19.over.chain'
#     shell:'''
#         Rscript summstat_finemap_v2g_summariser.R {params.pheno} {params.basepath} {params.summstats_basepath} {params.chainfilepath}
#     '''