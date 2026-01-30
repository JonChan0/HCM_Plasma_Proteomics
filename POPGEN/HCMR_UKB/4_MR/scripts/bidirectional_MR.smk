''' Snakemake script to run Bidirectional Mendelian Randomisation 
Author: Jonathan Chan
Date: 2024-05-01 (ish)
'''
#configfile: '../input/configs/cicp_and_hcm_config.yaml'
rule all:
    input: 
     config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_mainline.tsv',
     tsv_output = config['output_path']+config['d2_exposure']+'_to_'+config['d2_outcome']+'/'+config['d2_exposure']+'_mr_results_mainline.tsv'

rule direction1_mr:
    input:
        exposure_path = config['d1_exposure_path'],
        outcome_path = config['d1_outcome_path']
    params:
        selected_instruments = config['d1_instruments'],
        output_path = config['output_path'],
        exposure = config['d1_exposure'],
        outcome = config['d1_outcome'],
        ld_clumping = config['d1_ld_clump'],
        ld_eur_bed_file = config['ld_eur_bed_file'],
        plink_binary_path = config['plink_binary_path'],
        multivariate= 'FALSE'
    resources:
        mem_mb = 32000
    output:
        tsv_output = config['output_path']+config['d1_exposure']+'_to_'+config['d1_outcome']+'/'+config['d1_exposure']+'_mr_results_mainline.tsv'
    shell:
        '''
        module purge 
        module load R/4.3.2-gfbf-2023a
        Rscript TwoSampleMR.R {input.exposure_path} {input.outcome_path} {params.selected_instruments} {params.output_path} {params.exposure} {params.outcome} {params.ld_clumping} {params.ld_eur_bed_file} {params.plink_binary_path} {params.multivariate}
        '''

rule direction2_mr:
    input:
        exposure_path = config['d2_exposure_path'],
        outcome_path = config['d2_outcome_path']
    params:
        selected_instruments = config['d2_instruments'],
        output_path = config['output_path'],
        exposure = config['d2_exposure'],
        outcome = config['d2_outcome'],
        ld_clumping = config['d2_ld_clump'],
        ld_eur_bed_file = config['ld_eur_bed_file'],
        plink_binary_path = config['plink_binary_path'],
        multivariate = 'FALSE'
    output:
        tsv_output = config['output_path']+config['d2_exposure']+'_to_'+config['d2_outcome']+'/'+config['d2_exposure']+'_mr_results_mainline.tsv'
    resources:
        mem_mb = 32000
    shell:
        '''
        module purge 
        module load R/4.3.2-gfbf-2023a
        Rscript TwoSampleMR.R {input.exposure_path} {input.outcome_path} {params.selected_instruments} {params.output_path} {params.exposure} {params.outcome} {params.ld_clumping} {params.ld_eur_bed_file} {params.plink_binary_path} {params.multivariate}
        '''
