'''
Script to run colocalisation as a sensitivity analysis for +ve results from Mendelian randomisation
Author: Jonathan Chan
Date: 2024-10-07
'''

if config['susie'] == 'TRUE':
    coloc_output = expand(config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+'{locusname}/susie_{rsid}_{chr}_{lead_snp_position}_summaryresults.txt', zip,locusname=config['locusname'],rsid=config['rsid'],chr=config['chr'],lead_snp_position=config['lead_snp_position']) 
else:
    coloc_output = coloc_output =   expand(config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+'{locusname}/abf_{rsid}_{chr}_{lead_snp_position}_summaryresults.txt', zip,locusname=config['locusname'],rsid=config['rsid'],chr=config['chr'],lead_snp_position=config['lead_snp_position']) 

rule all:
    input: 
        coloc_output, 
        config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+config['traitname_A']+'_'+config['traitname_B']+'_summary_results.tsv'


#ARCHIVED!
# rule ld_generator: #This generates the LD file using plink from the QCed hard-called genotype files from UKB - deprecated because needs to be LD-specific
#     input:
#         bgen=ukb_config['ukb_imputed_path'] + '{chr}' + '_v3.bgen',
#         sample=config['ukb_imputed_sample'],
#         chr='{chr}',
#         region_start='{lead_snp_position}'-500000,
#         region_end = '{lead_snp_position}'+500000
#     output:
#         ld_info = config['output_path']+'/LD/'+'chr{chr}_' +'{lead_snp_position}_{locusname}_ld.bin'
#     params:
#         output_prefix= config['output_path']+'/LD/'+'chr{chr}_' +'{lead_snp_position}_{locusname}'
#     resources:
#         mem_mb = 64000
#     shell:
#         '''
#         module purge 
#         /well/PROCARDIS/jchan/bin/plink2 \
#             --memory {resources.mem_mb} \
#             --bgen {input.bgen} ref-first --sample {input.sample} \
#             --chr {input.chr} --from-bp {input.region_start} --to-bp {input.region_end} \
#             --mach-r2-filter 0.7 \
#             --r2-phased square bin \
#             --out {params.output_prefix}
#         '''

rule coloc:
    input:
        summstatA = config['pathA'],
        summstatB = config['pathB']
    output: 
        coloc_output = config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+'{locusname}/susie_{rsid}_{chr}_{lead_snp_position}_summaryresults.txt' if config['susie'] == 'TRUE' else config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+'{locusname}/abf_{rsid}_{chr}_{lead_snp_position}_summaryresults.txt'
    params:
        chr='{chr}',
        lead_snp_position='{lead_snp_position}',
        locusname = '{locusname}',
        rsid='{rsid}',
        traittype_A = config['traittype_A'],
        traittype_B = config['traittype_B'],
        traitsd_A = config['traitsd_A'],
        traitsd_B = config['traitsd_B'] ,
        output_path = config['output_path'],
        traitname_A = config['traitname_A'],
        traitname_B = config['traitname_B'],
        N_A = config['N_A'],
        N_B = config['N_B'],
        ld_path = config['ld_path'],
        susie=config['susie']
    resources:
        mem_mb=32000
    shell:
        '''
        module purge
        module load R/4.3.2-gfbf-2023a
        Rscript coloc_run.R {input.summstatA} {input.summstatB} {params.chr} {params.lead_snp_position} {params.rsid} {params.locusname} {params.traittype_A} {params.traittype_B} {params.traitsd_A} {params.traitsd_B} {params.output_path} {params.traitname_A} {params.traitname_B} {params.N_A} {params.N_B} {params.ld_path} {params.susie}
        '''

rule coloc_summariser: #This summarises across all genomic regions for a single pair of traits
    input: 
        rules.coloc.output.coloc_output,
        results_dir = dir(config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/')
    output:
        summary_results_tsv = config['output_path']+'plots/'+config['traitname_A']+'_'+config['traitname_B']+'/'+config['traitname_A']+'_'+config['traitname_B']+'_summary_results.tsv'
    resources:
        mem_mb = 8000
    shell:
        '''
        module purge
        module load R/4.3.2-gfbf-2023a
        Rscript coloc_summariser.R {input.results_dir}
        '''