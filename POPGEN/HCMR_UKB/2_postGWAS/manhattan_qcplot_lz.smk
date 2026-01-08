'''
Snakemake pipeline to perform Manhattan and other QC plotting + LocusZoom plotting
Author: Jonathan Chan
Date: 2024-05-06

Input: ****_manhattan_rsid.tsv
Output: ****_manhattan & other. pngs
'''

from os import listdir
from os.path import isfile, join
from re import findall

#configfile: 'config.yaml'

tsv_basepath = config['tsv_basepath']
phenos = [f for f in listdir(tsv_basepath) if isfile(join(tsv_basepath, f))]
phenotype_list_list = [findall('(.+)_manhattan_rsid.tsv', p) for p in phenos] 
phenotype = [item for sublist in phenotype_list_list for item in sublist] 

#Define conditional outputs
outputs=[expand(config['plot_output_path'] + '{pheno}/' + '{pheno}_manhattan.png', pheno=phenotype),
    expand(config['plot_output_path'] + '{pheno}/' + '{pheno}_qq.png', pheno=phenotype),
    expand(config['plot_output_path'] + '{pheno}/' + '{pheno}_pvalhist.png', pheno=phenotype)
    ,
    expand(config['output_path'] + 'lz_params/'+'{pheno}_lz_params.tsv', pheno=phenotype),
    expand(config['plot_output_path'] + '{pheno}/' +'lz/', pheno=phenotype)
    ]

if config['yaml_output'] == 'True':
    outputs = outputs+ [expand(config['output_yaml_path'] + '{pheno}_config.yaml', pheno=phenotype)]


rule all:
    input:
        outputs

rule manhattan_qcplot:
    input: 
        manhattan_tsv = tsv_basepath + '{pheno}_manhattan_rsid2.tsv' if 'hcmr' in config['basepath'] else tsv_basepath + '{pheno}_manhattan_rsid.tsv'
    output:
        manhattan_plot = config['plot_output_path'] + '{pheno}/' + '{pheno}_manhattan.png',
        qqplot = config['plot_output_path'] + '{pheno}/' + '{pheno}_qq.png',
        pvalhist = config['plot_output_path'] + '{pheno}/' + '{pheno}_pvalhist.png'
    resources:
        mem_mb=64000
    params:
        output_plot_path = config['plot_output_path'] + '{pheno}/',
        pvalue_format=config['pvalue_format']
    shell:'''
        module purge 
        module use -a /apps/eb/2022b/skylake/modules/all
        module load R/4.2.2-foss-2022b
        Rscript scripts/1_manhattan_qcplot.R {input.manhattan_tsv} {params.output_plot_path} {params.pvalue_format}
    '''
#Note that this step hangs on SLURM so should not be run! - You need to run locally
# rule lz_params_writing:
#     input:
#         manhattan_tsv =  tsv_basepath + '{pheno}_manhattan_rsid2.tsv' if 'hcmr' in config['basepath'] else tsv_basepath + '{pheno}_manhattan_rsid.tsv',
#         lz_tsv = tsv_basepath + '{pheno}_lz.tsv'
#     output:
#         lzparams = config['output_path'] + 'lz_params/' + '{pheno}_lz_params.tsv'
#     resources:
#         mem_mb=32000
#     params:
#         output_path = config['output_path'] + 'lz_params/' ,
#         pvalcutoff = config['pvalcutoff'],
        #   online_formatting = config['online_formatting']
#     shell:'''
#         module purge 
#         module use -a /apps/eb/2022b/skylake/modules/all
#         module load R/4.2.2-foss-2022b
#         Rscript scripts/2_lzparams_writer.R {input.manhattan_tsv} {input.lz_tsv} {params.output_path} {params.pvalcutoff} {params.online_formatting}
#     '''

rule locuzoom:
    input:
        config['output_path'] + 'lz_params/'+'{pheno}_lz_params.tsv'
    output:
        directory(config['plot_output_path'] + '{pheno}/' +'lz/')
    resources:
        mem_mb=16000
    conda: 'locuszoom'
    params:
        output_plot_path = config['plot_output_path'] + '{pheno}/',
        lz_plotting_script = config['lz_plot_script']
    shell:'''
        mkdir {params.output_plot_path}'/lz/'
        cd {params.output_plot_path}'/lz/'

        tail -n +2 {input}| while IFS=$'\t' read -r line; do
        {params.lz_plotting_script} $(echo "$line" | cut -f1-5)
        done
    '''

if config['yaml_output'] == 'True':
    rule yaml_output: #This writes out the .yaml file required as input into FINEMAP
        input:
            lz_params = config['output_path'] + 'lz_params/' + '{pheno}_lz_params.tsv'
        output:
            yaml_output = config['output_yaml_path'] + '{pheno}_config.yaml'
        resources:
            mem_mb=32000
        params:
            output_yaml_path = config['output_yaml_path'],
            bgen_filepath=config['bgen_filepath'],
            incl_filepath=config['incl_filepath'],
            sample_filepath=config['sample_filepath'],
            tsv_basepath=config['tsv_basepath'],
            n_samples=config['n_samples'],
            basepath=config['basepath']
        shell:'''
            module purge 
            module use -a /apps/eb/2022b/skylake/modules/all
            module load R/4.2.2-foss-2022b
            Rscript scripts/3_lz_yaml_writer.R {input.lz_params} {params.output_yaml_path} {params.bgen_filepath} {params.incl_filepath} {params.sample_filepath} {params.tsv_basepath} {params.n_samples} {params.basepath}
        '''