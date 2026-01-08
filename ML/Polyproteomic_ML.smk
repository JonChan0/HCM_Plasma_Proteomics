'''
Snakemake pipeline for investigating ML models for using plasma proteomic data for predictive utility
Author: Jonathan Chan
Date: 2025-01-23
'''

import os

rule all:
    input:
        model_pkl_files = expand(config['base_model_folder'] + config['folder_suffix']+"{model_name}_"+config['filesuffix']+"_best_model.pkl", model_name=config['model_names']),
        cv_results = expand(config['base_model_folder'] + config['folder_suffix']+"{model_name}_"+config['filesuffix']+"_cv_results.csv", model_name=config['model_names']),
        cv_crossmodel_comp = config['base_plot_folder'] + 'cv_performance/crossmodel_comparison/model_comparison.png'

#This outlines the training over the various different model types to compare different model classes 
rule train:
    input:
        X_train_data_path=config['base_data_folder'] + config['folder_suffix'] + config['X_train_data_filename'],
        y_train_data_path=config['base_data_folder'] + config['folder_suffix'] + config['y_train_data_filename']
    output:
        model_pkl_file = config['base_model_folder'] + config['folder_suffix']+"{model_name}_"+config['filesuffix']+"_best_model.pkl",
        cv_results = config['base_model_folder'] + config['folder_suffix']+"{model_name}_"+config['filesuffix']+"_cv_results.csv"
    threads: 8
    conda: 'python3.11_ml'
    resources:
        mem_mb = lambda wildcards, threads: threads * 8000  # Example: 4000 MB per thread
    params:
        model_output_folder=config['base_model_folder'] + config['folder_suffix'],
        plot_output_folder=lambda wildcards: config['base_model_folder'] + config['folder_suffix'] + config['model_output_names'][config['model_names'].index(wildcards.model_name)],
        features_to_bypass=config['features_to_bypass'], #This refers to the quantitative features you should NOT pass to feature selection step
        features_to_select=config['features_to_select'],  #This refers to the quantitative features you SHOULD pass to feature selection step
        target_variable=config['target_variable'],
        feature_selection=lambda wildcards: config['feature_selection'][config['model_names'].index(wildcards.model_name)],
        filesuffix=config['filesuffix']
    shell:'''
        # Set the number of threads for Python
        export OMP_NUM_THREADS=8

        python 1_Polyproteomic_ML_Train.py \
        --model_name {wildcards.model_name} \
        --plot_output_folder "{params.plot_output_folder}" \
        --model_output_folder "{params.model_output_folder}" \
        --feature_selection {params.feature_selection} \
        --features_to_bypass_fs {params.features_to_bypass} --features_to_select_fs {params.features_to_select} \
        --target_variable {params.target_variable} \
        --X_train_data {input.X_train_data_path} --y_train_data {input.y_train_data_path} \
        --filesuffix {params.filesuffix}
    '''

#This compares the validation-fold performance of the various model types over one script run

rule cv_performance_comparison:
    input:
        expand(config['base_model_folder'] + config['folder_suffix']+"{model_name}_"+config['filesuffix']+"_cv_results.csv", model_name=config['model_names'])
    output:
        cv_crossmodel_comp = config['base_plot_folder'] + 'cv_performance/crossmodel_comparison/model_comparison.png'
    threads: 2
    conda: 'R_4.3.3'
    resources:
        mem_mb = lambda wildcards, threads: threads * 8000  # Example: 4000 MB per thread
    params:
        input_folder = config['base_model_folder'] + config['folder_suffix']+'/',
        plot_output_folder=config['base_plot_folder'] + 'cv_performance/',
        model_class=config['model_class_of_interest']
    shell:'''
        # Set the number of threads for Python
        export OMP_NUM_THREADS=2

        Rscript 2_Polyproteomic_CV_Performance_Plotter.R {params.input_folder} {params.plot_output_folder} {params.model_class}
    '''


#This outputs the SHAP values of feature importance for a chosen model class of interest
# rule shap:
#     input:
#         model_pkl = rules.train.output.model_pkl_file,
#         X_train_preprocessed = config['base_data_folder'] + config['folder_suffix'] + 'X_train_preprocessed_' + '{model_name}_'+ config['folder_suffix'] + '.csv',
#         y_train = config['base_data_folder'] + config['folder_suffix'] + config['y_train_data_filename']
#     output:
#         config['base_plot_folder'] + 'feature_importance/' + config['folder_suffix'] + {wildcards.model_name}+'/shap/'
#     threads: 8
#     conda: 'python3.11_ml'
#     params:
#         plot_out  = config['base_plot_folder'] + 'feature_importance/' + config['folder_suffix'] + {wildcards.model_name}+'/shap/',
#         model_name = {wildcards.model_name},
#         pp_names = config['allpp_names_csv'],
#         bootstrap = config['shap_bootstrap']
#     shell:'''
#         python 3_Polyproteomic_SHAP.py \
#           --model_pkl_file {input.model_pkl} \
#           --X_train_preprocessed_path {input.X_train_preprocessed} \
#           --y_train_path {input.y_train} \
#           --plot_output_path {params.plot_out} \
#           --model_name {params.model_name} \
#           --pp_names_file {params.pp_names} \
#           --bootstrap_or_not {params.bootstrap}
#         '''