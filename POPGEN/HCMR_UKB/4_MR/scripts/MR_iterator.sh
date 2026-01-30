#This script iterates over all the config files in ../input/configs to run the mendelian_randomisation.smk pipeline for each pair of traits

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate gms

echo $(date)

dataset=$1

if [ $dataset == 'hcmr' ]
then

    #This applies to running the bidirectional MR with HCMR GWAS of plasma proteins with HCM
    for file in ../input/configs/hcmr/*.yaml
    do
        echo "Carrying out bidirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile bidirectional_MR.smk --configfile $file
    done

elif [ $dataset == 'crosspheno' ]
then

    #This applies to running only the unidirectional MR of HCMR GWAS of plasma proteins with CMR phenotypes
    for file in ../input/configs/pp/*.yaml
    do
        echo "Carrying out unidirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile unidirectional_MR.smk --configfile $file 

    done

    #This applies to running only the unidirectional multivariate MR for CMR -> NTproBNP or TnT
    for file in ../input/configs/unidirectional/multivariate/*.yaml
    do
        echo "Carrying out unidirectional multivariate MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile multivariate_MR.smk --configfile $file 

    done

elif [ $dataset == 'ukb_exclHCM' ]
then

    #This applies to running the bidirectional MR with UKB GWAS of plasma proteins with HCM
    for file in ../input/configs/ukb_pp_noHarper_exclHCM/*.yaml
    do
        echo "Carrying out bidirectional  MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile bidirectional_MR.smk --configfile $file 
    done

elif [ $dataset == 'ukb_inclHCM' ]
then

    #This applies to running the bidirectional MR with UKB GWAS of plasma proteins with HCM
    for file in ../input/configs/ukb_pp_noHarper_inclHCM/*.yaml
    do
        echo "Carrying out bidirectional  MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile bidirectional_MR.smk --configfile $file 
    done

elif [ $dataset == 'ukb_inclHCM_multiphenoMR' ]
then

    #This applies to running the bidirectional MR with UKB GWAS of plasma proteins with HCM
    for file in ../input/configs/ukb_pp_noHarper_inclHCM/multiphenoMR_nonMTAG_HCMSNPs/*.yaml
    do
        echo "Carrying out bidirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile bidirectional_MR.smk --configfile $file 
    done

    for file in ../input/configs/ukb_pp_noHarper_inclHCM/multiphenoMR_nonMTAG_HCMSNPs/unidirectional/*.yaml
    do
        echo "Carrying out UNIdirectional MR on ${file}"
        snakemake --profile bmrc_profile_smk5  --snakefile unidirectional_MR.smk --configfile $file 
    done

else
    echo 'No dataset specified'

fi









