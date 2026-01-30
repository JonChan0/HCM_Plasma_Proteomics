#This script iterates over all the config files in ../input/configs to run the finemap.smk pipeline for each trait

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate gms

echo $(date)

for file in ../input/configs/hcmr/*.yaml #For HCMR version
do
    echo "Carrying out finemapping on ${file}"
    snakemake --profile bmrc_profile_smk5  --snakefile finemap_v2g.smk --configfile $file

done