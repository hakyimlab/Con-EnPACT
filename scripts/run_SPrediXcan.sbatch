#!/bin/bash
#SBATCH --job-name=SPrediXcan_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=26
#SBATCH --time=06:00:00
#SBATCH --mem=52GB
    
module load python
conda activate imlabtools

SPREDIXCAN_PATH \
--model_db_path MODEL_DB_PATH \
--covariance COVARIANCE_PATH \
--gwas_file GWAS_FILE \
--model_db_snp_key varID \
--keep_non_rsid \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--pvalue_column pvalue \
--verbosity 1 \
--throw \
--overwrite \
--output_file OUTPUT_FILE
