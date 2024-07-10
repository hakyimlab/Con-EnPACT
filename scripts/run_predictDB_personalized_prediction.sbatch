#!/bin/bash
#SBATCH --job-name=PredictDB_personalized_prediction_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --mem=42G


source ~/.bashrc
conda activate /project2/haky/festus/miniconda3/envs/imlab

python3 /beagle3/haky/users/charles/project/singleXcanDL/MetaXcan/software/Predict.py \
--model_db_path PATH_TO_DB \
--model_db_snp_key varID \
--text_genotypes TEXT_GENOTYPE_FILE \
--text_sample_ids TEXT_SAMPLE_FILE \
--prediction_output OUTPUT_FILE \
--prediction_summary_output OUTPUT_SUMMARY_FILE \
--verbosity 9 \
--throw