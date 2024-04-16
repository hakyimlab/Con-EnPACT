#!/bin/bash
#SBATCH --job-name=EnPACT_training_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=06:00:00
#SBATCH --mem=80GB
#SBATCH --error=ERROR_LOG
#SBATCH --output=OUTPUT_LOG

source ~/.bashrc
conda activate /beagle3/haky/users/saideep/envs/enformer

# Path to EnPACT git repository with source code
path_to_enpact_repo="/beagle3/haky/users/saideep/github_repos/Con-EnPACT"

config_file=CONFIG_FILE

# Generate training data
python $path_to_enpact_repo/scripts/1.1_generate_enpact_training_data.py \
    $config_file

module load R

# Train EnPACT
python $path_to_enpact_repo/scripts/1.2_train_EnPACT_elastic_net.py \
    $config_file

# Evaluate EnPACT
python $path_to_enpact_repo/scripts/1.3.1_evaluate_training.py \
    $config_file