#!/bin/bash
#SBATCH --job-name=EnPACT_linearization_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=06:00:00
#SBATCH --mem=40GB
#SBATCH --error=ERROR_LOG
#SBATCH --output=OUTPUT_LOG

source ~/.bashrc
conda activate /project2/haky/festus/miniconda3/envs/predictdb_nf
module load R

# Path to EnPACT git repository with source code
path_to_enpact_repo="/beagle3/haky/users/saideep/github_repos/Con-EnPACT"

config_file=CONFIG_FILE

# Generate training data
python $path_to_enpact_repo/scripts/1.4_linearization.py \
    $config_file
