#!/bin/bash
#SBATCH --job-name=EnPACT_personalized_setup_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=01:00:00
#SBATCH --mem=2GB
#SBATCH --error=ERROR_LOG
#SBATCH --output=OUTPUT_LOG

source ~/.bashrc
conda activate /beagle3/haky/users/saideep/envs/enformer
module load R

# Path to EnPACT git repository with source code
path_to_enpact_repo="/beagle3/haky/users/saideep/github_repos/Con-EnPACT"

config_file=CONFIG_FILE

# Do predictDB standard

python $path_to_enpact_repo/scripts/2.1_predictDB_standard.py \
    $config_file False



