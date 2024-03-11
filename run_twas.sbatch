#!/bin/bash
#SBATCH --job-name=TWAS_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --error=ERROR_LOG
#SBATCH --output=OUTPUT_LOG

module load python
conda activate imlabtools

path_to_enpact_repo="/beagle3/haky/users/saideep/github_repos/Con-EnPACT"

config_file=CONFIG_FILE

# Generate training data
python $path_to_enpact_repo/scripts/1.5_TWAS.py \
    $config_file