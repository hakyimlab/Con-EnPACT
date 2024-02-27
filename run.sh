conda activate /beagle3/haky/users/saideep/envs/enformer

path_to_enpact_repo="/beagle3/haky/users/saideep/github_repos/Con-EnPACT"

# Generate training data
# python $path_to_enpact_repo/scripts/generate_enpact_training_data.py \
    # /beagle3/haky/users/saideep/projects/Con_EnPACT/models/flu_4bin_1milscaling/flu_4bin_1milscaling.json

module load R

# Train EnPACT
# python $path_to_enpact_repo/scripts/train_EnPACT_elastic_net.py \
    # /beagle3/haky/users/saideep/projects/Con_EnPACT/models/flu_4bin_1milscaling/flu_4bin_1milscaling.json

# Evaluate EnPACT
python $path_to_enpact_repo/scripts/evaluate_training.py \
    /beagle3/haky/users/saideep/projects/Con_EnPACT/models/flu_4bin_1milscaling/flu_4bin_1milscaling.json