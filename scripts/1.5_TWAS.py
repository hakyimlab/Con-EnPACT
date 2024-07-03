import random
import time

import subprocess

import json

import os,sys

#################################################################
# 0.) Load required input variables, and create directories
#################################################################



print("Loading input variables")

json_file = sys.argv[1]
parameters = json.load(open(json_file))

project_dir = parameters["general_parameters"]["project_directory"]

intermediates_dir_generate_data = os.path.join(project_dir, "intermediates",
                       "generate_enpact_training_data")
intermediates_dir_train_enpact = os.path.join(project_dir, "intermediates",
                       "train_enpact")
intermediates_dir_linearization = os.path.join(project_dir, "intermediates",
                        "linearization")
intermediates_dir = os.path.join(project_dir, "intermediates",
                        "TWAS")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

twas_parameters = parameters["TWAS"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)


#################################################################
# 1.) Build submission scripts for TWAS analsis
#################################################################

print("Building submission scripts for TWAS analysis")

gwas_files = twas_parameters["GWAS_data"]

if twas_parameters["TWAS_method"] == "SPrediXcan":

    twas_template_path = os.path.join(script_directory, "run_SPrediXcan.sbatch")

    model_path = os.path.join(intermediates_dir_linearization,"predictDB","filtered_db","predict_db_Model_training_filtered.db")
    cov_path = os.path.join(intermediates_dir_linearization,"predictDB","filtered_db","predict_db_Model_training_filtered.txt.gz")

    for gwas in gwas_files.keys():
        twas_template = open(twas_template_path, "r").read()

        twas_template = twas_template.replace("JOBNAME", gwas)

        twas_template = twas_template.replace("SPREDIXCAN_PATH", twas_parameters["path_to_TWAS_software"])

        twas_template = twas_template.replace("MODEL_DB_PATH", model_path)
        twas_template = twas_template.replace("COVARIANCE_PATH", cov_path)

        twas_template = twas_template.replace("GWAS_FILE", gwas_files[gwas])

        os.makedirs(os.path.join(intermediates_dir, gwas), exist_ok=True)
        twas_template = twas_template.replace("OUTPUT_FILE", os.path.join(intermediates_dir, gwas, gwas+".txt"))

        with open(os.path.join(intermediates_dir, gwas, gwas+".sbatch"), "w") as o:
            o.write(twas_template)

        subprocess.run(["sbatch", os.path.join(intermediates_dir, gwas, gwas+".sbatch")],
                       cwd=os.path.join(intermediates_dir, gwas))


