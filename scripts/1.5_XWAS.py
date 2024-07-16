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
                        "XWAS")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

xwas_parameters = parameters["XWAS"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)


#################################################################
# 1.) Build submission scripts for XWAS analsis
#################################################################

print("Building submission scripts for XWAS analysis")

xwas_datasets = xwas_parameters["XWAS_datasets"]

xwas_template_path = os.path.join(script_directory, "run_SPrediXcan.sbatch")

for gwas in xwas_datasets.keys():

    if xwas_parameters["XWAS_method"] == "SPrediXcan":

        current_linearization_dataset = xwas_datasets[gwas]["linearization_dataset"]

        model_path = os.path.join(intermediates_dir_linearization,current_linearization_dataset,"predictDB","filtered_db","predict_db_Model_training_filtered.db")
        cov_path = os.path.join(intermediates_dir_linearization,current_linearization_dataset,"predictDB","filtered_db","predict_db_Model_training_filtered.txt.gz")

        xwas_template = open(xwas_template_path, "r").read()

        xwas_template = xwas_template.replace("JOBNAME", gwas)

        xwas_template = xwas_template.replace("SPREDIXCAN_PATH", xwas_parameters["path_to_XWAS_software"])

        xwas_template = xwas_template.replace("MODEL_DB_PATH", model_path)
        xwas_template = xwas_template.replace("COVARIANCE_PATH", cov_path)

        xwas_template = xwas_template.replace("GWAS_FILE", xwas_datasets[gwas]["GWAS_sum_stats"])

        os.makedirs(os.path.join(intermediates_dir, gwas), exist_ok=True)
        xwas_template = xwas_template.replace("OUTPUT_FILE", os.path.join(intermediates_dir, gwas, gwas+".txt"))

        with open(os.path.join(intermediates_dir, gwas, gwas+".sbatch"), "w") as o:
            o.write(xwas_template)

        subprocess.run(["sbatch", os.path.join(intermediates_dir, gwas, gwas+".sbatch")],
                    cwd=os.path.join(intermediates_dir, gwas))


