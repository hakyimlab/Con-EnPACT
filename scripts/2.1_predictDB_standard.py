import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import h5py
import time

import subprocess

import json

import kipoiseq

import multiprocessing as mp

import pyBigWig 
import os,sys

import h5py

import pickle as pkl

import epigenome_utils
import logging_utils

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
intermediates_dir = os.path.join(project_dir, "intermediates",
                        "predictDB_standard")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

predictdb_parameters = parameters["predictDB_standard"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)

run_sbatch = sys.argv[2]


#################################################################
# 1.) Build PredictDB submission script
#################################################################

print("Building PredictDB submission script")


with open(os.path.join(script_directory,"run_predictDB.sbatch"), "r") as f:
    predictDB_script = f.read()

    predictDB_script = predictDB_script.replace("JOBNAME", "PredictDB_standard")

    predictDB_script = predictDB_script.replace("PREDICTDB_REPO", predictdb_parameters["path_to_predictDB_pipeline"])
    predictDB_script = predictDB_script.replace("EXPRESSION", parameters["generate_enpact_training_data"]["input_files"]["normalized_expression_data"])
    predictDB_script = predictDB_script.replace("OUTPUT_DIR", os.path.join(intermediates_dir,"predictDB"))
    predictDB_script = predictDB_script.replace("GENE_ANNOTATION", predictdb_parameters["feature_annotation_file"])
    predictDB_script = predictDB_script.replace("SNP_ANNOTATION", predictdb_parameters["snp_annotation_file"])
    predictDB_script = predictDB_script.replace("GENOTYPE", predictdb_parameters["genotype_file"])
    predictDB_script = predictDB_script.replace("NFOLDS", str(predictdb_parameters["nfolds"]))

    os.makedirs(os.path.join(intermediates_dir, "predictDB"), exist_ok=True)
    with open(os.path.join(intermediates_dir, "run_predictDB_standard.sbatch"), "w") as rp:
        rp.write(predictDB_script)

    if sys.argv[2] == "True":
        subprocess.run(["sbatch", os.path.join(intermediates_dir, "run_predictDB_standard.sbatch")],
                       cwd= os.path.join(intermediates_dir, "predictDB"))
        

#################################################################
# 2.) Build PredictDB personalized prediction submission script
#################################################################

print("Building PredictDB personalized prediction submission script")

sample_file = os.path.join(intermediates_dir, "prediction_samples.txt")
with open(predictdb_parameters["genotype_file"], "r") as f:
    samples = f.readline().strip().split("\t")
    with open(sample_file, "w") as f:
        for sample in samples[1:]:
            f.write(sample + "\t" + sample + "\n")

# chromosome variant_id position allele1 allele2 MAF
predictdb_dosage_file = os.path.join(intermediates_dir, "prediction_dosage.txt")
with open(predictdb_parameters["genotype_file"], "r") as f:
    with open(predictdb_dosage_file, "w") as d:
        next(f)
        for line in f:

            dosages = line.strip().split("\t")[1:]
            variant_id = line.strip().split("\t")[0]
            chrom = variant_id.split("_")[0]
            pos = variant_id.split("_")[1]
            allele1 = variant_id.split("_")[2]
            allele2 = variant_id.split("_")[3]
            maf = 0.5

            write_list = [chrom, variant_id, pos, allele1, allele2, str(maf)]+dosages

            d.write("\t".join(write_list) + "\n")
            

with open(os.path.join(script_directory,"run_predictDB_personalized_prediction.sbatch"), "r") as f:
    predictDB_script = f.read()

    predictDB_script = predictDB_script.replace("JOBNAME", "PredictDB_personalized_predictions")

    predictDB_script = predictDB_script.replace("PATH_TO_DB", os.path.join(intermediates_dir,"predictDB","filtered_db","predict_db_Model_training_filtered.db"))
    predictDB_script = predictDB_script.replace("TEXT_GENOTYPE_FILE", predictdb_dosage_file)
    predictDB_script = predictDB_script.replace("TEXT_SAMPLE_FILE", sample_file)
    predictDB_script = predictDB_script.replace("OUTPUT_FILE", os.path.join(intermediates_dir,"personalized_predictions.txt"))
    predictDB_script = predictDB_script.replace("OUTPUT_SUMMARY_FILE", os.path.join(intermediates_dir,"personalized_predictions_summary.txt"))

    os.makedirs(os.path.join(intermediates_dir, "predictDB"), exist_ok=True)
    with open(os.path.join(intermediates_dir, "run_predictDB_personalized_predictions.sbatch"), "w") as rp:
        rp.write(predictDB_script)

    # if sys.argv[2] == "True":
    #     subprocess.run(["sbatch", os.path.join(intermediates_dir, "run_predictDB_personalized_prediction.sbatch")],
    #                    cwd= os.path.join(intermediates_dir, "run_predictDB_personalized_prediction"))