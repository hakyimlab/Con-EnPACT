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

import enapct_utils
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
                        "linearization")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

linearization_parameters = parameters["linearization"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)

#################################################################
# 1.) Check predicted epigenome for each linearization dataset
#################################################################

print("Preparing features to linearize EnPACT model on")

linearization_datasets = linearization_parameters["linearization_datasets"]

# for ld in linearization_datasets.keys():


#TODO Add in checks here
#TODO Need to check that all the files exist, including mean summarized enformer predictions

#################################################################
# 2.) Make EnPACT predictions from personalized epigenome
#################################################################

print("Making EnPACT predictions")

window_size = parameters["generate_enpact_training_data"]["num_bins"]
context = parameters["general_parameters"]["context"]

model_path = os.path.join(intermediates_dir_train_enpact,
                        f"trained_enpact_eln_{context}.linear.rds")

for ld in linearization_datasets.keys():
    inds_file = linearization_datasets[ld]["individuals"]
    epigenome_pred_dir = linearization_datasets[ld]["epigenome_pred_dir"]
    with open(inds_file, "r") as inds_f:
        inds = inds_f.read().split()

    cur_lin_dir = os.path.join(intermediates_dir, ld)
    os.makedirs(os.path.join(cur_lin_dir, "enpact_preds"), exist_ok=True)

    pool = mp.Pool(processes=mp.cpu_count())

    arguments_for_starmap = []
    output_enpact_predictions = []
    for ind in inds:
        cur_enformer_pred = os.path.join(epigenome_pred_dir, ind+f"_{window_size}.txt")
        if not os.path.exists(cur_enformer_pred):
            print(cur_enformer_pred, " does not exist")
            continue
        cur_enpact_pred = os.path.join(cur_lin_dir, "enpact_preds", f"{ind}_{window_size}_enpact.txt")
        arguments_for_starmap.append((cur_enformer_pred, cur_enpact_pred, model_path))
        

#################################################################
# 3.) Format EnPACT predictions for PredictDB and prepare inputs
#################################################################

print("Formatting EnPACT predictions for PredictDB")

for ld in linearization_datasets.keys():

    cur_lin_dir = os.path.join(intermediates_dir, ld)

    paste_com = [
        "paste",
        "-d'\t'"
    ]

    inds_include = []
    for ind in inds_pred:
        if os.path.exists(os.path.join(cur_lin_dir, "enpact_preds",f"{ind}_{window_size}_enpact.txt")):
            paste_com.append(os.path.join(cur_lin_dir, "enpact_preds",f"{ind}_{window_size}_enpact.txt"))
            inds_include.append(ind)

    paste_com.append("|")
    paste_com.append("tail")
    paste_com.append("-n")
    paste_com.append("+2")

    paste_com.append(">")
    paste_com.append(os.path.join(cur_lin_dir, "enpact_preds", f"combined_predictions.txt"))

    print(" ".join(paste_com))

    subprocess.run(paste_com)


#################################################################
# 4.) Create PredictDB sbatch scripts
#################################################################

print("Creating PredictDB sbatch scripts")

for ld in linearization_datasets.keys():

    cur_lin_dir = os.path.join(intermediates_dir, ld)

    predictDB_dir = os.path.join(cur_lin_dir,"predictDB")
    os.makedirs(predictDB_dir, exist_ok=True)

    with open(os.path.join(predictDB_template), "r") as f:
        predictDB_script = f.read()

        predictDB_script = predictDB_script.replace("JOBNAME", "EnPACT_Linearization")

        predictDB_script = predictDB_script.replace("PREDICTDB_REPO", linearization_parameters["path_to_predictDB_pipeline"])
        predictDB_script = predictDB_script.replace("EXPRESSION", os.path.join(cur_lin_dir, "enpact_preds", f"combined_predictions.txt"))
        predictDB_script = predictDB_script.replace("OUTPUT_DIR", predictDB_dir)
        predictDB_script = predictDB_script.replace("GENE_ANNOTATION", )
        predictDB_script = predictDB_script.replace("SNP_ANNOTATION", linearization_datasets[ld]["snp_annotaton_file"])
        predictDB_script = predictDB_script.replace("GENOTYPE", linearization_datasets[ld]["genotype_file"])
        predictDB_script = predictDB_script.replace("NFOLDS", "10")


        with open(os.path.join(predictDB_dir,"run_predictDB.sbatch"), "w") as rp:
            rp.write(predictDB_script)
